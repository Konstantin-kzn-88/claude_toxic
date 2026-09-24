"""Fuel mass within fixed mass-concentration limits, integrated in 3D.

Uses the existing ambient-equivalent kg/m3 limits, NOT local mole fractions
in nonisothermal or overlapping independent clouds. No explosion participation factor.
"""
import csv
import json
from functools import lru_cache
from math import gamma
from pathlib import Path
import numpy as np
from scipy.special import gammainc
from scipy.optimize import brentq

DEFAULTS = dict(time_step_s=5., nx=96, ny=64)
FIELDS = [
    ('Взрывоопасная масса','flammable_mass','time_step_s','Шаг поиска максимума массы, с',1,0,'number'),
    ('Взрывоопасная масса','flammable_mass','nx','Узлы интегрирования X (16–256)',1,0,'integer'),
    ('Взрывоопасная масса','flammable_mass','ny','Узлы интегрирования Y (16–128)',1,0,'integer'),
]


def defaults(d):
    for k,v in DEFAULTS.items():d.setdefault('flammable_mass',{}).setdefault(k,v)
    return d


def limits(d):
    """Only explicitly named limits are used; never infer missing UFL from a catalog."""
    found={}
    for name,value in zip(d.get('threshold_labels',[]),d.get('thresholds_kg_m3',[])):
        key=name.strip().upper()
        if key in ('НКПР','ВКПР'):
            if key in found:raise ValueError('Повторный порог '+key)
            found[key]=value
    if len(found)!=2:return None
    lo,hi=found['НКПР'],found['ВКПР']
    if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not np.isfinite(v) for v in (lo,hi)) or not 0<lo<hi:
        raise ValueError('Для массы требуется 0 < НКПР < ВКПР, кг/м³')
    return float(lo),float(hi)


def validate(d):
    cfg=defaults(d)['flammable_mass']
    dt=cfg['time_step_s']
    if isinstance(dt,bool) or not isinstance(dt,(float,int)) or not np.isfinite(dt) or dt<.1:
        raise ValueError('Шаг поиска массы должен быть >= 0.1 с')
    for k,upper in [('nx',256),('ny',128)]:
        v=cfg[k]
        if isinstance(v,bool) or not isinstance(v,(int,float)) or not np.isfinite(v) or int(v)!=v or not 16<=v<=upper:
            raise ValueError(f'{k}: целое число от 16 до {upper}')
    limits(d)


@lru_cache(None)
def quadrature(n):
    return np.polynomial.legendre.leggauss(int(n))


def column_mass(amplitudes, scales, beta, lower, upper):
    """Exact vertical integral between monotone-profile crossing heights (z >= 0)."""
    amplitudes,scales=np.broadcast_arrays(np.asarray(amplitudes,float),np.asarray(scales,float))
    total=amplitudes.sum(axis=0)
    if len(amplitudes)==1:
        a=amplitudes[0];sz=scales[0]
        u_low=np.maximum(np.log(np.maximum(a,lower)/lower),0.)
        u_high=np.maximum(np.log(np.maximum(a,upper)/upper),0.)
        return a*sz*gamma(1/beta)/beta*(gammainc(1/beta,u_low)-gammainc(1/beta,u_high))
    def height(level):
        # If every component <= level/n, their sum <= level.
        top=np.max(scales*np.maximum(np.log(np.maximum(amplitudes*len(amplitudes)/level,1.)),0.)**(1/beta),axis=0)
        bottom=np.zeros_like(top)
        for _ in range(42):
            mid=(bottom+top)/2
            above=np.sum(amplitudes*np.exp(-(mid/scales)**beta),axis=0)>level
            bottom=np.where(above,mid,bottom);top=np.where(above,top,mid)
        return np.where(total>level,(bottom+top)/2,0.)
    z0,z1=height(upper),height(lower)
    return np.sum(amplitudes*scales*gamma(1/beta)/beta*
        (gammainc(1/beta,(z1/scales)**beta)-gammainc(1/beta,(z0/scales)**beta)),axis=0)


class MassIntegrator:
    def __init__(self,model,lower,upper):
        self.p=getattr(model,'primary',None);self.s=getattr(model,'secondary',None)
        if self.p is None and self.s is None:
            if hasattr(model,'feed'):self.s=model
            else:self.p=model
        self.end=float(model.time_end if hasattr(model,'time_end') else
                       self.s.options.observation_s if self.s is not None else self.p.times[-1])
        self.beta=(self.p if self.p is not None else self.s).beta
        self.lower,self.upper=lower,upper

    def source_extent(self,t):
        s=self.s
        if s is None:return None,False
        tau_end=float(s.state(s.x[-1],s.solution.sol(s.x[-1]))['travel_time_s'])
        back=max(0.,t-s.feed.duration_s)
        def position(age):
            if age<=0:return 0.
            if age>=tau_end:return float(s.x[-1])
            return brentq(lambda x:float(s.state(x,s.solution.sol(x))['travel_time_s'])-age,0.,float(s.x[-1]))
        clipped=t>tau_end
        a,b=position(back),position(t)
        return ((a,b) if b>a else None),clipped

    def mass(self,t,nx=96,ny=64):
        if not 0<=t<=self.end:raise ValueError('Время вне общего интервала')
        p=self.p;s=self.s;floor=self.lower/(2 if p is not None and s is not None else 1)
        ps=p.state(p.solution.sol(t)) if p is not None else None
        extent,clipped=self.source_extent(t);segments=[]
        if ps is not None and ps['centre_concentration_kg_m3']>=floor:
            r=np.sqrt(ps['core_radius_m']**2+ps['sy2_m2']*np.log(ps['centre_concentration_kg_m3']/floor))
            a,b=ps['centre_x_m']-r,ps['centre_x_m']+r
            if s is not None and b>s.x[-1]:clipped=True;b=min(b,float(s.x[-1]))
            if b>a:segments.append((a,b))
        if extent is not None:segments.append(extent)
        if not segments:return 0.,clipped
        edges=sorted(set(v for pair in segments for v in pair))
        xn,xw=quadrature(nx);yn,yw=quadrature(ny);result=0.
        for left,right in zip(edges[:-1],edges[1:]):
            if not any(a<=(left+right)/2<=b for a,b in segments):continue
            x=(left+right)/2+(right-left)/2*xn
            y_max=np.zeros_like(x)
            if ps is not None and ps['centre_concentration_kg_m3']>=floor:
                r2=ps['core_radius_m']**2+ps['sy2_m2']*np.log(ps['centre_concentration_kg_m3']/floor)
                y_max=np.sqrt(np.maximum(r2-(x-ps['centre_x_m'])**2,0.))
            ss=None
            if extent is not None and extent[0]<=(left+right)/2<=extent[1]:
                ss=s.state(x,s.solution.sol(x));c=ss['centre_concentration_kg_m3']
                symax=np.maximum(ss['core_half_width_m'],0.)+ss['sy_m']*np.sqrt(np.maximum(np.log(c/floor),0.))
                y_max=np.maximum(y_max,np.where(c>=floor,symax,0.))
            y=y_max[:,None]*(yn+1)/2;xx=x[:,None];amps=[];scales=[]
            if ps is not None:
                extra=np.maximum((xx-ps['centre_x_m'])**2+y*y-ps['core_radius_m']**2,0.)
                horizontal=np.exp(-extra/ps['sy2_m2']) if ps['sy2_m2']>0 else (extra==0)
                amps.append(ps['centre_concentration_kg_m3']*horizontal);scales.append(np.full_like(y,ps['sz_m']))
            if ss is not None:
                sy=ss['sy_m'][:,None];extra=np.maximum(y-np.maximum(ss['core_half_width_m'][:,None],0.),0.)
                horizontal=np.exp(-np.divide(extra**2,sy**2,out=np.full_like(y,np.inf),where=sy>0))
                horizontal=np.where((sy==0)&(extra==0),1.,horizontal)
                amps.append(ss['centre_concentration_kg_m3'][:,None]*horizontal);scales.append(np.broadcast_to(ss['sz_m'][:,None],y.shape))
            if amps:
                columns=column_mass(amps,scales,self.beta,self.lower,self.upper)
                # y>=0 half plane, reflected by factor 2; z>=0 integrated analytically.
                result+=float(np.sum(xw*(right-left)/2*y_max*np.sum(columns*yw,axis=1)))
        return result,clipped


def save_flammable_mass(model,data,output):
    validate(data);out=Path(output);out.mkdir(parents=True,exist_ok=True);bounds=limits(data)
    if bounds is None:
        for name in ('flammable_mass.csv','flammable_mass.png'):(out/name).unlink(missing_ok=True)
        (out/'flammable_mass.json').write_text(json.dumps(dict(status='not_calculated',reason='Требуются оба порога с названиями НКПР и ВКПР, кг/м³'),ensure_ascii=False,indent=2),encoding='utf-8');return
    m=MassIntegrator(model,*bounds);cfg=data['flammable_mass'];step=cfg['time_step_s']
    if m.end/step>10000:raise ValueError('Не более 10000 временных отсчётов массы')
    times=list(np.arange(0.,m.end,step))+[m.end]
    if m.s is not None and m.s.feed.duration_s<=m.end:times.append(m.s.feed.duration_s)
    times=np.unique(times);rows=[]
    for t in times:
        coarse,clipped=m.mass(float(t),cfg['nx'],cfg['ny'])
        fine,clipped2=m.mass(float(t),2*cfg['nx'],2*cfg['ny'])
        rows.append(dict(time_s=float(t),fuel_mass_kg=fine,coarse_mass_kg=coarse,
                         spatial_difference_kg=abs(fine-coarse),domain_truncated=clipped or clipped2))
    peak=max(rows,key=lambda r:r['fuel_mass_kg'])
    error=max(r['spatial_difference_kg'] for r in rows)
    result=dict(status='experimental_ambient_equivalent_limits',quantity='fuel_mass_not_air_fuel_mixture',
        sampled_max_mass_kg=peak['fuel_mass_kg'],time_of_sampled_max_s=peak['time_s'],time_end_s=m.end,
        time_step_s=step,lower_kg_m3=bounds[0],upper_kg_m3=bounds[1],
        threshold_basis=data.get('thresholds_basis','Заданные массовые концентрации, кг/м³'),
        integration='Gauss quadrature X/Y, analytic vertical integral with numerical roots for sum',
        fine_nx=2*cfg['nx'],fine_ny=2*cfg['ny'],max_spatial_difference_kg=error,
        spatial_difference_relative_to_peak=error/peak['fuel_mass_kg'] if peak['fuel_mass_kg'] else 0.,
        domain_truncated=any(r['domain_truncated'] for r in rows),peak_at_time_boundary=peak['time_s'] in (0.,m.end),
        warnings=['Максимум по отсчётам времени; уменьшите шаг для проверки временной сходимости.',
                  'НКПР/ВКПР — фиксированные массовые пороги; для неизотермического облака это не проверка локальных объёмных долей.',
                  'При domain_truncated учтена только рассчитанная часть пространства; за ней масса неизвестна.',
                  'Коэффициент участия во взрыве не применён.'],history=rows)
    (out/'flammable_mass.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    with (out/'flammable_mass.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    from matplotlib import pyplot as plt
    fig,ax=plt.subplots(figsize=(9,5),constrained_layout=True)
    ax.plot(times,[r['fuel_mass_kg'] for r in rows],color='#126782')
    ax.scatter([peak['time_s']],[peak['fuel_mass_kg']],color='#dc2626',label=f"Максимум по отсчётам: {peak['fuel_mass_kg']:.4g} кг; {peak['time_s']:.4g} с")
    ax.set(xlabel='Время, с',ylabel='Масса горючего вещества, кг',title=f'Масса в диапазоне НКПР–ВКПР • 0–{m.end:g} с')
    ax.legend();ax.grid(alpha=.2);ax.set_ylim(bottom=0)
    note='Фиксированные пороги в кг/м³ • Без коэффициента участия во взрыве'
    if result['domain_truncated']:note+='\nОбласть ограничена: масса за её пределами неизвестна.'
    fig.supxlabel(note,fontsize=9)
    fig.savefig(out/'flammable_mass.png',dpi=180);plt.close(fig)
