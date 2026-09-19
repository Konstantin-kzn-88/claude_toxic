"""Dry isothermal HEGADAS-based plume, HGSYSTEM chapter 7.A + 7.B.
Restricted implementation, not the complete HGSYSTEM/HEGADAS-T software.
See SECONDARY.md for equation mapping, boundary conditions and limitations.
"""
from dataclasses import dataclass,asdict
from math import pi,sqrt,gamma,isfinite
from types import SimpleNamespace
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from cloud import Gas,Vessel,Atmosphere,Options,validate,friction_velocity,R_GAS,G,_positive

A_EDGE = sqrt(pi)/2

@dataclass(frozen=True)
class GasFeed:
    rate_kg_s: float
    duration_s: float

@dataclass(frozen=True)
class PlumeOptions:
    distance_m: float=500.
    observation_s: float=300.
    space_step_m: float=.5
    max_step_m: float=1.
    rtol: float=1e-7
    atol: float=1e-9
    gaussian_core_fraction: float=1e-6

class SecondaryCloud:
    def __init__(self,gas,feed,met,options=PlumeOptions()):
        validate(gas,Vessel(1.,met.pressure_pa,met.temperature_k),met,Options())
        for obj in (feed,options):
            for k,v in asdict(obj).items():_positive(k,v)
        if not 1e-9<=options.gaussian_core_fraction<=1e-3:
            raise ValueError('Допуск исчезновения ядра: от 1e-9 до 1e-3')
        if met.temperature_k != met.surface_temperature_k:
            raise ValueError('Температуры газа, воздуха и поверхности должны совпадать')
        if options.distance_m>10000:raise ValueError('Расстояние не более 10000 м')
        if options.distance_m/options.space_step_m>50000:raise ValueError('Не более 50000 пространственных шагов')
        self.gas,self.feed,self.met,self.options=gas,feed,met,options
        self.q=feed.rate_kg_s;self.beta=1+met.alpha_wind
        self.hfactor=gamma(1/self.beta)/self.beta
        self.wind_k=met.wind_speed_m_s/gamma(1/self.beta)/(self.hfactor*met.wind_height_m)**met.alpha_wind
        self.rho_air=met.pressure_pa*met.air_molar_mass_kg_mol/(R_GAS*met.temperature_k)
        self.rho_gas=met.pressure_pa*gas.molar_mass_kg_mol/(R_GAS*met.temperature_k)
        if self.rho_gas<=self.rho_air:
            raise ValueError(
                f'{gas.name}: при температуре {met.temperature_k-273.15:g} °C '
                f'и давлении {met.pressure_pa/1000:g} кПа плотность газа '
                f'{self.rho_gas:.3f} кг/м³ не выше плотности воздуха {self.rho_air:.3f} кг/м³. '
                'Текущая вторичная модель рассчитывает только тяжёлый однофазный газ '
                'при температуре воздуха. Охлаждение источника, капли, фазовые переходы '
                'и подъём лёгкого газа пока не реализованы. '
                'Плотность вычисляется по молярной массе, температуре и давлению; '
                'плотность жидкости подставлять вместо неё нельзя.')
        self.u_star=friction_velocity(met)
        self.delta=dict(A=.22,B=.16,C=.11,D=.08,E=.06,F=.04)[met.stability]
        # Explicit inlet closure: pure-gas equivalent section B=H, Sy=0.
        # HEGADAS 2.6 without the Order-385 minimum-height correction.
        self.initial_b=brentq(lambda h:2*h*h*self.rho_gas*self.speed(h)-self.q,1e-10,1e6)
        if self.initial_b>met.alpha_height_limit_m:raise ValueError('Начальная высота выше диапазона alpha')

    def speed(self,height):return self.wind_k*np.asarray(height)**self.met.alpha_wind

    def sigma(self,x):
        x=np.asarray(x)
        return self.delta*x/np.sqrt(1+1e-4*x)

    def inverse_sigma(self,sigma):
        v=np.asarray(sigma)**2/self.delta**2
        return .5*(1e-4*v+np.sqrt((1e-4*v)**2+4*v))

    def diffusivity(self,width):
        # 7.A.2.3: k_y^e(W)=sigma * sigma_e'(x_e(sigma)), sigma=sqrt(2/pi)*W.
        sigma=sqrt(2/pi)*np.maximum(width,0.)
        xe=self.inverse_sigma(sigma)
        return sigma*self.delta*(1+.5e-4*xe)/(1+1e-4*xe)**1.5

    def state(self,x,y):
        total,b_eff,tau,sy2=y
        total=np.maximum(total,self.q)
        moles=self.q/self.gas.molar_mass_kg_mol+(total-self.q)/self.met.air_molar_mass_kg_mol
        density=total*self.met.pressure_pa/(moles*R_GAS*self.met.temperature_k)
        height=(total/(2*b_eff*density*self.wind_k))**(1/self.beta)
        speed=self.speed(height)
        sy=np.sqrt(np.maximum(sy2,0.))
        core=b_eff-A_EDGE*sy
        return dict(total_rate_kg_s=total,half_width_m=b_eff,core_half_width_m=core,sy_m=sy,
                    height_m=height,sz_m=height/self.hfactor,speed_m_s=speed,density_kg_m3=density,
                    travel_time_s=tau,centre_concentration_kg_m3=self.q/(2*b_eff*height*speed))

    def richardson(self,s):
        h=s['height_m'];rho=s['density_kg_m3']
        ri=G*h*np.maximum(1-self.rho_air/rho,0)/self.u_star**2
        rist=G*h*np.maximum(rho/self.rho_air-1,0)/self.u_star**2
        return ri,rist

    def collapse_margin(self,x,y):
        s=self.state(x,y);ri,rist=self.richardson(s)
        # 2.12b, multiplied through to avoid division by Ri at neutral buoyancy.
        return float(s['half_width_m']/s['height_m']-8/(3*.41)*np.sqrt(ri)*np.sqrt(1+.8*rist))

    def rhs(self,x,y,regime='gravity'):
        s=self.state(x,y);h=float(s['height_m']);rho=float(s['density_kg_m3']);u=float(s['speed_m_s'])
        _,rist=self.richardson(s)
        top=.41*self.u_star*self.beta/sqrt(1+.8*float(rist))
        b=y[1];sy=float(s['sy_m']);sy2dot=4*float(self.diffusivity(b))
        if regime=='gravity':
            growth=1.15/u*sqrt(max(G*h*(1-self.rho_air/rho),0))
        elif regime=='mixing':
            # 7.B (3): k_y*(Sy)=a^2 k_y^e(a Sy), a=sqrt(pi)/2.
            # This substitution follows from the prescribed Gaussian matching condition.
            growth=2*A_EDGE**2*float(self.diffusivity(A_EDGE*sy))/b
        else:
            growth=2*A_EDGE**2*float(self.diffusivity(b))/b
        # HEGADAS 2.9 uses standard molar volume V0=22.4 m3/kmol.
        # Air entrainment in mol/(m*s); no side entrainment during gravity spreading.
        molar_top=2*b*top/.0224
        moles=self.q/self.gas.molar_mass_kg_mol+(y[0]-self.q)/self.met.air_molar_mass_kg_mol
        molar_side=0. if regime=='gravity' else moles*growth/b  # 7.B 2.9*
        return [self.met.air_molar_mass_kg_mol*(molar_top+molar_side),growth,1/u,sy2dot]

    def run(self):
        self.segments=[];self.transitions=[]
        current=0.;y=np.array([self.q,self.initial_b,0.,0.])
        regime='gravity' if self.collapse_margin(0,y)<0 else 'mixing'
        if regime=='mixing':self.transitions.append(dict(type='gravity_collapse',x_m=0.))
        stop='requested_distance'
        while current<self.options.distance_m:
            def height(x,v):return self.met.alpha_height_limit_m-float(self.state(x,v)['height_m'])
            events=[height];names=['alpha_height_exceeded']
            if regime=='gravity':
                def transition(x,v):return -self.collapse_margin(x,v)
                events.append(transition);names.append('gravity_collapse')
            elif regime=='mixing':
                def transition(x,v):
                    return float(self.state(x,v)['core_half_width_m'])/v[1]-self.options.gaussian_core_fraction
                events.append(transition);names.append('gaussian_transition')
            for f in events:f.terminal=True;f.direction=-1
            sol=solve_ivp(lambda x,v:self.rhs(x,v,regime),(current,self.options.distance_m),y,
                dense_output=True,events=events,method='DOP853',max_step=self.options.max_step_m,
                rtol=self.options.rtol,atol=self.options.atol)
            if not sol.success:raise RuntimeError(sol.message)
            self.segments.append((current,float(sol.t[-1]),regime,sol.sol))
            stop=next((n for n,v in zip(names,sol.t_events) if len(v)),'requested_distance')
            if stop in ('requested_distance','alpha_height_exceeded'):break
            current=float(sol.t[-1]);y=sol.y[:,-1].copy()
            record=dict(type=stop,x_m=current)
            if stop=='gravity_collapse':regime='mixing'
            else:
                record['core_fraction_before_projection']=float(self.state(current,y)['core_half_width_m']/y[1])
                # Retain B, flux, height and centre concentration. Change Sy by <= tolerance.
                y[3]=(y[1]/A_EDGE)**2
                record['virtual_origin_m']=float(self.inverse_sigma(sqrt(y[3]/2))-current)
                regime='gaussian'
            self.transitions.append(record)
        end=self.segments[-1][1]
        self.solution=SimpleNamespace(sol=self._dense,t=np.array([0.,end]),success=True)
        self.x=np.unique(np.r_[np.arange(0,end,self.options.space_step_m),end,
                              [r['x_m'] for r in self.transitions]])
        self.profile=self.state(self.x,self._dense(self.x));self.stop_reason=stop
        warnings=['Ограниченная изотермическая реализация HEGADAS; независимая валидация не выполнена.',
                  'Сухой идеальный газ, однородная плотность воздуха; без термодинамики HF и фазовых переходов.',
                  'Входное сечение B=H, чистый газ. Струя и опорожнение ёмкости не рассчитываются.',
                  'Время осреднения боковой диффузии фиксировано: 600 с; это не длительность подачи.',
                  'Конечная подача моделируется переносом стационарного профиля, не полным HEGADAS-T; продольная диффузия отсутствует.']
        if stop!='requested_distance':warnings.append('Дальнейший участок не вычислен: '+stop)
        if np.min(self.profile['height_m'])<self.met.roughness_m:warnings.append('Эффективная высота ниже шероховатости: применимость профиля ветра ограничена.')
        return dict(status='experimental_hegadas_isothermal',model='HGSYSTEM HEGADAS 7.A + 7.B (restricted)',
                    stop_reason=stop,warnings=warnings,transitions=self.transitions,
                    averaging_time_s=600.,friction_velocity_m_s=self.u_star,
                    feed=asdict(self.feed),released_mass_kg=self.q*self.feed.duration_s,
                    initial_half_width_m=self.initial_b,initial_height_m=self.initial_b,
                    x_m=self.x.tolist(),profile={k:np.asarray(v).tolist() for k,v in self.profile.items()})

    def _dense(self,x):
        arr=np.asarray(x,float);flat=arr.ravel();out=np.empty((4,flat.size))
        if np.any((flat<0)|(flat>self.segments[-1][1])):raise ValueError('X вне рассчитанного диапазона')
        for i,(start,end,regime,fn) in enumerate(self.segments):
            mask=(flat>=start)&((flat<end) if i<len(self.segments)-1 else (flat<=end))
            if np.any(mask):out[:,mask]=fn(flat[mask])
        return out.reshape((4,)+arr.shape)

    def concentration(self,x,y=0.,z=0.,time_s=None):
        """time=None: maximum over [0, observation_s]. Beyond solved x is NaN, not zero."""
        x,y,z=np.broadcast_arrays(np.asarray(x,float),np.asarray(y,float),np.asarray(z,float))
        if not all(np.all(np.isfinite(a)) for a in [x,y,z]) or np.any(z<0):raise ValueError('Конечные координаты, Z >= 0')
        if time_s is not None and (not isfinite(time_s) or time_s<0):raise ValueError('Время >= 0')
        clipped=np.clip(x,0,self.x[-1])
        unique,inverse=np.unique(clipped.ravel(),return_inverse=True)
        profiles=self.state(unique,self.solution.sol(unique))
        s={k:np.asarray(v)[inverse].reshape(x.shape) for k,v in profiles.items()}
        sy=s['sy_m'];outside=np.maximum(abs(y)-np.maximum(s['core_half_width_m'],0),0)
        h=np.exp(-np.divide(outside**2,sy**2,out=np.full_like(outside,np.inf),where=sy>0))
        h=np.where((sy==0)&(outside==0),1.,h)
        c=s['centre_concentration_kg_m3']*h*np.exp(-(z/s['sz_m'])**self.beta)
        tau=s['travel_time_s']
        active=tau<=self.options.observation_s if time_s is None else (tau<=time_s)&(time_s<tau+self.feed.duration_s)
        return np.where(x>self.x[-1],np.nan,np.where((x>=0)&active,c,0.))

    def zone(self,threshold,height_m=0.):
        _positive('Порог',threshold)
        if not isfinite(height_m) or height_m<0:raise ValueError('Z >= 0')
        p=self.profile
        logratio=np.log(p['centre_concentration_kg_m3']/threshold)-(height_m/p['sz_m'])**self.beta
        active=(logratio>=0)&(p['travel_time_s']<=self.options.observation_s)
        half=np.where(active,np.maximum(p['core_half_width_m'],0)+p['sy_m']*np.sqrt(np.maximum(logratio,0)),0.)
        inds=np.flatnonzero(active)
        integrate=getattr(np,'trapezoid',None)
        if integrate is None:integrate=np.trapz
        return dict(threshold_kg_m3=threshold,section_height_m=height_m,
                    x_max_m=float(self.x[inds[-1]]) if len(inds) else None,
                    max_width_m=float(2*half.max()),area_m2=float(integrate(2*half,self.x)),
                    touches_solved_boundary=bool(active[-1]),spatial_step_m=self.options.space_step_m,
                    x_m=self.x.tolist(),half_width_m=half.tolist(),active=active.tolist())
