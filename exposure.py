"""Linear inhalation dose: integral C dt; kg/m³·s -> mg·min/l (1000/60)."""
import json
import csv
from pathlib import Path
import numpy as np

TOXICITY_FIELDS = [
    ('Токсодоза','toxicity','pct50_mg_min_l','Пороговая PCt50, мг·мин/л (0 — не задана)',1,0,'number'),
    ('Токсодоза','toxicity','lct50_mg_min_l','Смертельная LCt50, мг·мин/л (0 — не задана)',1,0,'number'),
]

def defaults(d):
    for key in ('pct50_mg_min_l','lct50_mg_min_l'):d.setdefault('toxicity',{}).setdefault(key,0.)
    return d

def validate(d):
    defaults(d)
    for value in d['toxicity'].values():
        if isinstance(value,bool) or not isinstance(value,(float,int)) or not np.isfinite(value) or value<0:
            raise ValueError('Токсодоза: конечное число >= 0; 0 означает отсутствие критерия')
    a,b=(d['toxicity'][k] for k in ('pct50_mg_min_l','lct50_mg_min_l'))
    if a and b and b<a:raise ValueError('LCt50 должна быть не меньше PCt50')

def integrate_primary(model,x,y,z,end):
    times=np.unique(np.append(model.times[model.times<end],end))
    previous=model.concentration(0.,x,y,z);dose=np.zeros_like(previous,dtype=float)
    for t0,t1 in zip(times[:-1],times[1:]):
        current=model.concentration(float(t1),x,y,z)
        dose+=(previous+current)*.5*(t1-t0);previous=current
    return dose*1000/60

def integrate_secondary(model,x,y,z,end):
    x,y,z=np.broadcast_arrays(np.asarray(x,float),np.asarray(y,float),np.asarray(z,float))
    clipped=np.clip(x,0,model.x[-1]);unique,inverse=np.unique(clipped,return_inverse=True)
    arrival=np.asarray(model.state(unique,model.solution.sol(unique))['travel_time_s'])[inverse].reshape(x.shape)
    # Stationary amplitude and exact duration of the rectangular finite-feed pulse.
    amplitude=model.concentration(x,y,z,None)
    duration=np.maximum(0.,np.minimum(model.feed.duration_s,end-arrival))
    return amplitude*duration*1000/60

def dose(model,x,y=0.,z=0.,end=None):
    if hasattr(model,'primary'):
        end=model.time_end if end is None else min(end,model.time_end)
        return integrate_primary(model.primary,x,y,z,end)+integrate_secondary(model.secondary,x,y,z,end)
    if hasattr(model,'feed'):
        return integrate_secondary(model,x,y,z,model.options.observation_s if end is None else end)
    return integrate_primary(model,x,y,z,float(model.times[-1]) if end is None else end)

def save_exposure(model,data,output):
    validate(data)
    from matplotlib import pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.lines import Line2D
    out=Path(output)
    p=getattr(model,'primary',None);s=getattr(model,'secondary',None)
    if p is None and not hasattr(model,'feed') and not hasattr(model,'primary'):p=model
    if s is None and hasattr(model,'feed'):s=model
    end=model.time_end if hasattr(model,'primary') else s.options.observation_s if s is not None else float(p.times[-1])
    radius=max(st['core_radius_m']+6*np.sqrt(st['sy2_m2']) for st in p.states) if p is not None else 0.
    width=max(10.,radius,float(np.max(np.maximum(s.profile['core_half_width_m'],0)+6*s.profile['sy_m'])) if s is not None else 0.)
    xmax=float(s.x[-1]) if s is not None else max(st['centre_x_m'] for st in p.states)+radius
    x=np.linspace(-radius,xmax,301);y=np.linspace(-width,width,201);xx,yy=np.meshgrid(x,y)
    z=data.get('section_height_m',0.);field=dose(model,xx,yy,z)
    np.savez_compressed(out/'toxic_dose_fields.npz',x_m=x,y_m=y,z_m=z,dose_mg_min_l=field,time_end_s=end)
    criteria=[(label,data['toxicity'][key]) for label,key in [('PCt50','pct50_mg_min_l'),('LCt50','lct50_mg_min_l')] if data['toxicity'][key]>0]
    fig,ax=plt.subplots(figsize=(10,5.5),constrained_layout=True)
    peak=float(np.nanmax(field))
    floor=min(t for _,t in criteria)/10 if criteria else max(peak/1000,1e-12)
    visible=np.isfinite(field)&(field>=floor)&(field>0)
    cmap=plt.get_cmap('viridis').copy();cmap.set_bad('white')
    mesh=ax.pcolormesh(x,y,np.ma.masked_where(~visible,field),shading='auto',
        cmap=cmap,norm=LogNorm(vmin=floor,vmax=max(peak,floor*1.01)))
    if visible.any():
        left=min(0.,float(xx[visible].min()));right=max(0.,float(xx[visible].max()))
        half=max(float(np.max(np.abs(yy[visible]))),float(y[1]-y[0]))
        margin=.08*max(right-left,2*half,1.)
        ax.set_xlim(max(x[0],left-margin),min(x[-1],right+margin))
        ax.set_ylim(max(y[0],-half-margin),min(y[-1],half+margin))
    else:
        ax.text(.5,.5,f'Нет значений выше цветового порога {floor:.3g} мг·мин/л',
                transform=ax.transAxes,ha='center',fontsize=10)
    zones=[];handles=[]
    colors={'PCt50':'#f59e0b','LCt50':'#ef4444'}
    for label,threshold in criteria:
        mask=field>=threshold
        if np.nanmin(field)<threshold<np.nanmax(field):
            cs=ax.contour(x,y,field,levels=[threshold],colors=[colors[label]],linewidths=1.6)
            ax.clabel(cs,fmt={threshold:label},fontsize=9)
            handles.append(Line2D([0],[0],color=colors[label],label=f'{label} = {threshold:g} мг·мин/л'))
        touches=bool(mask[0].any() or mask[-1].any() or mask[:,0].any() or mask[:,-1].any())
        zones.append(dict(label=label,threshold_mg_min_l=threshold,reached=bool(mask.any()),
            x_min_m=float(xx[mask].min()) if mask.any() else None,x_max_m=float(xx[mask].max()) if mask.any() else None,
            touches_grid_boundary=touches))
    source=ax.scatter([0],[0],marker='x',color='black',label='Источник',zorder=5)
    ax.annotate('Ветер',xy=(.96,.93),xytext=(.78,.93),xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->'),va='center')
    ax.set(title=f'Накопленная токсодоза • 0–{end:.1f} с • Z={z:g} м',xlabel='X — по ветру, м',ylabel='Y — поперёк ветра, м')
    ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.15)
    ax.legend(handles=[source]+handles,loc='lower right',fontsize=8)
    fig.colorbar(mesh,ax=ax,shrink=.85,label='Токсодоза, мг·мин/л (логарифмическая шкала)')
    fig.supxlabel(f'Цветовой порог: {floor:.3g} мг·мин/л • Только рассчитанный интервал экспозиции',fontsize=9)
    fig.savefig(out/'toxic_dose_xy.png',dpi=180,bbox_inches='tight');plt.close(fig)
    observations=[]
    for point in data.get('receptors',[]):
        value=float(dose(model,*[point.get(k,0.) for k in ('x_m','y_m','z_m')]))
        observations.append(dict(point=point,dose_mg_min_l=value if np.isfinite(value) else None,
            status='computed_interval_only' if np.isfinite(value) else 'outside_solved_domain',
            criteria=[dict(label=label,threshold_mg_min_l=t,reached=bool(value>=t) if np.isfinite(value) else None) for label,t in criteria]))
    with (out/'toxic_dose_receptors.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.writer(f);writer.writerow(['name','x_m','y_m','z_m','dose_mg_min_l','time_end_s'])
        for row in observations:writer.writerow([row['point'].get('name',''),*[row['point'].get(k,0.) for k in ('x_m','y_m','z_m')],row['dose_mg_min_l'],end])
    result=dict(time_start_s=0.,time_end_s=end,units='mg min/l',definition='integral C dt; 1 kg s/m3 = 1000/60 mg min/l',
        method='primary: sampled trapezoids; secondary: exact rectangular pulse; combined: sum over common interval',
        zones=zones,receptors=observations,criteria=data['toxicity'],
        display_min_mg_min_l=floor,color_scale='log',
        warnings=['Только рассчитанный интервал: после остановки модели экспозиция не экстраполируется.',
                  'Сравнение с PCt50/LCt50, не расчёт вероятности по пробит-функции.',
                  'Границы зон приближённые по пространственной сетке; проверьте touches_grid_boundary.'])
    (out/'toxic_dose.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
