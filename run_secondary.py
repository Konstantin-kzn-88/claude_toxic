"""CLI for the simple isothermal secondary plume."""
import argparse,json,csv
from dataclasses import asdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from secondary import SecondaryCloud,GasFeed,PlumeOptions
from cloud import Gas,Atmosphere


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,default=Path(__file__).parent/'examples/secondary_gas.json')
    parser.add_argument('--output',type=Path,default=Path(__file__).parent/'results_secondary')
    args=parser.parse_args();d=json.loads(args.input.read_text(encoding='utf-8-sig'))
    m=SecondaryCloud(Gas(**d['gas']),GasFeed(**d['feed']),Atmosphere(**d['atmosphere']),PlumeOptions(**d['plume_options']))
    r=m.run();r['input']=d;args.output.mkdir(exist_ok=True,parents=True)
    height=d.get('section_height_m',0.);thresholds=d.get('thresholds_kg_m3',[])
    zones=[m.zone(t,height) for t in thresholds]
    (args.output/'result.json').write_text(json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    (args.output/'zones.json').write_text(json.dumps(dict(status=r['status'],warnings=r['warnings'],zones=zones),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    with (args.output/'profile.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.writer(f);writer.writerow(['x_m']+list(m.profile))
        writer.writerows([float(x)]+[float(m.profile[k][i]) for k in m.profile] for i,x in enumerate(m.x))
    fields=['threshold_kg_m3','x_max_m','max_width_m','area_m2','touches_solved_boundary']
    with (args.output/'zones.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(zones)
    floor=min(thresholds)/10 if thresholds else m.rho_gas/1000
    low=m.zone(floor,height)
    x_right=max(10.,(low['x_max_m'] or 0)*1.1)
    # Include all displayed snapshots above floor, even after shut-off.
    x_right=min(x_right,m.x[-1]);y_half=max(5.,low['max_width_m']*.65)
    x=np.linspace(0,x_right,401);y=np.linspace(-y_half,y_half,201)
    xx,yy=np.meshgrid(x,y,indexing='xy');labels=d.get('threshold_labels',[])
    snapshots=d.get('snapshot_times_s',[30.,60.,120.])
    if any(not isinstance(t,(int,float)) or not np.isfinite(t) or t<0 or t>m.options.observation_s for t in snapshots):
        raise ValueError('Время снимка должно лежать в интервале наблюдения')
    arrays=dict(x_m=x,y_m=y,z_m=np.array(height));manifest=[]
    for index,time in enumerate([None]+snapshots):
        values=m.concentration(xx,yy,height,time);arrays['max_kg_m3' if time is None else f'snapshot_{index}_kg_m3']=values
        fig,ax=plt.subplots(figsize=(10,5.5),constrained_layout=True)
        mesh=ax.pcolormesh(x,y,np.ma.masked_less(values,floor),shading='auto',cmap='viridis',norm=LogNorm(floor,max(m.rho_gas,floor*1.01)))
        for i,threshold in enumerate(thresholds):
            if np.nanmin(values)<threshold<np.nanmax(values):
                contours=ax.contour(x,y,values,levels=[threshold],colors=[['#ef4444','#f59e0b','#06b6d4'][i%3]],linewidths=1.5)
                ax.clabel(contours,fmt={threshold:labels[i] if i<len(labels) else f'{threshold:.4g}'},fontsize=8)
        if not np.any(np.isfinite(values)&(values>=floor)):
            ax.text(.5,.5,f'Нет значений выше цветового порога {floor:.3g} кг/м³',transform=ax.transAxes,ha='center',fontsize=10)
        ax.scatter([0],[0],marker='x',c='black',label='Вход в шлейф')
        ax.annotate('Ветер',xy=(.96,.92),xytext=(.77,.92),xycoords='axes fraction',arrowprops=dict(arrowstyle='->'),va='center')
        ax.set(xlabel='X — по ветру, м',ylabel='Y — поперёк ветра, м',
               title=f'Вторичное облако • '+(f'максимум за 0–{m.options.observation_s:g} с' if time is None else f't = {time:g} с')+f' • Z = {height:g} м')
        ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.15);ax.legend(loc='lower right',fontsize=8)
        fig.colorbar(mesh,ax=ax,shrink=.85,label='Концентрация, кг/м³ (логарифмическая шкала)')
        fig.supxlabel('HEGADAS 7.A + 7.B · изотермический шлейф · экспериментальная реализация.',fontsize=9)
        name='secondary_xy_max.png' if time is None else f'secondary_xy_t{time:g}s.png'
        fig.savefig(args.output/name,dpi=180);plt.close(fig);manifest.append(dict(file=name,time_s=time))
    np.savez_compressed(args.output/'secondary_fields.npz',**arrays)
    (args.output/'maps.json').write_text(json.dumps(dict(status=r['status'],array_order='[y_index,x_index]',time_end_s=m.options.observation_s,display_min_kg_m3=floor,images=manifest),ensure_ascii=False,indent=2),encoding='utf-8')
    fronts=[]
    for time in sorted(set([0.,m.feed.duration_s,m.options.observation_s]+snapshots)):
        if time>m.options.observation_s:continue
        age_back=max(time-m.feed.duration_s,0.)
        tau=m.profile['travel_time_s']
        fronts.append(dict(time_s=time,front_x_m=float(np.interp(time,tau,m.x)),back_x_m=float(np.interp(age_back,tau,m.x)),
                           front_beyond_domain=bool(time>tau[-1]),back_beyond_domain=bool(age_back>tau[-1])))
    (args.output/'fronts.json').write_text(json.dumps(fronts,ensure_ascii=False,indent=2),encoding='utf-8')
    observations=[]
    for point in d.get('receptors',[]):
        x0,y0,z0=[point.get(k,0.) for k in ['x_m','y_m','z_m']]
        if x0>m.x[-1]:
            observations.append(dict(point=point,status='outside_solved_domain'));continue
        if x0<0:
            observations.append(dict(point=point,status='upwind_zero_in_this_model'));continue
        s=m.state(x0,m.solution.sol(x0));arrival=float(s['travel_time_s']);departure=arrival+m.feed.duration_s
        peak=float(m.concentration(x0,y0,z0,None))
        observations.append(dict(point=point,arrival_s=arrival,last_exit_s=departure if departure<=m.options.observation_s else None,
                                 duration_in_window_s=max(0.,min(departure,m.options.observation_s)-arrival),
                                 max_concentration_kg_m3=peak,status=r['status']))
    (args.output/'receptors.json').write_text(json.dumps(observations,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({k:r[k] for k in ['status','stop_reason','released_mass_kg','initial_half_width_m','warnings']},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
