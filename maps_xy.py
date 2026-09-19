"""Horizontal X-Y concentration maps. Array layout: concentration[y_index, x_index]."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from zones import envelope, section


def concentration_grid(model, x, y, height_m, time_s=None):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if x.ndim != 1 or y.ndim != 1 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError('Оси X и Y должны быть конечными одномерными массивами')
    if len(x) < 2 or len(y) < 2 or np.any(np.diff(x)<=0) or np.any(np.diff(y)<=0):
        raise ValueError('Каждая ось должна содержать минимум две возрастающие координаты')
    if not np.isfinite(height_m) or height_m < 0: raise ValueError('Высота Z >= 0')
    xx, yy = np.meshgrid(x, y, indexing='xy')
    if time_s is not None:
        return model.concentration(time_s, xx, yy, height_m)
    # Stream over time: no four-dimensional time-space array is allocated.
    result = np.zeros(xx.shape)
    for s in model.states:
        cz = s['centre_concentration_kg_m3'] * np.exp(-(height_m/s['sz_m'])**model.beta)
        extra = np.maximum((xx-s['centre_x_m'])**2 + yy*yy - s['core_radius_m']**2, 0.)
        horizontal = np.exp(-extra/s['sy2_m2']) if s['sy2_m2'] > 0 else (extra == 0)
        np.maximum(result, cz*horizontal, out=result)
    return result


def save_xy_maps(model, data, output, zones):
    output = Path(output)
    config = data.get('maps_xy', {})
    height = data.get('section_height_m', 0.)
    thresholds = data.get('thresholds_kg_m3', [])
    positive_peak = max(s['centre_concentration_kg_m3'] * np.exp(-(height/s['sz_m'])**model.beta) for s in model.states)
    floor = config.get('display_min_kg_m3', min(thresholds)/10 if thresholds else max(positive_peak/1000, 1e-12))
    if isinstance(floor, bool) or not np.isfinite(floor) or floor <= 0:
        raise ValueError('display_min_kg_m3 > 0')
    nx, ny = config.get('nx', 301), config.get('ny', 201)
    if any(isinstance(n, bool) or not isinstance(n, int) or n < 21 or n > 1001 for n in (nx, ny)):
        raise ValueError('nx, ny: целые от 21 до 1001')
    requested = config.get('snapshot_times_s', [10., 30., 60.])
    if not isinstance(requested, list) or len(requested) > 20:
        raise ValueError('Допускается не более 20 времён снимков')
    if any(isinstance(t, bool) or not isinstance(t, (int,float)) or not np.isfinite(t) or t < 0 for t in requested):
        raise ValueError('Время снимка должно быть конечным и неотрицательным')
    skipped = [t for t in requested if t > model.times[-1]]
    snapshots = sorted(set(t for t in requested if t <= model.times[-1]))
    # Fit the window to the union at the minimum displayed concentration.
    low_zone = envelope(model, floor, height)
    if low_zone['x_min_m'] is None:
        extent = [-20., 20., -20., 20.]
    else:
        xmin, xmax = min(0., low_zone['x_min_m']), max(0., low_zone['x_max_m'])
        half_width = max(1., low_zone['max_width_m']/2)
        margin = .08 * max(xmax-xmin, 2*half_width, 1.)
        extent = [xmin-margin, xmax+margin, -half_width-margin, half_width+margin]
    x = np.linspace(extent[0], extent[1], nx)
    y = np.linspace(extent[2], extent[3], ny)
    maximum = concentration_grid(model, x, y, height)
    vmax = max(positive_peak, floor*1.01)
    norm = LogNorm(vmin=floor, vmax=vmax)
    cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('#f6f7f9')
    labels = data.get('threshold_labels', [])
    colors = ['#ef4444', '#f59e0b', '#ec4899', '#22d3ee']
    meta = dict(axes={'X':'по ветру, м', 'Y':'поперёк ветра, м', 'Z_m':height},
                array_order='concentration[y_index, x_index]', extent_m=extent,
                nx=nx, ny=ny, display_min_kg_m3=floor, color_scale='logarithmic',
                color_max_kg_m3=vmax, time_end_s=float(model.times[-1]),
                maximum_time_step_s=model.options.output_step_s,
                status='experimental_unvalidated', skipped_snapshot_times_s=skipped,
                note='Цвет ниже нижнего предела скрыт; это не нулевая концентрация. Максимум — за вычисленный интервал.',
                images=[])
    arrays = dict(x_m=x, y_m=y, z_m=np.array(height), maximum_kg_m3=maximum,
                  time_end_s=np.array(model.times[-1]))

    def draw(values, filename, time=None):
        fig, ax = plt.subplots(figsize=(10, 6.5), constrained_layout=True)
        mesh = ax.pcolormesh(x, y, np.ma.masked_less(values, floor), cmap=cmap, norm=norm, shading='auto', rasterized=True)
        handles=[]
        for i, threshold in enumerate(thresholds):
            color=colors[i%len(colors)]
            label=labels[i] if i<len(labels) else f'{threshold:.4g} кг/м³'
            if time is None:
                zone=zones[i]
                if not zone['x_m']: continue
                zx=np.asarray(zone['x_m']); zy=np.asarray(zone['half_width_m'])
                ax.plot(zx,zy,c=color,lw=1.7); ax.plot(zx,-zy,c=color,lw=1.7)
            else:
                disk=section(model,time,threshold,height)
                if disk is None: continue
                theta=np.linspace(0,2*np.pi,361)
                ax.plot(disk['centre_x_m']+disk['radius_m']*np.cos(theta),disk['radius_m']*np.sin(theta),c=color,lw=1.7)
            handles.append(Line2D([0],[0],color=color,lw=1.7,label=label))
        ax.scatter([0],[0],c='#0f172a',marker='x',s=45,zorder=5)
        handles.append(Line2D([0],[0],marker='x',color='#0f172a',lw=0,label='Место выброса'))
        ax.annotate('Ветер',xy=(.95,.94),xytext=(.77,.94),xycoords='axes fraction',
                    arrowprops=dict(arrowstyle='->',color='#0f172a'),color='#0f172a',ha='left',va='center')
        title=f'Максимальная концентрация за 0–{model.times[-1]:g} с' if time is None else f'Облако в момент t = {time:g} с'
        ax.set(title=title+f' • Z = {height:g} м',xlabel='X — по направлению ветра, м',ylabel='Y — поперёк ветра, м',xlim=extent[:2],ylim=extent[2:])
        ax.set_aspect('equal',adjustable='box'); ax.grid(alpha=.15)
        ax.legend(handles=handles,loc='lower right',fontsize=8,framealpha=.95)
        bar=fig.colorbar(mesh,ax=ax,pad=.025,shrink=.86);bar.set_label('Концентрация, кг/м³ • логарифмическая шкала')
        fig.supxlabel('Экспериментальная модель. Контуры — пороги; серый фон — ниже цветовой шкалы.',fontsize=9)
        fig.savefig(output/filename,dpi=200);plt.close(fig)
        meta['images'].append(dict(file=filename,type='maximum_over_time' if time is None else 'snapshot',time_s=time))

    draw(maximum,'cloud_xy_max.png')
    for index,t in enumerate(snapshots):
        values=concentration_grid(model,x,y,height,t)
        arrays[f'snapshot_{index}_kg_m3']=values
        arrays[f'snapshot_{index}_time_s']=np.array(t)
        draw(values,f'cloud_xy_t{t:g}s.png',t)
    np.savez_compressed(output/'cloud_xy_fields.npz',**arrays)
    (output/'cloud_xy_metadata.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    return meta
