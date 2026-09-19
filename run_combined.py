"""Sum primary/secondary concentrations in a common time and coordinate system."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from combined import CombinedCloud


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding='utf-8-sig'))
    a, b = data['primary'], data['secondary']
    m = CombinedCloud.from_inputs(a, b)
    out = args.output; out.mkdir(parents=True, exist_ok=True)
    height = a.get('section_height_m', 0.)
    requested = sorted(set(a.get('maps_xy', {}).get('snapshot_times_s', []) + b.get('snapshot_times_s', [])))
    snapshots = [t for t in requested if 0 <= t <= m.time_end]
    times = np.unique(np.concatenate((m.times, snapshots)))
    thresholds = sorted(set(a.get('thresholds_kg_m3', []) + b.get('thresholds_kg_m3', [])))
    # Use the full solved downwind domain; undefined secondary values are not zeroed.
    radius = max(s['core_radius_m'] + 3*np.sqrt(s['sy2_m2']) for s in m.primary.states)
    width = max(10., radius, float(np.max(np.maximum(m.secondary.profile['core_half_width_m'], 0.) + 3*m.secondary.profile['sy_m'])))
    x = np.linspace(-radius, m.secondary.x[-1], 301)
    y = np.linspace(-width, width, 201)
    xx, yy = np.meshgrid(x, y)
    peak = m.maximum(xx, yy, height, times)
    arrays = dict(x_m=x, y_m=y, z_m=np.array(height), time_s=times, sampled_max_kg_m3=peak)
    images = []
    for i, t in enumerate([None] + snapshots):
        if t is None:
            values = peak; label = f'Максимум суммы по сетке времени 0–{m.time_end:g} с'; filename = 'combined_xy_max.png'
        else:
            c1, c2, values = m.concentration(t, xx, yy, height)
            arrays[f'primary_{i}_kg_m3'] = c1; arrays[f'secondary_{i}_kg_m3'] = c2; arrays[f'total_{i}_kg_m3'] = values
            label = f'Суммарная концентрация, t = {t:g} с'; filename = f'combined_xy_t{t:g}s.png'
        fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
        mesh = ax.pcolormesh(x, y, values, shading='auto', cmap='viridis')
        for threshold in thresholds:
            if np.nanmin(values) < threshold < np.nanmax(values):
                cs = ax.contour(x, y, values, levels=[threshold], colors='red')
                ax.clabel(cs, fmt='%g')
        ax.set(title=label + f' • Z = {height:g} м', xlabel='X по ветру, м', ylabel='Y, м')
        fig.colorbar(mesh, ax=ax, label='кг/м³')
        fig.supxlabel('Суперпозиция независимых облаков; общее начало координат и время старта.', fontsize=9)
        fig.savefig(out/filename, dpi=150); plt.close(fig)
        images.append(dict(file=filename, time_s=t, array_key='sampled_max_kg_m3' if t is None else f'total_{i}_kg_m3'))
    np.savez_compressed(out/'combined_fields.npz', **arrays)
    points = []
    for p in a.get('receptors', []) + b.get('receptors', []):
        if p not in points: points.append(p)
    observations = []
    with (out/'receptors.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['name', 'x_m', 'y_m', 'z_m', 'time_s', 'primary_kg_m3', 'secondary_kg_m3', 'total_kg_m3'])
        for index, p in enumerate(points):
            coords = [p.get(k, 0.) for k in ('x_m', 'y_m', 'z_m')]
            if coords[0] > m.secondary.x[-1]:
                observations.append(dict(point=p, status='outside_common_domain')); continue
            ts = m.receptor_times(coords[0], snapshots)
            values = np.array([m.concentration(float(t), *coords) for t in ts], dtype=float)
            writer.writerows([p.get('name', ''), *coords, t, *row] for t, row in zip(ts, values))
            fig, ax = plt.subplots(constrained_layout=True)
            for j, label in enumerate(['Первичное', 'Вторичное', 'Сумма']): ax.plot(ts, values[:, j], label=label)
            ax.set(title=p.get('name', '') or str(coords), xlabel='Время, с', ylabel='Концентрация, кг/м³'); ax.legend(); ax.grid(alpha=.2)
            fig.savefig(out/f'receptor_{index+1}.png', dpi=150); plt.close(fig)
            observations.append(dict(point=p, status='experimental_superposition', sampled_max_kg_m3=float(values[:, 2].max())))
    write(out/'receptors.json', observations)
    write(out/'result.json', dict(status='experimental_superposition', stop_reason='common_valid_interval',
        time_end_s=m.time_end, x_max_m=float(m.secondary.x[-1]), section_height_m=height,
        images=images, skipped_snapshot_times_s=[t for t in requested if t not in snapshots],
        max_time_step_s=float(np.max(np.diff(times))) if len(times)>1 else 0.,
        warnings=['Сумма независимых полей без взаимодействия облаков; одинаковое вещество, общий источник и старт t=0.',
                  'Максимум оценён по сетке времени, возможен пропуск пика между отсчётами.',
                  'Общий расчёт ограничен временем первичной модели и наблюдения вторичной; вне области сумма неизвестна.',
                  'Массы двух источников задаются пользователем без двойного учёта.'] ))

    from exposure import save_exposure
    exposure_data=dict(a,receptors=points)
    save_exposure(m,exposure_data,out)

if __name__ == '__main__': main()
