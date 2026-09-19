"""CLI: python run.py --input examples/chloromethane.json --output results"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cloud import Gas, Vessel, Atmosphere, Options, PrimaryCloud
from zones import envelope, receptor, validate_threshold
from maps_xy import save_xy_maps


def main():
    parser = argparse.ArgumentParser(description='Первичное облако: экспериментальная однофазная модель №385')
    parser.add_argument('--input', type=Path, default=Path(__file__).parent / 'examples/chloromethane.json')
    parser.add_argument('--output', type=Path, default=Path(__file__).parent / 'results')
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding='utf-8'))
    for threshold in data.get('thresholds_kg_m3', []): validate_threshold(threshold)
    model = PrimaryCloud(Gas(**data['gas']), Vessel(**data['vessel']), Atmosphere(**data['atmosphere']), Options(**data.get('options', {})))
    result = model.run()
    args.output.mkdir(parents=True, exist_ok=True)
    result['input'] = data
    (args.output / 'result.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    with (args.output / 'trajectory.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['time_s'] + list(result['states'][0]))
        writer.writeheader()
        writer.writerows(dict(time_s=t, **s) for t, s in zip(result['time_s'], result['states']))
    x = np.linspace(0, 1000, 1001)
    c = model.maximum_concentration(x)
    with (args.output / 'concentration.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['distance_m', 'max_concentration_kg_m3_in_simulated_interval'])
        writer.writerows(zip(x, c))
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    fig.suptitle('Первичное облако • экспериментальная модель, не валидирована', fontsize=12)
    ax[0].semilogy(x, np.maximum(c, 1e-10), color='#176b87', lw=2)
    ax[0].set(xlabel='Расстояние по ветру, м', ylabel='Максимальная концентрация, кг/м³', ylim=(1e-3, 3), title=f'Максимум за 0–{model.times[-1]:g} с')
    for threshold in data.get('thresholds_kg_m3', []):
        if not isinstance(threshold, (float, int)) or not np.isfinite(threshold) or threshold <= 0:
            raise ValueError('Пороги должны быть положительными конечными числами')
        ax[0].axhline(threshold, ls='--', lw=1, label=f'{threshold:.4g} кг/м³')
    if data.get('thresholds_kg_m3'): ax[0].legend(fontsize=8)
    ax[1].plot(model.times, [s['area_m2'] for s in model.states], color='#176b87', lw=2)
    ax[1].set(xlabel='Время, с', ylabel='Эффективная площадь, м²', title='Площадь πR² (не площадь опасной зоны)')
    for a in ax: a.grid(alpha=.25)
    fig.savefig(args.output / 'primary_cloud.png', dpi=160)
    plt.close(fig)
    zones = [envelope(model, threshold, data.get('section_height_m', 0.)) for threshold in data.get('thresholds_kg_m3', [])]
    zone_document = dict(model_status=result['status'], stop_reason=result['stop_reason'], warnings=result['warnings'], zones=zones)
    (args.output / 'zones.json').write_text(json.dumps(zone_document, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    fields = ['threshold_kg_m3', 'section_height_m', 'status', 'x_min_m', 'x_max_m', 'max_width_m', 'footprint_area_m2', 'time_end_s', 'active_at_end']
    with (args.output / 'zones.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader(); writer.writerows(zones)
    if zones:
        fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
        colors = ['#126782', '#d97706', '#577590', '#6a4c93']
        for i, zone in enumerate(zones):
            if not zone['x_m']: continue
            zx = np.array(zone['x_m']); zy = np.array(zone['half_width_m'])
            ax.fill_between(zx, -zy, zy, color=colors[i % len(colors)], alpha=.15)
            ax.plot(zx, zy, color=colors[i % len(colors)], label=f"{zone['threshold_kg_m3']:.4g} кг/м³")
            ax.plot(zx, -zy, color=colors[i % len(colors)])
        ax.scatter([0], [0], c='black', marker='x', label='Место выброса')
        ax.set(xlabel='X — по направлению ветра, м', ylabel='Y — поперёк ветра, м', title=f'Зоны достижения концентраций за 0–{model.times[-1]:g} с\nЭкспериментальный расчёт')
        ax.set_aspect('equal', adjustable='box'); ax.grid(alpha=.25); ax.legend(fontsize=8)
        fig.savefig(args.output / 'zones.png', dpi=160); plt.close(fig)
    save_xy_maps(model, data, args.output, zones)
    observations = []
    for point in data.get('receptors', []):
        coords = {k: point[k] for k in ['x_m', 'y_m', 'z_m'] if k in point}
        observation = receptor(model, **coords, thresholds=data.get('thresholds_kg_m3', []))
        observation['name'] = point.get('name', '')
        observations.append(observation)
    (args.output / 'receptors.json').write_text(json.dumps(dict(model_status=result['status'], stop_reason=result['stop_reason'], warnings=result['warnings'], receptors=observations), ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    if observations:
        with (args.output / 'receptors.csv').open('w', encoding='utf-8-sig', newline='') as f:
            writer = csv.writer(f); writer.writerow(['name', 'x_m', 'y_m', 'z_m', 'time_s', 'concentration_kg_m3'])
            for o in observations:
                writer.writerows([o['name'], o['x_m'], o['y_m'], o['z_m'], t, c] for t, c in zip(o['time_s'], o['concentration_kg_m3']))
        fig, ax = plt.subplots(figsize=(10, 4.5), constrained_layout=True)
        for o in observations: ax.plot(o['time_s'], o['concentration_kg_m3'], label=o['name'] or f"x={o['x_m']} м")
        for threshold in data.get('thresholds_kg_m3', []): ax.axhline(threshold, ls='--', lw=.8, c='grey')
        ax.set(xlabel='Время, с', ylabel='Концентрация, кг/м³', title='Концентрации в контрольных точках • экспериментальный расчёт')
        ax.legend(); ax.grid(alpha=.25)
        fig.savefig(args.output / 'receptors.png', dpi=160); plt.close(fig)
    print(json.dumps(dict(initial=result['initial'], stop_reason=result['stop_reason'], warnings=result['warnings']), ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
