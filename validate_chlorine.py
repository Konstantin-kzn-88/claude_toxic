"""Reproduce the supplied chlorine comparison without fitting physical equations.

python validate_chlorine.py --output results_chlorine_validation
All 68 reference rows are retained; rows after model termination are marked.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from cloud import Gas, Vessel, Atmosphere, Options, PrimaryCloud
from primary_calibration import metadata

ROOT = Path(__file__).resolve().parent
CASE = ROOT / 'validation' / 'primary_chlorine'
FIELDS = ('centre_x_m', 'speed_m_s', 'total_mass_kg', 'radius_m',
          'core_radius_m', 'height_m', 'sy_m', 'sz_m', 'density_kg_m3',
          'temperature_k', 'centre_concentration_kg_m3')


def compare(output):
    data = json.loads((CASE / 'input.json').read_text(encoding='utf-8'))
    with (CASE / 'reference.csv').open(newline='') as stream:
        reference = [{k: float(v) for k, v in row.items()}
                     for row in csv.DictReader(stream)]
    model = PrimaryCloud(Gas(**data['gas']), Vessel(**data['vessel']),
                         Atmosphere(**data['atmosphere']), Options(**data['options']))
    result = model.run()
    end = float(model.times[-1])
    rows = []
    fit = json.loads((CASE / 'fit_parameters.json').read_text(encoding='utf-8'))
    for row in reference:
        available = row['time_s'] <= end
        state = model.state(model.solution.sol(row['time_s'])) if available else None
        if available:
            state['sy_m'] = float(np.sqrt(state['sy2_m2']))
        for key in FIELDS:
            expected = row[key]
            actual = float(state[key]) if available else None
            group = ('fit' if int(row['row']) in fit['train_rows'] else
                     'control' if int(row['row']) in fit['holdout_rows'] else 'startup')
            rows.append(dict(row=int(row['row']), time_s=row['time_s'], group=group, parameter=key,
                             reference=expected, model=actual,
                             difference=actual-expected if available else None,
                             percent=100*(actual/expected-1) if available and expected else None,
                             status='compared' if available else 'outside_model_interval'))
    output.mkdir(parents=True, exist_ok=True)
    with (output / 'comparison.csv').open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    first, second = reference[:2]
    # Independent consistency check on supplied coordinates, not on our solver.
    initial_secant = ((second['centre_x_m']-first['centre_x_m']) /
                      (second['time_s']-first['time_s']))
    summary = dict(
        status='comparison_only_not_validated', reference_rows=len(reference),
        compared_rows=sum(r['time_s'] <= end for r in reference),
        excluded_rows=sum(r['time_s'] > end for r in reference),
        model_end_s=end, stop_reason=result['stop_reason'],
        calibration=metadata(),
        initial_reference_dx_dt_m_s=initial_secant,
        initial_reference_effective_speed_m_s=first['speed_m_s'],
        initial_reference_speed_ratio=initial_secant/first['speed_m_s'],
        reference_cp_j_kg_k=480., reference_cv_j_kg_k=362.736,
        model_cp_j_kg_k=model.cp_g, model_cv_j_kg_k=model.cv_g,
        max_absolute_percent={k: max(abs(r['percent']) for r in rows
                                   if r['parameter'] == k and r['percent'] is not None)
                              for k in FIELDS},
        control_max_absolute_percent={k: max(abs(r['percent']) for r in rows
                                   if r['group'] == 'control' and r['parameter'] == k
                                   and r['percent'] is not None) for k in FIELDS},
        caveats=['Reference program/version and certificate not identified.',
                 f"{sum(r['time_s'] <= end for r in reference)}/{len(reference)} times overlap the model interval.",
                 'Pressure is interpreted as absolute; surface temperature is assumed 313.15 K.',
                 'Air constants use model defaults, not fully specified reference inputs.',
                 'Alpha is held at 0.65 up to 40 m explicitly for the empirical reference case.',
                 'Use validate_chlorine_zones.py for the 1800 s dose-zone comparison.',
                 'Control rows did not enter numerical fitting, but belong to the same single accident.'])
    (output / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False,
                                                   indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results_chlorine_validation')
    args = parser.parse_args()
    print(json.dumps(compare(args.output), ensure_ascii=False, indent=2))
