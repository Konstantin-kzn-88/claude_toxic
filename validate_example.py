"""Reproducible diagnostic comparison, NOT a certification of the model.
Reference concentrations are approximate readings from figure 9-1.
"""
import json
from pathlib import Path
import numpy as np
from cloud import Gas, Vessel, Atmosphere, Options, PrimaryCloud
from zones import envelope


def main():
    root = Path(__file__).parent
    data = json.loads((root/'examples/chloromethane.json').read_text(encoding='utf-8'))
    x = np.array([100., 200., 400., 600., 1000.])
    reference = np.array([.20, .08, .022, .010, .0032])
    models = []
    for sample_step, solver_step in [(.1, .5), (.05, .25)]:
        m = PrimaryCloud(Gas(**data['gas']), Vessel(**data['vessel']), Atmosphere(**data['atmosphere']),
                         Options(**dict(data['options'], output_step_s=sample_step, max_step_s=solver_step)))
        m.run(); models.append(m)
    values = [m.maximum_concentration(x) for m in models]
    error = (values[0]/reference-1)*100
    report = dict(working_tolerance_percent=15., reference_quality='approximate_manual_reading_of_figure_9_1',
                  criterion_scope='relative_deviation_of_maximum_concentration_at_five_points',
                  concentration_comparison_pass=bool(np.all(abs(error)<=15)),
                  maximum_numerical_change_percent=float(np.max(abs(values[1]/values[0]-1))*100),
                  points=[dict(distance_m=float(xx), reference_approx_kg_m3=float(ref), computed_kg_m3=float(v),
                               relative_error_percent=float(err), within_15_percent=bool(abs(err)<=15))
                          for xx,ref,v,err in zip(x,reference,values[0],error)])
    report['zone_temporal_convergence'] = []
    for threshold in data['thresholds_kg_m3']:
        coarse, fine = [envelope(m,threshold) for m in models]
        report['zone_temporal_convergence'].append(dict(threshold_kg_m3=threshold,
            changes_percent={key:(fine[key]/coarse[key]-1)*100 for key in ['x_max_m','max_width_m','footprint_area_m2']}))
    (root/'results/validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__': main()
