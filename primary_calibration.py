"""Default empirical correction fitted to one supplied chlorine TOXI-3 case.

This is a calibrated model, not a literal implementation of Order 385.
No accuracy claim is made for other substances or atmospheric/source conditions.
See validation/primary_chlorine/CALIBRATION.md for data and holdout checks.
"""
DRIFT_GAIN = 1.1882598077389053
LATERAL_GAIN = 1.1826850306787646
LATERAL_POWER = 0.9992350415915767
LATERAL_DENOMINATOR = 0.17876511213296598
LATERAL_LATE_POWER = -0.11924361595098754
REFERENCE_TIME_S = 600.0
INITIAL_SY_M = 1e-5
MAX_DURATION_S = 1800.0


def lateral_factor(t):
    """Dimensionless empirical multiplier on d(Sy²)/dt, t >= 0.

    The late branch already represents the fitted time dependence; a separate
    (t/600)**0.2 averaging multiplier must not also be applied to this correction.
    """
    early = min(max(t, 0.) / REFERENCE_TIME_S, 1.)
    late = max(t / REFERENCE_TIME_S, 1.)
    return (LATERAL_GAIN * early**LATERAL_POWER /
            (1 + LATERAL_DENOMINATOR * early) * late**LATERAL_LATE_POWER)


def metadata():
    return dict(id='chlorine_toxi3_20260926_v1', type='empirical_single_case',
                default_for_primary=True, independently_validated=False,
                drift_gain=DRIFT_GAIN, lateral_gain=LATERAL_GAIN,
                lateral_power=LATERAL_POWER, lateral_denominator=LATERAL_DENOMINATOR,
                lateral_late_power=LATERAL_LATE_POWER, reference_time_s=REFERENCE_TIME_S,
                initial_sy_m=INITIAL_SY_M, max_duration_s=MAX_DURATION_S,
                reference_conditions=dict(substance='chlorine', mass_kg=1000.,
                    pressure_abs_pa=110000., gas_temperature_k=313.15,
                    air_temperature_k=313.15, wind_m_s=1., wind_height_m=10.,
                    roughness_m=.55, stability='F', alpha_wind=.65),
                other_conditions_accuracy='not_established')
