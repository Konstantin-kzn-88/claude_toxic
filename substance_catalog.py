"""Reference properties transcribed from the supplied Order 385, Appendix 7."""
import copy

SOURCE = 'Приказ Ростехнадзора №385 от 02.11.2022, приложение 7, таблица 7-1'
# M [g/mol], gamma, normal boiling temperature [°C]. Not a phase envelope.
_ROWS = [
    ('Хлор', 70.9, 1.30, -34.1),
    ('Сероводород', 34.1, 1.30, -60.4),
    ('Пропан', 44., 1.13, -42.),
    ('Бутан', 58., 1.10, -.5),
    ('Аммиак', 17., 1.34, -33.4),
]
CATALOG = {name: dict(molar_mass_g_mol=m, adiabatic_index=g,
                     boiling_temperature_c=t, source=SOURCE)
           for name, m, g, t in _ROWS}
CATALOG['Хлорметан'] = dict(molar_mass_g_mol=51., adiabatic_index=1.25,
    boiling_temperature_c=None, source='Приказ №385, приложение 9, пример 1: молярная масса и показатель адиабаты')

TERRAIN_LABEL = 'Центры малых городов'

def city_f_defaults(data):
    """Fresh GUI defaults only: loading saved projects must preserve their weather."""
    d = copy.deepcopy(data)
    d['atmosphere'].update(stability='F', roughness_m=.55, alpha_wind=.655,
        alpha_height_limit_m=20., alpha_basis='Приказ №385, приложение 7: таблица 7-3 — центры малых городов, z0=0.55 м; '
        'таблица 7-5, F, высоты до 20 м: линейная интерполяция между z0=0.5 (α=0.65) и z0=0.6 (α=0.66); α=0.655. Интерполяция — решение реализации.')
    d['gas']['phase_and_ideal_gas_assumptions_confirmed'] = False
    return d
