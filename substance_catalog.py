"""Substance presets with reference or user-supplied property provenance."""
import copy
from cloud import R_GAS

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

# PCt50 and LCt50 [mg min/l]; LFL/UFL [% vol], Appendix 7 table 7-1.
for _name,_pc,_lc,_limits in [
    ('Хлор',.60,6.,None),('Сероводород',1.,15.,(4.3,45.)),
    ('Аммиак',15.,150.,(16.,25.)),('Пропан',0.,0.,(2.,9.5)),
    ('Бутан',0.,0.,(1.5,9.)),('Хлорметан',0.,0.,(8.1,17.4))]:
    CATALOG[_name].update(pct50_mg_min_l=_pc,lct50_mg_min_l=_lc,flammability_percent=_limits)
CATALOG['Хлорметан']['limits_source']='https://www.cdc.gov/niosh/npg/npgd0403.html'

_OIL_SOURCE = ('Пользовательские данные для паровой фазы нефти: M = 150 кг/кмоль; '
    'НКПР/ВКПР = 2,9/15 % об.; Cp ≈ 1,96 кДж/(кг·К) при 300 °C '
    '(температура интерпретирована как °C). k = Cp/(Cp − R/M), идеальный газ.')
CATALOG['Нефть — паровая фаза'] = dict(
    molar_mass_g_mol=150., adiabatic_index=1960./(1960.-R_GAS/.150),
    cp_j_kg_k=1960., cp_temperature_k=573.15,
    boiling_temperature_c=None, pct50_mg_min_l=0., lct50_mg_min_l=0.,
    flammability_percent=(2.9,15.), source=_OIL_SOURCE, limits_source=_OIL_SOURCE,
    input_note='Cp и рассчитанный k относятся к 300 °C. Применимость постоянного k '
        'при другой температуре требует обоснования. Задайте нижнюю границу однофазности '
        'для выбранных условий; 300 °C не является этой границей. '
        'Во вторичном облаке вводится расход уже испарившегося вещества, кг/с.')


def flammable_thresholds(name,molar_mass,temperature,pressure):
    import math
    from cloud import R_GAS
    if any(not math.isfinite(v) or v<=0 for v in (molar_mass,temperature,pressure)):
        raise ValueError('Для пересчёта НКПР/ВКПР нужны положительные M, T и P')
    limits=CATALOG[name]['flammability_percent']
    return [] if limits is None else [v/100*pressure*molar_mass/(R_GAS*temperature) for v in limits]


def refresh_collected_limits(d,base,text):
    name=d['gas']['name']
    if name==base.get('auto_limits_name') and text==base.get('auto_limits_text') and name in CATALOG:
        d['thresholds_kg_m3']=flammable_thresholds(name,d['gas']['molar_mass_kg_mol'],d['atmosphere']['temperature_k'],d['atmosphere']['pressure_pa'])
        d['threshold_labels']=['НКПР','ВКПР'] if d['thresholds_kg_m3'] else []
        d['thresholds_basis']='Объёмные пределы справочника, пересчёт при окружающих T/P; не локальные доли при переменной температуре. '+CATALOG[name].get('limits_source',SOURCE)
    return d
