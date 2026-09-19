"""User input schema and reversible conversion between displayed units and SI JSON."""
import copy
import math
from cloud import Gas, Vessel, Atmosphere, Options, PrimaryCloud

# tab, section, key, label, scale (SI/display), offset (SI), kind
FIELDS = [
 ('Вещество и ёмкость','gas','name','Вещество',1,0,'text'),
 ('Вещество и ёмкость','gas','molar_mass_kg_mol','Молярная масса, г/моль',.001,0,'number'),
 ('Вещество и ёмкость','gas','adiabatic_index','Показатель адиабаты',1,0,'number'),
 ('Вещество и ёмкость','gas','minimum_valid_temperature_k','Нижняя граница однофазной модели, °C',1,273.15,'number'),
 ('Вещество и ёмкость','vessel','volume_m3','Объём газовой полости, м³',1,0,'number'),
 ('Вещество и ёмкость','vessel','pressure_abs_pa','АБСОЛЮТНОЕ давление газа, МПа',1e6,0,'number'),
 ('Вещество и ёмкость','vessel','temperature_k','Температура газа, °C',1,273.15,'number'),
 ('Погода','atmosphere','temperature_k','Температура воздуха, °C',1,273.15,'number'),
 ('Погода','atmosphere','surface_temperature_k','Температура поверхности, °C',1,273.15,'number'),
 ('Погода','atmosphere','pressure_pa','Атмосферное давление, кПа',1e3,0,'number'),
 ('Погода','atmosphere','wind_speed_m_s','Скорость ветра, м/с',1,0,'number'),
 ('Погода','atmosphere','wind_height_m','Высота измерения ветра, м',1,0,'number'),
 ('Погода','atmosphere','roughness_m','Шероховатость поверхности, м',1,0,'number'),
 ('Погода','atmosphere','stability','Устойчивость атмосферы',1,0,'class'),
 ('Погода','atmosphere','alpha_wind','Коэффициент профиля ветра α',1,0,'number'),
 ('Погода','atmosphere','alpha_basis','Обоснование коэффициента α',1,0,'text'),
 ('Погода','atmosphere','alpha_height_limit_m','Предельная высота для выбранного α, м',1,0,'number'),
 ('Расчёт и карты','options','duration_s','Продолжительность, с (не более 600)',1,0,'number'),
 ('Расчёт и карты','options','output_step_s','Шаг выходных результатов, с',1,0,'number'),
 ('Расчёт и карты','options','max_step_s','Максимальный шаг интегратора, с',1,0,'number'),
 ('Расчёт и карты',None,'section_height_m','Высота сечения Z, м',1,0,'number'),
 ('Расчёт и карты','maps_xy','nx','Сетка X, число узлов (21–1001)',1,0,'integer'),
 ('Расчёт и карты','maps_xy','ny','Сетка Y, число узлов (21–1001)',1,0,'integer'),
 ('Дополнительно','options','rtol','Относительная точность интегратора',1,0,'number'),
 ('Дополнительно','options','atol','Абсолютная точность интегратора',1,0,'number'),
 ('Дополнительно','atmosphere','air_molar_mass_kg_mol','Молярная масса воздуха, г/моль',.001,0,'number'),
 ('Дополнительно','atmosphere','air_cp_j_kg_k','Теплоёмкость воздуха Cp, Дж/(кг·K)',1,0,'number'),
]


def number(text, label):
    try: value=float(str(text).strip().replace(',','.'))
    except (ValueError,TypeError): raise ValueError(f'{label}: введите число') from None
    if not math.isfinite(value): raise ValueError(f'{label}: требуется конечное число')
    return value


def completed(data):
    d=copy.deepcopy(data)
    for key, defaults in [('options',vars(Options())),('maps_xy',{'nx':301,'ny':201,'snapshot_times_s':[10.,30.,60.]})]:
        d.setdefault(key,{})
        for name,value in defaults.items(): d[key].setdefault(name,value)
    for name,value in [('alpha_height_limit_m',20.),('air_molar_mass_kg_mol',.02896),('air_cp_j_kg_k',1005.)]:
        d['atmosphere'].setdefault(name,value)
    d.setdefault('section_height_m',0.)
    d.setdefault('thresholds_kg_m3',[]);d.setdefault('threshold_labels',[]);d.setdefault('receptors',[])
    return d


def display_values(data):
    d=completed(data);out={}
    for tab,sec,key,label,scale,offset,kind in FIELDS:
        value=(d[sec] if sec else d)[key]
        out[(sec,key)]=str(value) if kind in ('text','class') else format((value-offset)/scale,'.15g')
    return out


def collect(base, values, confirmed, thresholds_text, receptors_text, snapshots_text):
    d=completed(base)
    for tab,sec,key,label,scale,offset,kind in FIELDS:
        raw=values[(sec,key)].strip()
        if kind in ('text','class'):
            if not raw: raise ValueError(f'{label}: заполните поле')
            value=raw
        else:
            value=number(raw,label)*scale+offset
            if kind=='integer':
                if value != int(value): raise ValueError(f'{label}: требуется целое число')
                value=int(value)
        (d[sec] if sec else d)[key]=value
    d['gas']['phase_and_ideal_gas_assumptions_confirmed']=bool(confirmed)
    thresholds=[];labels=[]
    for line in thresholds_text.splitlines():
        if not line.strip(): continue
        parts=line.split(';')
        if len(parts)!=2: raise ValueError('Пороги: каждая строка должна иметь вид Название; концентрация, кг/м³')
        name=parts[0].strip();value=number(parts[1],'Порог')
        if not name or value<=0: raise ValueError('Нужны название порога и концентрация > 0')
        labels.append(name);thresholds.append(value)
    if not thresholds: raise ValueError('Добавьте хотя бы один порог концентрации')
    d['thresholds_kg_m3']=thresholds;d['threshold_labels']=labels
    d['thresholds_basis']='Заданы пользователем в кг/м³ через окно исходных данных.'
    points=[]
    for line in receptors_text.splitlines():
        if not line.strip(): continue
        parts=line.split(';')
        if len(parts)!=4: raise ValueError('Точки: каждая строка должна иметь вид Название; X; Y; Z (м)')
        if not parts[0].strip(): raise ValueError('Укажите название контрольной точки')
        xyz=[number(t,'Координата точки') for t in parts[1:]]
        if xyz[2]<0: raise ValueError('Высота контрольной точки Z >= 0')
        points.append(dict(name=parts[0].strip(),x_m=xyz[0],y_m=xyz[1],z_m=xyz[2]))
    d['receptors']=points
    times=[number(t,'Время снимка') for t in snapshots_text.split(';') if t.strip()]
    if len(times)>20 or any(t<0 or t>d['options']['duration_s'] for t in times):
        raise ValueError('Не более 20 снимков; время каждого — от 0 до продолжительности расчёта')
    d['maps_xy']['snapshot_times_s']=sorted(set(times))
    if any(not 21<=d['maps_xy'][k]<=1001 for k in ['nx','ny']): raise ValueError('Сетка X и Y: от 21 до 1001 узла')
    if d['section_height_m']<0: raise ValueError('Высота сечения Z >= 0')
    if 'display_min_kg_m3' in d['maps_xy']:
        val=d['maps_xy']['display_min_kg_m3']
        if isinstance(val,bool) or not isinstance(val,(int,float)) or not math.isfinite(val) or val<=0:
            raise ValueError('Некорректный нижний предел цветовой шкалы в JSON')
    PrimaryCloud(Gas(**d['gas']),Vessel(**d['vessel']),Atmosphere(**d['atmosphere']),Options(**d['options']))
    return d
