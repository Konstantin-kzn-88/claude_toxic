"""Secondary GUI fields, SI JSON adapter."""
import copy
from input_data import FIELDS as PRIMARY_FIELDS,number
from secondary import SecondaryCloud,GasFeed,PlumeOptions
from cloud import Gas,Atmosphere
FIELDS=[f for f in PRIMARY_FIELDS if f[1] in ('gas','atmosphere') and f[2]!='surface_temperature_k']
FIELDS=[(tab,sec,key,'Общая температура газа, воздуха и поверхности, °C' if key=='temperature_k' else label,scale,offset,kind) for tab,sec,key,label,scale,offset,kind in FIELDS]
FIELDS += [('Вещество и ёмкость','feed','rate_kg_s','Расход газа на входе в шлейф, кг/с',1,0,'number'),
 ('Вещество и ёмкость','feed','duration_s','Длительность подачи газа, с',1,0,'number'),
 ('Расчёт и карты','plume_options','distance_m','Расчётное расстояние, м (до 10000)',1,0,'number'),
 ('Расчёт и карты','plume_options','observation_s','Время наблюдения, с',1,0,'number'),
 ('Расчёт и карты','plume_options','space_step_m','Шаг пространственной сетки, м',1,0,'number'),
 ('Расчёт и карты','plume_options','max_step_m','Шаг интегратора по X, м',1,0,'number'),
 ('Расчёт и карты',None,'section_height_m','Высота сечения Z, м',1,0,'number'),
 ('Дополнительно','plume_options','rtol','Относительная точность',1,0,'number'),
 ('Дополнительно','plume_options','gaussian_core_fraction','Допуск доли ядра для гауссовского перехода',1,0,'number'),
 ('Дополнительно','plume_options','atol','Абсолютная точность',1,0,'number')]

def completed(data):
    d=copy.deepcopy(data)
    for k,v in vars(PlumeOptions()).items():d.setdefault('plume_options',{}).setdefault(k,v)
    d.setdefault('threshold_labels',[]);d.setdefault('receptors',[]);d.setdefault('section_height_m',0.)
    d.setdefault('snapshot_times_s',[30.,60.,120.])
    return d

def display_values(data):
    d=completed(data);out={}
    for tab,sec,key,label,scale,offset,kind in FIELDS:
        v=(d[sec] if sec else d)[key]
        out[(sec,key)]=str(v) if kind in ('text','class') else format((v-offset)/scale,'.15g')
    return out

def collect(base,values,confirmed,thresholds_text,receptors_text,snapshots_text):
    d=completed(base)
    for tab,sec,key,label,scale,offset,kind in FIELDS:
        text=values[(sec,key)].strip()
        if kind in ('text','class'):
            if not text:raise ValueError(label+': заполните поле')
            value=text
        else:value=number(text,label)*scale+offset
        (d[sec] if sec else d)[key]=value
    d['gas']['phase_and_ideal_gas_assumptions_confirmed']=bool(confirmed)
    d['atmosphere']['surface_temperature_k']=d['atmosphere']['temperature_k']
    d['thresholds_kg_m3']=[];d['threshold_labels']=[]
    for line in thresholds_text.splitlines():
        if not line.strip():continue
        p=line.split(';')
        if len(p)!=2 or not p[0].strip():raise ValueError('Порог: Название; концентрация, кг/м³')
        c=number(p[1],'Порог')
        if c<=0:raise ValueError('Порог > 0')
        d['threshold_labels'].append(p[0].strip());d['thresholds_kg_m3'].append(c)
    if not d['thresholds_kg_m3']:raise ValueError('Нужен хотя бы один порог')
    d['receptors']=[]
    for line in receptors_text.splitlines():
        if not line.strip():continue
        p=line.split(';')
        if len(p)!=4 or not p[0].strip():raise ValueError('Точка: Название; X; Y; Z')
        x,y,z=[number(v,'Координата') for v in p[1:]]
        if z<0:raise ValueError('Z >= 0')
        d['receptors'].append(dict(name=p[0].strip(),x_m=x,y_m=y,z_m=z))
    times=[number(t,'Время снимка') for t in snapshots_text.split(';') if t.strip()]
    if len(times)>20 or any(t<0 or t>d['plume_options']['observation_s'] for t in times):raise ValueError('Время снимка вне интервала наблюдения (не более 20 снимков)')
    d['snapshot_times_s']=sorted(set(times))
    if d['section_height_m']<0:raise ValueError('Z >= 0')
    build_model(d)
    return d

def build_model(d):return SecondaryCloud(Gas(**d['gas']),GasFeed(**d['feed']),Atmosphere(**d['atmosphere']),PlumeOptions(**d['plume_options']))
