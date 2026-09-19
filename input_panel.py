"""Tkinter input editor; starts calculation in a separate process, no GUI freezing."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from datetime import datetime
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from cloud import Gas, Vessel, Atmosphere, Options, PrimaryCloud

from substance_catalog import CATALOG, TERRAIN_LABEL, city_f_defaults

ROOT=Path(__file__).resolve().parent

class InputPanel(ttk.Frame):
    def __init__(self,master,secondary=False):
        super().__init__(master)
        self.secondary=secondary
        if secondary: import secondary_input as inputs
        else: import input_data as inputs
        self.inputs=inputs
        self.example=ROOT/("examples/secondary_gas.json" if secondary else "examples/chloromethane.json")
        self.runner=ROOT/("run_secondary.py" if secondary else "run.py")
        self.process=None;self.log_handle=None;self.last_output=None;self.loading=False
        self.vars={};self.base={};self.current_file=None
        self.columnconfigure(0,weight=1);self.rowconfigure(1,weight=1)
        self.tabs=ttk.Notebook(self);self.tabs.grid(row=1,column=0,sticky='nsew',padx=12)
        for title in dict.fromkeys(f[0] for f in self.inputs.FIELDS):
            page=ttk.Frame(self.tabs,padding=12);self.tabs.add(page,text='Вещество и источник' if secondary and title=='Вещество и ёмкость' else title);page.columnconfigure(1,weight=1)
            for row,f in enumerate([f for f in self.inputs.FIELDS if f[0]==title]):
                _,sec,key,label,_,_,kind=f
                ttk.Label(page,text=label).grid(row=row,column=0,sticky='w',padx=(0,12),pady=7)
                var=tk.StringVar();self.vars[(sec,key)]=var
                w=ttk.Combobox(page,textvariable=var,values=list('ABCDEF'),state='readonly') if kind=='class' else ttk.Entry(page,textvariable=var)
                w.grid(row=row,column=1,sticky='ew',pady=7)
                if sec in ('gas','vessel','atmosphere','feed'):
                    var.trace_add('write',self.changed_physics)
            row+=1
            if title=='Вещество и ёмкость':
                self.catalog_choice=tk.StringVar(value='Выбрать из справочника')
                selector=ttk.Combobox(page,textvariable=self.catalog_choice,values=list(CATALOG),state='readonly')
                selector.grid(row=row,column=0,columnspan=2,sticky='ew')
                selector.bind('<<ComboboxSelected>>',self.apply_substance)
                row+=1
                self.confirmed=tk.BooleanVar()
                ttk.Checkbutton(page,text='Однофазность и применимость идеального газа для этих условий проверены',variable=self.confirmed).grid(row=row,column=0,columnspan=2,sticky='w',pady=15)
                ttk.Label(page,text='Расход задаётся после струевого участка, при атмосферном давлении.\nТемпература газа равна температуре воздуха и поверхности.' if secondary else 'Pабс = Pизб + Pатм. Свойства газа и пороги вводятся для выбранного вещества.\nСмена названия вещества сама по себе не меняет его свойства.',wraplength=850).grid(row=row+1,column=0,columnspan=2,sticky='w')
            if title=='Погода':
                ttk.Button(page,text='F · Центры малых городов',command=self.apply_city_weather).grid(row=row,column=0,columnspan=2,sticky='w')
                row+=1
                self.terrain_text=tk.StringVar()
                ttk.Label(page,textvariable=self.terrain_text).grid(row=row,column=0,columnspan=2,sticky='w')
                self.vars[('atmosphere','roughness_m')].trace_add('write',self.update_terrain)
                row+=1
                ttk.Label(page,text='При смене класса устойчивости или шероховатости проверьте α и его обоснование.\nАвтоматический подбор α из таблицы методики пока не реализован.',wraplength=850).grid(row=row,column=0,columnspan=2,sticky='w',pady=12)
            if title=='Расчёт и карты':
                self.snapshots=tk.StringVar()
                ttk.Label(page,text='Снимки облака, с (через ;)').grid(row=row,column=0,sticky='w',pady=7)
                ttk.Entry(page,textvariable=self.snapshots).grid(row=row,column=1,sticky='ew')
        p=ttk.Frame(self.tabs,padding=12);self.tabs.add(p,text='Пороги и точки');p.columnconfigure(0,weight=1);p.rowconfigure(1,weight=1);p.rowconfigure(3,weight=1)
        ttk.Label(p,text='Пороги: Название; концентрация в кг/м³ — один порог на строку').grid(row=0,column=0,sticky='w')
        self.thresholds=tk.Text(p,height=6,wrap='none');self.thresholds.grid(row=1,column=0,sticky='nsew',pady=8)
        ttk.Label(p,text='Контрольные точки: Название; X; Y; Z — координаты в метрах, X по ветру').grid(row=2,column=0,sticky='w')
        self.receptors=tk.Text(p,height=7,wrap='none');self.receptors.grid(row=3,column=0,sticky='nsew',pady=8)
        self.status=tk.StringVar()
        self.load(self.example, fresh=True)

    def apply_city_weather(self):
        defaults=city_f_defaults(self.base)['atmosphere']
        for key in ('stability','roughness_m','alpha_wind','alpha_height_limit_m','alpha_basis'):
            self.vars[('atmosphere',key)].set(str(defaults[key]))
        self.confirmed.set(False)

    def update_terrain(self,*_):
        try: city=abs(float(self.vars[('atmosphere','roughness_m')].get().replace(',', '.'))-.55)<1e-12
        except ValueError: city=False
        self.terrain_text.set(TERRAIN_LABEL+' (z₀ = 0,55 м)' if city else 'Шероховатость задана вручную')

    def apply_substance(self,*_):
        name=self.catalog_choice.get();item=CATALOG[name]
        self.vars[('gas','name')].set(name)
        self.vars[('gas','molar_mass_kg_mol')].set(str(item['molar_mass_g_mol']))
        self.vars[('gas','adiabatic_index')].set(str(item['adiabatic_index']))
        # The normal boiling point does not establish a valid phase envelope.
        self.vars[('gas','minimum_valid_temperature_k')].set('')
        self.confirmed.set(False)
        self.thresholds.delete('1.0','end')
        self.base.pop('thresholds_basis',None)
        note='Свойства загружены. Задайте нижнюю границу однофазной модели и пороги для выбранного вещества.'
        if item['boiling_temperature_c'] is not None:
            note+=f"\nСправочная температура кипения: {item['boiling_temperature_c']:g} °C; она не заменяет проверку фазового состояния при заданном давлении."
        if item['molar_mass_g_mol']<28.96:
            note+='\nПри температуре воздуха газ легче воздуха: текущая вторичная модель его не поддерживает.'
        messagebox.showinfo('Справочник веществ',note+'\n'+item['source'])

    def changed_physics(self,*_):
        if not self.loading and hasattr(self,'confirmed'):
            self.confirmed.set(False)

    def load(self,path,fresh=False):
        try:
            d=self.inputs.completed(json.loads(Path(path).read_text(encoding='utf-8-sig')))
            if fresh:d=city_f_defaults(d)
            self.catalog_choice.set('Выбрать из справочника')
            values=self.inputs.display_values(d) # Parse before touching the form.
            self.loading=True
            self.base=d;self.current_file=Path(path)
            for k,v in values.items(): self.vars[k].set(v)
            self.confirmed.set(d['gas'].get('phase_and_ideal_gas_assumptions_confirmed',False))
            self.snapshots.set('; '.join(str(t) for t in (d['snapshot_times_s'] if self.secondary else d['maps_xy']['snapshot_times_s'])))
            self.thresholds.delete('1.0','end')
            self.thresholds.insert('1.0','\n'.join(f"{d['threshold_labels'][i] if i<len(d['threshold_labels']) else 'Порог '+str(i+1)}; {c:.15g}" for i,c in enumerate(d['thresholds_kg_m3'])))
            self.receptors.delete('1.0','end')
            self.receptors.insert('1.0','\n'.join(f"{p.get('name','Точка')}; {p['x_m']}; {p.get('y_m',0)}; {p.get('z_m',0)}" for p in d['receptors']))
            self.status.set('Загружено: '+str(path))
        except Exception as e: messagebox.showerror('Не удалось открыть',str(e))
        finally: self.loading=False

    def data(self):
        return self.inputs.collect(self.base,{k:v.get() for k,v in self.vars.items()},self.confirmed.get(),self.thresholds.get('1.0','end'),self.receptors.get('1.0','end'),self.snapshots.get())

