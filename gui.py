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

ROOT=Path(__file__).resolve().parent

class App(tk.Tk):
    def __init__(self,secondary=False):
        super().__init__()
        self.secondary=secondary
        if secondary: import secondary_input as inputs
        else: import input_data as inputs
        self.inputs=inputs
        self.example=ROOT/("examples/secondary_gas.json" if secondary else "examples/chloromethane.json")
        self.runner=ROOT/("run_secondary.py" if secondary else "run.py")
        self.title(('Вторичное' if secondary else 'Первичное')+' облако — исходные данные');self.geometry('1080x780');self.minsize(880,620)
        self.process=None;self.log_handle=None;self.last_output=None;self.loading=False
        self.vars={};self.base={};self.current_file=None
        self.columnconfigure(0,weight=1);self.rowconfigure(1,weight=1)
        head=ttk.Frame(self,padding=12);head.grid(row=0,column=0,sticky='ew')
        ttk.Label(head,text='Вторичное облако · постоянная подача газа' if secondary else 'Первичное облако · ёмкость с газом',font=('Segoe UI',16,'bold')).pack(anchor='w')
        ttk.Label(head,text='HEGADAS · изотермический шлейф · экспериментальная реализация.' if secondary else 'Полное разрушение, без жидкой фазы и подпитки. Карты X–Y. Экспериментальная модель.').pack(anchor='w',pady=(4,0))
        ttk.Button(head,text='Открыть первичное облако' if secondary else 'Открыть вторичное облако',command=lambda: subprocess.Popen([sys.executable,str(ROOT/'gui.py')]+([] if secondary else ['--secondary']))).pack(anchor='e')
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
                self.confirmed=tk.BooleanVar()
                ttk.Checkbutton(page,text='Однофазность и применимость идеального газа для этих условий проверены',variable=self.confirmed).grid(row=row,column=0,columnspan=2,sticky='w',pady=15)
                ttk.Label(page,text='Расход задаётся после струевого участка, при атмосферном давлении.\nТемпература газа равна температуре воздуха и поверхности.' if secondary else 'Pабс = Pизб + Pатм. Свойства газа и пороги вводятся для выбранного вещества.\nСмена названия вещества сама по себе не меняет его свойства.',wraplength=850).grid(row=row+1,column=0,columnspan=2,sticky='w')
            if title=='Погода':
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
        foot=ttk.Frame(self,padding=12);foot.grid(row=2,column=0,sticky='ew');foot.columnconfigure(1,weight=1)
        ttk.Label(foot,text='Папка результатов').grid(row=0,column=0,sticky='w')
        self.output=tk.StringVar(value=str(ROOT/('results_secondary_user' if secondary else 'results_user')))
        ttk.Entry(foot,textvariable=self.output).grid(row=0,column=1,sticky='ew',padx=8)
        ttk.Button(foot,text='Выбрать…',command=self.choose_output).grid(row=0,column=2)
        buttons=ttk.Frame(foot);buttons.grid(row=1,column=0,columnspan=3,sticky='ew',pady=10)
        for title,command in [('Открыть JSON',self.open_input),('Сохранить JSON',self.save_input),('Пример',self.load_example),('Проверить',self.check_input)]:
            ttk.Button(buttons,text=title,command=command).pack(side='left',padx=(0,6))
        self.run_button=ttk.Button(buttons,text='Рассчитать',command=self.calculate);self.run_button.pack(side='left',padx=6)
        self.stop_button=ttk.Button(buttons,text='Остановить',command=self.stop,state='disabled');self.stop_button.pack(side='left')
        ttk.Button(buttons,text='Открыть результаты',command=self.open_results).pack(side='right')
        self.status=tk.StringVar(value='');ttk.Label(foot,textvariable=self.status,wraplength=1000).grid(row=2,column=0,columnspan=3,sticky='w')
        self.progress=ttk.Progressbar(foot,mode='indeterminate');self.progress.grid(row=3,column=0,columnspan=3,sticky='ew',pady=(8,0))
        self.protocol('WM_DELETE_WINDOW',self.close)
        self.load(self.example)

    def changed_physics(self,*_):
        if not self.loading and hasattr(self,'confirmed'):
            self.confirmed.set(False)

    def load(self,path):
        try:
            d=self.inputs.completed(json.loads(Path(path).read_text(encoding='utf-8-sig')))
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

    def open_input(self):
        p=filedialog.askopenfilename(filetypes=[('Исходные данные JSON','*.json')])
        if p:self.load(p)

    def load_example(self): self.load(self.example)

    def save_input(self):
        try:
            d=self.data()
            p=filedialog.asksaveasfilename(defaultextension='.json',initialfile='accident_input.json',filetypes=[('JSON','*.json')])
            if p:
                Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
                self.base=d;self.current_file=Path(p);self.status.set('Сохранено: '+p)
        except Exception as e:messagebox.showerror('Ошибка исходных данных',str(e))

    def check_input(self):
        try:
            d=self.data()
            if self.secondary:
                m=self.inputs.build_model(d)
                self.status.set(f'Проверка пройдена. Расход {m.q:g} кг/с; масса {m.q*m.feed.duration_s:g} кг; B = H = {m.initial_b:.3f} м.')
                return
            m=PrimaryCloud(Gas(**d['gas']),Vessel(**d['vessel']),Atmosphere(**d['atmosphere']),Options(**d['options']))
            a=m.initial
            self.status.set(f"Проверка пройдена. Масса {a['mass_kg']:.3f} кг; R = H = {a['radius_m']:.3f} м; начальная температура {a['temperature_k']-273.15:.2f} °C.")
        except Exception as e:messagebox.showerror('Ошибка исходных данных',str(e))

    def choose_output(self):
        p=filedialog.askdirectory()
        if p:self.output.set(p)

    def calculate(self):
        if self.process is not None:return
        try:
            d=self.data()
            if not self.output.get().strip():raise ValueError('Выберите папку результатов')
            parent=Path(self.output.get()).expanduser().resolve();parent.mkdir(parents=True,exist_ok=True)
            folder=Path(tempfile.mkdtemp(prefix=datetime.now().strftime('%Y%m%d_%H%M%S_'),dir=parent))
            inp=folder/'input.json';inp.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
            self.log_handle=(folder/'run.log').open('w',encoding='utf-8')
            env=dict(os.environ,PYTHONIOENCODING='utf-8')
            self.process=subprocess.Popen([sys.executable,str(self.runner),'--input',str(inp),'--output',str(folder)],stdout=self.log_handle,stderr=subprocess.STDOUT,env=env,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
            self.last_output=folder;self.stopped=False
            self.run_button.configure(state='disabled');self.stop_button.configure(state='normal');self.progress.start(12)
            self.status.set('Расчёт выполняется. Используется снимок введённых данных; изменения формы относятся к следующему запуску.')
            self.after(200,self.poll)
        except Exception as e:
            if self.log_handle:self.log_handle.close();self.log_handle=None
            messagebox.showerror('Не удалось запустить',str(e))

    def poll(self):
        code=self.process.poll()
        if code is None:self.after(200,self.poll);return
        self.log_handle.close();self.log_handle=None;self.process=None
        self.progress.stop();self.run_button.configure(state='normal');self.stop_button.configure(state='disabled')
        if self.stopped:
            self.status.set('Расчёт остановлен. В папке могут находиться неполные результаты: '+str(self.last_output))
        elif code:
            self.status.set('Ошибка расчёта. Подробности в run.log: '+str(self.last_output))
            messagebox.showerror('Ошибка расчёта',(self.last_output/'run.log').read_text(encoding='utf-8',errors='replace')[-5000:])
        else:
            r=json.loads((self.last_output/'result.json').read_text(encoding='utf-8'))
            self.status.set('Готово. '+('Досрочная остановка: '+r['stop_reason']+'. ' if r['stop_reason']!='requested_duration' else '')+str(self.last_output))

    def stop(self):
        if self.process is not None:self.stopped=True;self.process.terminate()

    def open_results(self):
        p=self.last_output or Path(self.output.get())
        if not p.exists():messagebox.showinfo('Результаты','Сначала выполните расчёт');return
        if os.name=='nt':os.startfile(str(p))
        else:subprocess.Popen(['open' if sys.platform=='darwin' else 'xdg-open',str(p)])

    def close(self):
        if self.process is not None:
            if not messagebox.askyesno('Закрыть','Остановить текущий расчёт и закрыть окно?'):return
            self.process.terminate();self.process.wait(timeout=5)
            if self.log_handle:self.log_handle.close()
        self.destroy()

if __name__=='__main__':App(secondary='--secondary' in sys.argv).mainloop()
