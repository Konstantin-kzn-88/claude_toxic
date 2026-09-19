"""One window for independent primary, secondary, or both source calculations."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from datetime import datetime
import tkinter as tk
from tkinter import ttk,filedialog,messagebox
from input_panel import InputPanel,ROOT
from project_data import SCHEMA,MODES,validate_project,validate_drafts

LABELS={'primary':'Первичное','secondary':'Вторичное','both':'Оба облака'}

class App(tk.Tk):
    def __init__(self,secondary=False,both=False):
        super().__init__()
        self.title('Газовые облака — единый проект');self.geometry('1120x900');self.minsize(900,780)
        self.process=None;self.log=None;self.last_output=None;self.queue=[]
        self.columnconfigure(0,weight=1);self.rowconfigure(1,weight=1)
        head=ttk.Frame(self,padding=10);head.grid(row=0,column=0,sticky='ew')
        self.mode=tk.StringVar(value='both' if both else 'secondary' if secondary else 'primary')
        ttk.Label(head,text='Рассчитать:').pack(side='left')
        for key,label in LABELS.items():ttk.Radiobutton(head,text=label,variable=self.mode,value=key,command=self.switch).pack(side='left',padx=12)
        self.tabs=ttk.Notebook(self);self.tabs.grid(row=1,column=0,sticky='nsew')
        self.panels={}
        for name in ('primary','secondary'):
            p=InputPanel(self.tabs,secondary=name=='secondary');self.panels[name]=p;self.tabs.add(p,text=LABELS[name]+' облако')
        foot=ttk.Frame(self,padding=10);foot.grid(row=2,column=0,sticky='ew');foot.columnconfigure(1,weight=1)
        ttk.Label(foot,text='Два источника рассчитываются отдельно. Концентрации и зоны не объединяются.\nМассы источников задавайте без двойного учёта одного запаса вещества.',wraplength=1050).grid(row=0,column=0,columnspan=3,sticky='w')
        ttk.Label(foot,text='Папка результатов').grid(row=1,column=0)
        self.output=tk.StringVar(value=str(ROOT/'results_projects'))
        ttk.Entry(foot,textvariable=self.output).grid(row=1,column=1,sticky='ew',padx=8)
        ttk.Button(foot,text='Выбрать',command=self.choose_output).grid(row=1,column=2)
        buttons=ttk.Frame(foot);buttons.grid(row=2,column=0,columnspan=3,sticky='w',pady=6)
        for label,fn in [('Открыть JSON',self.open_input),('Сохранить проект',self.save),('Пример',self.example),('Проверить',self.check)]:
            ttk.Button(buttons,text=label,command=fn).pack(side='left',padx=3)
        self.run_button=ttk.Button(buttons,text='Рассчитать',command=self.calculate);self.run_button.pack(side='left',padx=3)
        self.stop_button=ttk.Button(buttons,text='Остановить',command=self.stop,state='disabled');self.stop_button.pack(side='left',padx=3)
        ttk.Button(buttons,text='Результаты',command=self.open_results).pack(side='left',padx=3)
        ttk.Button(foot,text='Скопировать вещество и погоду: первичное → вторичное',command=self.copy_common).grid(row=3,column=0,columnspan=3,sticky='w')
        self.status=tk.StringVar();ttk.Label(foot,textvariable=self.status,wraplength=1050).grid(row=4,column=0,columnspan=3,sticky='w',pady=6)
        self.protocol('WM_DELETE_WINDOW',self.close);self.switch()

    def switch(self):
        active=MODES[self.mode.get()]
        for name,p in self.panels.items():self.tabs.tab(p,state='normal' if name in active else 'disabled')
        self.tabs.select(self.panels[active[0]])

    def project(self):
        return dict(schema=SCHEMA,mode=self.mode.get(),editors={name:dict(base=p.base,
            values=[dict(section=k[0],key=k[1],value=v.get()) for k,v in p.vars.items()],
            confirmed=p.confirmed.get(),thresholds=p.thresholds.get('1.0','end-1c'),
            receptors=p.receptors.get('1.0','end-1c'),snapshots=p.snapshots.get()) for name,p in self.panels.items()})

    def restore(self,d):
        validate_drafts(d)
        for name,p in self.panels.items():
            draft=d['editors'][name];p.loading=True
            try:
                p.base=draft['base']
                for v in draft['values']:p.vars[(v['section'],v['key'])].set(v['value'])
                p.confirmed.set(draft['confirmed']);p.snapshots.set(draft['snapshots'])
                for key in ('thresholds','receptors'):
                    widget=getattr(p,key);widget.delete('1.0','end');widget.insert('1.0',draft[key])
            finally:p.loading=False
        self.mode.set(d['mode']);self.switch()

    def open_input(self):
        path=filedialog.askopenfilename(filetypes=[('JSON','*.json')])
        if not path:return
        try:
            d=json.loads(Path(path).read_text(encoding='utf-8-sig'))
            if d.get('schema')==SCHEMA:self.restore(d)
            else:
                name='secondary' if 'feed' in d else 'primary'
                p=self.panels[name];p.inputs.display_values(p.inputs.completed(d))
                p.load(path);self.mode.set(name);self.switch()
            self.status.set('Загружено: '+path)
        except Exception as e:messagebox.showerror('Не удалось открыть',str(e))

    def save(self):
        path=filedialog.asksaveasfilename(defaultextension='.json',initialfile='cloud_project.json',filetypes=[('JSON','*.json')])
        if path:
            try:
                self.write(Path(path),self.project());self.status.set('Проект сохранён, включая поля отключённого режима. Проверка выполняется перед расчётом.')
            except Exception as e:messagebox.showerror('Не удалось сохранить',str(e))

    @staticmethod
    def write(path,data):path.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')

    def example(self):
        for name in MODES[self.mode.get()]:self.panels[name].load(self.panels[name].example)
        self.status.set('Пример загружен для выбранных режимов.')

    def copy_common(self):
        a,b=self.panels['primary'],self.panels['secondary']
        for k,v in a.vars.items():
            if k[0] in ('gas','atmosphere') and k in b.vars:b.vars[k].set(v.get())
        b.confirmed.set(False)
        self.status.set('Вещество и погода скопированы. Подтвердите однофазность во вторичном режиме; его температура поверхности равна температуре воздуха.')

    def check(self):
        try:
            data=validate_project(self.project());self.status.set('Проверка пройдена: '+', '.join(LABELS[n] for n in data)+'.')
        except Exception as e:messagebox.showerror('Исходные данные',str(e))

    def choose_output(self):
        p=filedialog.askdirectory()
        if p:self.output.set(p)

    def calculate(self):
        if self.process is not None:return
        try:
            project=self.project();data=validate_project(project)
            if not self.output.get().strip():raise ValueError('Выберите папку результатов')
            root=Path(self.output.get()).expanduser().resolve();root.mkdir(parents=True,exist_ok=True)
            self.last_output=Path(tempfile.mkdtemp(prefix=datetime.now().strftime('%Y%m%d_%H%M%S_'),dir=root))
            self.write(self.last_output/'project.json',project)
            self.manifest=dict(mode=project['mode'],status='running',combined_concentrations=False,calculations={})
            self.queue=list(data);self.stopped=False
            for name,d in data.items():
                folder=self.last_output/name;folder.mkdir();self.write(folder/'input.json',d)
                self.manifest['calculations'][name]=dict(status='pending',directory=name)
            self.run_button.configure(state='disabled');self.stop_button.configure(state='normal')
            self.next_job()
        except Exception as e:
            self.run_button.configure(state='normal');self.stop_button.configure(state='disabled')
            messagebox.showerror('Не удалось запустить',str(e))

    def next_job(self):
        if not self.queue:
            self.finish('completed');return
        self.current=self.queue.pop(0);folder=self.last_output/self.current
        self.manifest['calculations'][self.current]['status']='running';self.write(self.last_output/'summary.json',self.manifest)
        self.log=(folder/'run.log').open('w',encoding='utf-8')
        runner='run_secondary.py' if self.current=='secondary' else 'run.py'
        try:
            self.process=subprocess.Popen([sys.executable,str(ROOT/runner),'--input',str(folder/'input.json'),'--output',str(folder)],stdout=self.log,stderr=subprocess.STDOUT,env=dict(os.environ,PYTHONIOENCODING='utf-8'),creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        except Exception as e:
            self.log.close();self.log=None;self.manifest['calculations'][self.current]['status']='failed';self.finish('failed');messagebox.showerror('Запуск',str(e));return
        self.status.set('Рассчитывается '+LABELS[self.current].lower()+' облако. Изменения полей относятся к следующему запуску.')
        self.after(200,self.poll)

    def poll(self):
        code=self.process.poll()
        if code is None:self.after(200,self.poll);return
        self.log.close();self.log=None;self.process=None
        state='stopped' if self.stopped else 'failed' if code else 'completed'
        row=self.manifest['calculations'][self.current];row['status']=state;row['exit_code']=code
        if state=='completed':
            try:
                result=json.loads((self.last_output/self.current/'result.json').read_text(encoding='utf-8'))
                row['stop_reason']=result['stop_reason'];row['warnings']=result.get('warnings',[])
            except Exception as e:
                row['status']='failed';row['error']=str(e);self.finish('failed');return
            self.next_job()
        else:self.finish(state)

    def finish(self,state):
        for name in self.queue:self.manifest['calculations'][name]['status']='not_run'
        self.queue=[];self.manifest['status']=state;self.write(self.last_output/'summary.json',self.manifest)
        self.run_button.configure(state='normal');self.stop_button.configure(state='disabled')
        self.status.set({'completed':'Расчёты завершены. Результаты раздельные: ','failed':'Ошибка: смотрите summary.json и run.log. ','stopped':'Остановлено; результаты могут быть неполными. '}[state]+str(self.last_output))

    def stop(self):
        if self.process is not None:self.stopped=True;self.process.terminate()

    def open_results(self):
        p=self.last_output or Path(self.output.get())
        if not p.exists():messagebox.showinfo('Результаты','Сначала выполните расчёт');return
        if os.name=='nt':os.startfile(str(p))
        else:subprocess.Popen(['open' if sys.platform=='darwin' else 'xdg-open',str(p)])

    def close(self):
        if self.process is not None:
            if not messagebox.askyesno('Закрыть','Остановить текущий расчёт?'):return
            self.stop()
            try:self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:self.process.kill();self.process.wait()
            self.poll()
        self.destroy()

if __name__=='__main__':App(secondary='--secondary' in sys.argv,both='--both' in sys.argv).mainloop()
