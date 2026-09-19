import copy,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
import input_data,secondary_input
from project_data import validate_project,validate_drafts,SCHEMA
from gui import App
ROOT=Path(__file__).resolve().parents[1]

def example_project():
    p=dict(schema=SCHEMA,mode='both',editors={})
    for name,module,file in [('primary',input_data,'chloromethane.json'),('secondary',secondary_input,'secondary_gas.json')]:
        d=module.completed(json.loads((ROOT/'examples'/file).read_text()))
        p['editors'][name]=dict(base=d,values=[dict(section=s,key=k,value=v) for (s,k),v in module.display_values(d).items()],confirmed=True,
            thresholds='Порог; 0.1',receptors='',snapshots='30; 60')
    return p

class ProjectTests(unittest.TestCase):
    def test_all_modes_and_roundtrip(self):
        p=example_project()
        for mode,expected in [('primary',{'primary'}),('secondary',{'secondary'}),('both',{'primary','secondary'})]:
            p['mode']=mode;r=json.loads(json.dumps(p,ensure_ascii=False))
            self.assertEqual(validate_drafts(r),p);self.assertEqual(set(validate_project(r)),expected)

    def test_invalid_inactive_draft_preserved_but_not_calculated(self):
        p=example_project();p['mode']='primary';p['editors']['secondary']['values'][0]['value']=''
        self.assertEqual(validate_drafts(p),p);self.assertEqual(set(validate_project(p)),{'primary'})
        p['mode']='both'
        with self.assertRaisesRegex(ValueError,'Вторичное'):validate_project(p)

    def test_incomplete_project_rejected_before_restoring(self):
        p=example_project();p['editors']['secondary']['values'].pop()
        with self.assertRaises(ValueError):validate_drafts(p)

    def test_two_job_lifecycle_and_stop(self):
        # Exercise GUI orchestration without needing a window server.
        for stop in (False,True):
            with tempfile.TemporaryDirectory() as tmp:
                fake=SimpleNamespace(last_output=Path(tmp),queue=['primary','secondary'],stopped=False,
                    run_button=SimpleNamespace(configure=lambda **k:None),stop_button=SimpleNamespace(configure=lambda **k:None),
                    status=SimpleNamespace(set=lambda text:None),after=lambda *a:None,
                    manifest={'status':'running','calculations':{n:{'status':'pending'} for n in ['primary','secondary']}})
                fake.poll=lambda:App.poll(fake);fake.write=App.write;fake.next_job=lambda:App.next_job(fake);fake.finish=lambda state:App.finish(fake,state)
                for name in fake.queue:
                    folder=Path(tmp)/name;folder.mkdir();(folder/'result.json').write_text(json.dumps({'stop_reason':'requested_distance'}))
                with patch('gui.subprocess.Popen',return_value=SimpleNamespace(poll=lambda:0)) as popen:
                    fake.next_job();fake.stopped=stop;App.poll(fake)
                    if not stop:App.poll(fake)
                    self.assertEqual(popen.call_count,1 if stop else 2)
                    summary=json.loads((Path(tmp)/'summary.json').read_text())
                    self.assertEqual(summary['status'],'stopped' if stop else 'completed')
                    self.assertEqual(summary['calculations']['secondary']['status'],'not_run' if stop else 'completed')
