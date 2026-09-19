"""Lossless GUI drafts and validation of the selected calculation modes."""
import input_data
import secondary_input

MODES={'primary':('primary',),'secondary':('secondary',),'both':('primary','secondary')}
SCHEMA='cloud-project-v1'

def validate_project(project):
    if project.get('schema')!=SCHEMA or project.get('mode') not in MODES:
        raise ValueError('Неизвестный формат проекта или режим расчёта')
    outputs={}
    for name in MODES[project['mode']]:
        try:
            d=project['editors'][name]
            module=secondary_input if name=='secondary' else input_data
            values={(v['section'],v['key']):v['value'] for v in d['values']}
            outputs[name]=module.collect(d['base'],values,d['confirmed'],d['thresholds'],d['receptors'],d['snapshots'])
        except Exception as e:
            raise ValueError(('Первичное' if name=='primary' else 'Вторичное')+' облако: '+str(e)) from e
    return outputs

def validate_drafts(project):
    if project.get('schema')!=SCHEMA or project.get('mode') not in MODES:
        raise ValueError('Неизвестный формат проекта')
    # Structural validation only: incomplete input can be saved and reopened.
    for name,module in [('primary',input_data),('secondary',secondary_input)]:
        d=project['editors'][name]
        module.display_values(module.completed(d['base']))
        values={(v['section'],v['key']):v['value'] for v in d['values']}
        if set(values)!={(f[1],f[2]) for f in module.FIELDS} or not all(isinstance(v,str) for v in values.values()):
            raise ValueError('Неполный набор полей: '+name)
        if type(d['confirmed']) is not bool or not all(isinstance(d[k],str) for k in ['thresholds','receptors','snapshots']):
            raise ValueError('Некорректный формат редактора: '+name)
    return project
