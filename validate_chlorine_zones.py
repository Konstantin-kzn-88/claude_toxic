"""Compare dose-zone geometry at 1800 s; independent of fitting targets."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from scipy.integrate import trapezoid
ROOT = Path(__file__).resolve().parent
from cloud import *
def main():
    d=json.loads((ROOT/'validation/primary_chlorine/input.json').read_text(encoding='utf-8'))
    m=PrimaryCloud(Gas(**d['gas']),Vessel(**d['vessel']),Atmosphere(**d['atmosphere']),Options(**d['options']));m.run()
    results=[]
    for dt in [.25,.125]:
     ts=np.arange(0,1800+dt/2,dt);s=[m.state(y) for y in m.solution.sol(ts).T]
     cc=np.array([a['centre_concentration_kg_m3'] for a in s]);xc=np.array([a['centre_x_m'] for a in s]);rc2=np.array([a['core_radius_m']**2 for a in s]);sy2=np.array([a['sy2_m2'] for a in s])
     def dose(x,y=0):return trapezoid(cc*np.exp(-np.maximum((x-xc)**2+y*y-rc2,0)/sy2),ts)
     for label,threshold,leftref,rightref,widthref in [('PCt',.036,-94.94,1396.48,287.65),('LCt',.36,-92.89,356.69,189.82)]:
      left=brentq(lambda x:dose(x)-threshold,-1000,0);right=brentq(lambda x:dose(x)-threshold,0,3000)
      def width(x):return brentq(lambda y:dose(x,y)-threshold,0,1000) if dose(x)>threshold else 0.
      # Scan to locate the widest connected portion, then refine that interval.
      xx=np.linspace(left,right,101);ww=np.array([width(x) for x in xx]);j=int(ww.argmax())
      opt=minimize_scalar(lambda x:-width(x),bounds=(xx[max(j-1,0)],xx[min(j+1,100)]),method='bounded')
      result=dict(dt=dt,label=label,left_m=left,right_m=right,half_width_m=-opt.fun,reference_left_m=leftref,reference_right_m=rightref,reference_half_width_m=widthref,right_error_percent=100*(right/rightref-1),width_error_percent=100*(-opt.fun/widthref-1));results.append(result);print(result)
    out=ROOT/'results_chlorine_validation';out.mkdir(exist_ok=True)
    (out/'zones.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8')

if __name__ == '__main__':
    main()
