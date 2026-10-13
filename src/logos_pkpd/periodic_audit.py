"""Full factorial and local/stress-grid regression audit."""
import json
from pathlib import Path
from .periodic import cycle, constant, search

def run(output):
    OUT=Path(output); OUT.mkdir(parents=True, exist_ok=True)
    checks=[]
    def check(name,actual,expected,tol):
        error=abs(float(actual)-float(expected))
        checks.append(dict(check=name,actual=float(actual),expected=float(expected),error=error,tolerance=tol,passed=error<=tol))
    reported=[(.08,.170721,.944709,.944598,.520267),(3.159490,.111507,.766421,.490731,1.649774),(4.442336,.085510,1.789967,.862802,2.247619),(2.870281,.091278,1.738049,1.124215,1.521377)]
    factorial=[]
    for case,(sig,dep),target in zip('ABCD',[(False,False),(True,False),(False,True),(True,True)],reported):
     v=cycle(target[0],sig,dep,N=40001)
     for name,expected in zip(['J','peak_Z','mean_Z','peak_C'],target[1:]):check('Periodic reference '+case+' '+name,v[name],expected,8e-7)
     check('Periodic '+case+' periodic closure',v['closure'],0,1e-10)
     choose,rows=search(sig,dep);best=choose(3.5)
     check('Periodic '+case+' search J',best['J'],target[1],8e-7)
     factorial.append(dict(case=case,constant=constant(sig,dep),at_reported_interval=v,new_search=best,scan=rows))
     print('Independent factorial',case,'complete',flush=True)
    for v, expected_J, expected_Z in zip(factorial, [.170732,.084337,.057260,.073656], [.945494,.945494,17.346653,17.346653]):
        check('Constant '+v['case']+' J',v['constant']['J'],expected_J,8e-7)
        check('Constant '+v['case']+' injury',v['constant']['Z'],expected_Z,8e-7)
    local=[]
    for r0 in [.56,.7,.84]:
     for rho0 in [.8,1.,1.2]:
      choose,rows=search(True,True,r0,rho0);ref=constant(True,True,r0,rho0)
      for cap in [1.,2.8,3.5,4.2,20.]:
       v=choose(cap);local.append(dict(r0=r0,rho0=rho0,cap=cap,constant=ref,constant_safe=ref['Z']<=cap,cycle=v))
      print('Independent local search',r0,rho0,'complete',flush=True)
    core=[v for v in local if v['cap'] in [2.8,3.5,4.2]]
    check('Periodic 27 local constant infeasible',sum(not v['constant_safe'] for v in core),27,0)
    check('Periodic 27 local selected feasible',sum(v['cycle'] is not None for v in core),27,0)
    for key,minimum,maximum,tol in [('J',.083470,.097223,8e-7),('peak_Z',1.583994,1.998673,2e-5),('mean_Z',.933695,1.413159,2e-5),('T',2.764844,3.012607,2e-5)]:
     vals=[v['cycle'][key] for v in core];check('Periodic local minimum '+key,min(vals),minimum,tol);check('Periodic local maximum '+key,max(vals),maximum,tol)
    check('Periodic cap one no cycle detected',sum(v['cycle'] is None for v in local if v['cap']==1),9,0)
    check('Periodic cap twenty constant feasible',sum(v['constant_safe'] for v in local if v['cap']==20),6,0)
    check('Periodic cap twenty cycle feasible',sum(v['cycle'] is not None for v in local if v['cap']==20),9,0)
    record=dict(all_passed=all(v['passed'] for v in checks),scope='Independent reconstruction from declared model equations, 121-point searches and high-resolution integrating-factor quadrature; not original-source execution or interval arithmetic.',checks=checks,factorial=factorial,local=local)
    (OUT/'periodic_reference.json').write_text(json.dumps(record,indent=2))
    print('Periodic audit',len(checks),'checks;',sum(not c['passed'] for c in checks),'failures',flush=True)
    if not record['all_passed']:
        raise RuntimeError('Independent periodic audit failed; see the saved checks.')
    return record
