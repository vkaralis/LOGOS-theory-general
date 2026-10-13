"""Run all public numerical examples and verification from a fresh directory."""
import argparse
import csv
import json
from pathlib import Path
import platform
import re
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from scipy.integrate import solve_ivp
from . import __version__, finite, tail, synthetic, periodic_audit, benchmarks, frontier, recovery, recovery_finite
from .periodic import rates

def save(fig, out, name):
    fig.tight_layout()
    for extension in ('png','svg'):
        path = out/f'{name}.{extension}'
        metadata = {'Software': None} if extension == 'png' else {'Creator': None, 'Date': None}
        fig.savefig(path, dpi=160, metadata=metadata)
        if extension == 'png':
            with Image.open(path) as image:
                image.save(path, format='PNG')
        if extension == 'svg':
            path.write_text(re.sub(r'<metadata>.*?</metadata>\s*', '', path.read_text(), flags=re.DOTALL))
    plt.close(fig)

def table(path, rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def reproduce(output):
    out=Path(output)
    if out.exists():
        raise FileExistsError('Choose a new output directory')
    out.mkdir(parents=True)
    result=finite.run()
    (out/'finite.json').write_text(json.dumps(result,indent=2)+'\n')
    table(out/'finite.csv',[{k:r[k] for k in ['regimen','J','adjacent_gain_percent','C_upper_bound','injury_upper_bound_full_tail']} for r in result['rows']])
    reference=np.array([.629606362035524,.672945396127757,.677043273327298,.679552690590195])
    if not np.allclose([r['J'] for r in result['rows']],reference,rtol=0,atol=5e-9):
        raise RuntimeError('Finite reference mismatch')
    # Independent integration split at all dose jumps and infusion transitions.
    errors=[]
    for row in result['rows']:
        if row['regimen']=='expanded':
            transitions=[0,row['maintenance'],3];state=np.array([row['concentration'],0.])
        else:
            transitions=sorted(set([0,*row['times'],3]));state=np.zeros(2)
        for a,b in zip(transitions[:-1],transitions[1:]):
            if row['regimen']!='expanded':
                state[0]+=sum(d for d,t in zip(row['doses'],row['times']) if t==a)
            u=row['concentration'] if row['regimen']=='expanded' and a<row['maintenance'] else 0.
            sol=solve_ivp(lambda t,y:[u-y[0],y[0]/(1+y[0])],(a,b),state,rtol=1e-11,atol=1e-13)
            if not sol.success:raise RuntimeError(sol.message)
            state=sol.y[:,-1]
        errors.append(abs(state[1]-row['J']))
    if max(errors)>2e-10:raise RuntimeError('Independent finite integration mismatch')
    t=np.linspace(0,3,3001);fig,axes=plt.subplots(1,2,figsize=(10,4))
    trajectory_rows=[]
    for row in result['rows']:
        c=finite.trajectory(row,t)
        # Exact accumulated benefit avoids trapezoids across concentration jumps.
        if row['regimen']=='expanded':
            hold=row['maintenance'];p=row['concentration']
            j=np.minimum(t,hold)*p/(1+p)+np.where(t>hold,np.log1p(p)-np.log1p(p*np.exp(-np.maximum(t-hold,0))),0)
        else:
            j=np.array([finite.objective([d for d,s in zip(row['doses'],row['times']) if s<=tt],[s for s in row['times'] if s<=tt],tt) for tt in t])
        axes[0].plot(t,c,label=row['regimen']);axes[1].plot(t,j,label=row['regimen'])
        trajectory_rows += [dict(regimen=row['regimen'],time=float(tt),C=float(cc),J=float(jj)) for tt,cc,jj in zip(t,c,j)]
    table(out/'finite_trajectories.csv',trajectory_rows)
    for ax,y in zip(axes,['Concentration','Cumulative Emax benefit']):ax.set(xlabel='Time',ylabel=y);ax.legend()
    save(fig,out,'finite_regimens')
    periodic=periodic_audit.run(out)
    table(out/'factorial.csv',[dict(case=v['case'],constant_J=v['constant']['J'],constant_Z=v['constant']['Z'],**v['new_search']) for v in periodic['factorial']])
    table(out/'sensitivity.csv',[dict(r0=v['r0'],rho0=v['rho0'],cap=v['cap'],constant_safe=v['constant_safe'],cycle_detected=v['cycle'] is not None,J=v['cycle']['J'] if v['cycle'] else '',interval=v['cycle']['T'] if v['cycle'] else '') for v in periodic['local']])
    comparison_rows=[]
    for v in periodic['factorial']:
        for strategy in ['constant','bolus']:
            r=v['constant'] if strategy=='constant' else v['new_search']
            comparison_rows.append(dict(case=v['case'],strategy=strategy,mean_benefit=r['J'],peak_injury=r['Z'] if strategy=='constant' else r['peak_Z'],mean_injury=r['Z'] if strategy=='constant' else r['mean_Z'],peak_C=.5 if strategy=='constant' else r['peak_C'],interval='' if strategy=='constant' else r['T'],feasible=r['Z']<=3.5 if strategy=='constant' else r['peak_Z']<=3.5 and r['peak_C']<=4,resource_rule='exact amount per cycle; mean input 0.5'))
    table(out/'periodic_comparison.csv',comparison_rows)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,v in zip(axes,[periodic['factorial'][0],periodic['factorial'][3]]):
        scan=v['scan'];tt=np.array([x['T'] for x in scan]);jj=np.array([x['J'] for x in scan]);safe=np.array([x['peak_Z']<=3.5 and x['peak_C']<=4 for x in scan])
        ax.plot(tt,jj,color='gray',label='Sampled bolus family');ax.plot(tt,np.where(safe,jj,np.nan),color='tab:blue',label='Periodic caps satisfied')
        ax.axhline(v['constant']['J'],ls='--',color='tab:orange',label='Constant (safe)' if v['constant']['Z']<=3.5 else 'Constant (unsafe)')
        ax.scatter([v['new_search']['T']],[v['new_search']['J']],color='tab:blue')
        ax.set(title='Concave, independent rates' if v['case']=='A' else 'Sigmoidal, dependent rates',xlabel='Cycle interval',ylabel='Mean benefit');ax.legend(fontsize=8)
    save(fig,out,'regular_bolus_scans')
    fig,axes=plt.subplots(2,2,figsize=(10,7))
    independent_periodic_errors=[]
    for ax,v,(sig,dep) in zip(axes.flat,periodic['factorial'],[(False,False),(True,False),(False,True),(True,True)]):
        scan=v['scan'];tt=np.array([x['T'] for x in scan]);jj=np.array([x['J'] for x in scan]);safe=np.array([x['peak_Z']<=3.5 and x['peak_C']<=4 for x in scan])
        ax.plot(tt,jj,color='gray');ax.plot(tt,np.where(safe,jj,np.nan),color='tab:blue',label='Feasible sampled cycle')
        ax.axhline(v['constant']['J'],ls='--',color='tab:orange',label='Constant (safe)' if v['constant']['Z']<=3.5 else 'Constant (unsafe)')
        best=v['new_search'];ax.scatter([best['T']],[best['J']],color='tab:blue');ax.set(title=f"Case {v['case']}",xlabel='Cycle interval',ylabel='Mean benefit');ax.legend(fontsize=8)
        T=best['T'];c0=.5*T/(-np.expm1(-T))
        # Independently solve affine maps as augmented ODEs, then integrate benefit/injury.
        def maps(t,y):
            q,r,d,drive,repair=[float(x[0]) for x in rates(np.array([c0*np.exp(-t)]),sig,dep,.7,1.)]
            return [-(r+d)*y[0],r-(r+d)*y[1],-repair*y[2],drive-repair*y[3]]
        m=solve_ivp(maps,(0,T),[1,0,1,0],rtol=2e-12,atol=2e-14)
        if not m.success:raise RuntimeError(m.message)
        A,B,Z,W=m.y[:,-1];initial=[B/(1-A),W/(1-Z),0,0]
        def dynamics(t,y):
            q,r,d,drive,repair=[float(x[0]) for x in rates(np.array([c0*np.exp(-t)]),sig,dep,.7,1.)]
            return [r-(r+d)*y[0],drive-repair*y[1],y[0]*q,y[1]]
        sol=solve_ivp(dynamics,(0,T),initial,rtol=2e-12,atol=2e-14,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        err=max(abs(sol.y[2,-1]/T-best['J']),abs(sol.y[3,-1]/T-best['mean_Z']))
        independent_periodic_errors.append(err)
    if max(independent_periodic_errors)>1e-7:raise RuntimeError('Periodic independent ODE mismatch')
    save(fig,out,'factorial_cycles')
    fig,ax=plt.subplots(figsize=(9,5));cells=periodic['local']; labels=[];values=[]
    for v in cells:
        labels.append(f"{v['r0']}/{v['rho0']}/{v['cap']}")
        values.append(int(v['constant_safe'])+2*int(v['cycle'] is not None))
    for row in range(9):
        for col in range(5):
            v=cells[row*5+col];text=f"{v['cycle']['J']:.4f}\nT={v['cycle']['T']:.3f}" if v['cycle'] else 'No cycle'
            ax.text(col,row,text,ha='center',va='center',fontsize=7,color='white' if values[row*5+col]<2 else 'black')
    image=ax.imshow(np.array(values).reshape(9,5),vmin=0,vmax=3,cmap=plt.get_cmap('viridis',4),aspect='auto')
    ax.set_xticks(range(5),['1 (stress)','2.8','3.5','4.2','20 (stress)']);ax.set_yticks(range(9),[f"r0={v['r0']}, rho0={v['rho0']}" for v in cells[::5]]);ax.set(xlabel='Injury cap',title='Periodic feasibility: sampled family')
    cb=fig.colorbar(image,ax=ax,ticks=[0,1,2,3]);cb.ax.set_yticklabels(['Neither detected','Constant only','Cycle only','Both'])
    save(fig,out,'sensitivity_map')
    t=np.linspace(0,10,2001);z=tail.injury(t);p=tail.peak()
    table(out/'cessation_tail.csv',[dict(time=float(a),injury=float(b)) for a,b in zip(t,z)])
    fig,ax=plt.subplots(figsize=(7,4));ax.plot(t,z,label='Injury');ax.axhline(.1,color='tab:red',ls='--',label='Cap 0.1');ax.scatter([p['time']],[p['value']],label='Exact future peak');ax.set(xlabel='Time since cessation',ylabel='Injury');ax.legend();save(fig,out,'cessation_tail')
    synthetic.generate(out/'synthetic_pk.csv')
    recovery_record=recovery.run(out)
    recovery_finite_record=recovery_finite.run(out)
    frontier_rows=frontier.run(out)
    fig,ax=plt.subplots(figsize=(7,4))
    ax.plot([v['target'] for v in frontier_rows],[v['minimum_dose_candidate'] for v in frontier_rows],'o-',label='Flexible numerical candidate')
    feasible=[v for v in frontier_rows if v['constant_feasible']]
    ax.plot([v['target'] for v in feasible],[v['constant_minimum_dose'] for v in feasible],'s--',label='Constant family (feasible targets)')
    ax.set(xlabel='Required quadratic-model benefit',ylabel='Administered amount');ax.legend();save(fig,out,'dose_frontier')
    bench=benchmarks.run(out)
    summary=dict(version=__version__,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,all_passed=True,finite_ode_max_error=max(errors),periodic_ode_max_error=max(independent_periodic_errors),periodic_checks=len(periodic['checks']),analytical_checks=len(bench['checks']),recovery_checks=len(recovery_record['checks'])+len(recovery_finite_record['checks']),frontier_checks=20,tail_peak=p,data_kind='synthetic; dimensionless; no clinical observations')
    (out/'verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    run=sub.add_parser('reproduce');run.add_argument('--output',type=Path,default=Path('results'))
    gen=sub.add_parser('generate-data');gen.add_argument('--output',type=Path,required=True);gen.add_argument('--subjects',type=int,default=24)
    args=parser.parse_args()
    if args.command=='reproduce':reproduce(args.output)
    else:synthetic.generate(args.output,args.subjects)

if __name__=='__main__':main()
