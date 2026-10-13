"""Compare newly recomputed results with the release reference at tolerances."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('output',type=Path)
args=parser.parse_args()
fresh=args.output;ref=ROOT/'reference'
checks=0

def close(actual,expected,tolerance,label):
    global checks
    checks+=1
    if not np.isfinite(actual) or abs(actual-expected)>tolerance:
        raise RuntimeError(f'{label}: {actual} != {expected} at {tolerance}')

for folder in [fresh,ref]:
    summary=json.loads((folder/'verification.json').read_text())
    if summary.get('all_passed') is not True:raise RuntimeError(f'Failed verification in {folder}')
    for name in ['periodic_reference.json','independent_benchmarks.json']:
        result=json.loads((folder/name).read_text())
        if result.get('all_passed') is not True:raise RuntimeError(f'Failed audit in {folder/name}')
        if not all(x['passed'] for x in result['checks']):raise RuntimeError('Failed individual audit check')
a=json.loads((fresh/'finite.json').read_text());b=json.loads((ref/'finite.json').read_text())
if a['seed']!=b['seed'] or a['timing_starts']!=24:raise RuntimeError('Finite seed/search mismatch')
for x,y in zip(a['rows'],b['rows'],strict=True):
    if x['regimen']!=y['regimen']:raise RuntimeError('Regimen mismatch')
    close(x['J'],y['J'],5e-9,x['regimen'])
a=json.loads((fresh/'periodic_reference.json').read_text());b=json.loads((ref/'periodic_reference.json').read_text())
for x,y in zip(a['factorial'],b['factorial'],strict=True):
    for key in ['J','peak_Z','mean_Z','peak_C']:
        close(x['new_search'][key],y['new_search'][key],2e-5,f"{x['case']} {key}")
    for xs,ys in zip(x['scan'],y['scan'],strict=True):
        for key in ['J','peak_Z','mean_Z','peak_C']:close(xs[key],ys[key],1e-7,f'Scan {key}')
for x,y in zip(a['local'],b['local'],strict=True):
    if (x['r0'],x['rho0'],x['cap'],x['constant_safe'],x['cycle'] is None)!=(y['r0'],y['rho0'],y['cap'],y['constant_safe'],y['cycle'] is None):
        raise RuntimeError('Sensitivity classification mismatch')
    if x['cycle']:
        for key in ['J','peak_Z','mean_Z','peak_C']:close(x['cycle'][key],y['cycle'][key],2e-5,f'Local {key}')
def recursive(actual,expected,label):
    if isinstance(expected,dict):
        if actual.keys()!=expected.keys():raise RuntimeError(f'{label}: key mismatch')
        for key in expected:recursive(actual[key],expected[key],label+'.'+key)
    elif isinstance(expected,list):
        for i,(x,y) in enumerate(zip(actual,expected,strict=True)):recursive(x,y,f'{label}[{i}]')
    elif isinstance(expected,bool) or expected is None or isinstance(expected,str):
        if actual!=expected:raise RuntimeError(f'{label}: metadata/classification mismatch')
    elif isinstance(expected,(int,float)):
        close(actual,expected,2e-5,label)

for name in ['A_constant_family.json','A_robustness.json','A_imperfect_information.json','A_finite_constant_family.json','B_dose_frontier.json','B_frontier_final_verification.json','A_mathematical_verification.json']:
    recursive(json.loads((fresh/name).read_text()),json.loads((ref/name).read_text()),name)

if (fresh/'synthetic_pk.csv').read_bytes()!=(ROOT/'data/synthetic_pk.csv').read_bytes():raise RuntimeError('Synthetic data mismatch')
for name in ['finite_regimens','regular_bolus_scans','factorial_cycles','sensitivity_map','cessation_tail','dose_frontier']:
    for ext in ['png','svg']:
        if (fresh/f'{name}.{ext}').stat().st_size==0:raise RuntimeError('Missing figure')
print(f'PASS: {checks} numerical comparisons, all audit checks and synthetic data identity')
