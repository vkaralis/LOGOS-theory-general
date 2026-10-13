"""Generate simulated subjects and noisy concentrations; no clinical records."""
import csv
from pathlib import Path
import numpy as np

SEED=20261003

def generate(output, subjects=24, seed=SEED):
    if subjects<1:
        raise ValueError('subjects must be positive')
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(seed)
    with out.open('w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['synthetic_id','time','dose','k','volume','true_concentration','simulated_measurement','is_synthetic'])
        for i in range(subjects):
            k=float(rng.lognormal(0,.15));v=float(rng.lognormal(0,.1));dose=1.
            for t in np.linspace(0,6,25):
                c=dose/v*np.exp(-k*t);obs=max(0.,c+float(rng.normal(0,.03)))
                writer.writerow([f'SYN-{i+1:03d}',f'{t:.8f}',dose,f'{k:.12f}',f'{v:.12f}',f'{c:.12f}',f'{obs:.12f}','true'])
