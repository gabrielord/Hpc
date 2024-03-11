import numpy as np
import h5py
from parse_h5_traces import parse, components
import matplotlib.pyplot as plt

def compute_misfit(trace_truth, trace_synthetic, m):
    """Compute the misfit between the synthetic and the grand truth traces for a given monitor."""
    misfit = 0
    for direction in ['x','y','z']:
        misfit += np.linalg.norm(trace_truth.displ_values(m, direction) - trace_synthetic.displ_values(m, direction))
    return misfit

def generate_misfit_files():
    # define options
    opt = {'syd':'./input_files/traces/',
           'trd':'./Uobs/',
           'fmt':'h5',
           'nam':['all'],
           'var':['Displ'],
           'rdr':['x','y','z'],
           'mon':[0, 1, 2]}
    
    # parse database
    opt_truth = opt.copy()
    opt_truth['wkd'] = opt['trd']
    opt_synthetic = opt.copy()
    opt_synthetic['wkd'] = opt['syd']

    stream_truth = parse(**opt_truth)
    print("Grand truth database parsed!")
    stream_synthetic = parse(**opt_synthetic)
    print("Synthetic database parsed!\n")

    # compute misift
    trace_truth = list(stream_truth.values())[0]
    trace_synthetic = list(stream_synthetic.values())[0]

    dt_gt = trace_truth.Time[1] - trace_truth.Time[0]
    dt_synth = trace_synthetic.Time[1] - trace_synthetic.Time[0]

    if dt_gt != dt_synth:
        raise ValueError("The time delay between two successive dates for grand truth and synthetic traces are different!")

    for m in range(3):
        for direction in ['x','y','z']:
            print(f"Writing misfit file for monitor {m} and direction {direction}", end='\t')
            misfit = trace_synthetic.displ_values(m, direction)[::-1] - trace_truth.displ_values(m, direction)[::-1]
            with (open (f"input_files/misfit_{m}_{direction}.txt", "w")) as f:
                for i in range(len(misfit)):
                    f.write(f"{trace_truth.Time[i]}, {misfit[i]}\n")
            print("Done!")


if __name__ == "__main__":
    generate_misfit_files()