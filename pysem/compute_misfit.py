import numpy as np
import h5py
from parse_h5_traces import parse, components
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import time
import os

def plot_misfit(m, direction, trace_truth, trace_synthetic, misfit):
    """
    Function to display the misfit for monitor m according to a given direction
    """
    plt.clf()
    fig, ax = plt.subplots(2,1)
    ax[0].set_title(f"Grand truth and synthetic traces for monitor {m} and direction {direction}")
    ax[0].plot([0,2], [0,0], "k--")
    ax[0].plot(trace_truth.Time, trace_truth.displ_values(m, direction), "-g", label='Grand truth', alpha=0.8)
    ax[0].plot(trace_truth.Time, trace_synthetic, "-r", label='Synthetic', alpha=0.8)

    ax[0].legend()

    ax[1].set_title(f"Misfit for monitor {m} and direction {direction}")
    ax[1].plot([0,2], [0,0], "k--")
    ax[1].plot(trace_truth.Time, misfit[::-1], "-r", label='Misfit', alpha=0.8)
    ax[1].legend()

    fig.tight_layout()

    os.makedirs("monitors_misfit", exist_ok=True)
    print(f"monitors_misfit/plot_{m}_{direction}.png")
    plt.savefig(f"monitors_misfit/plot_{m}_{direction}.png")




def generate_misfit_files():
    # define options
    opt = {'syd':'forward/traces/',
           'trd':'../Uobs/',
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
    stream_synthetic = parse(**opt_synthetic)

    # compute misift
    trace_truth = list(stream_truth.values())[0]
    trace_synthetic = list(stream_synthetic.values())[0]

    dt_gt = trace_truth.Time[1] - trace_truth.Time[0]
    dt_synth = trace_synthetic.Time[1] - trace_synthetic.Time[0]

    for m in range(3):
        for direction in ['x','y','z']:
            print(f"Writing misfit file for monitor {m} and direction {direction}", end='\t')
            if dt_gt != dt_synth:
                interpolator = interp1d(trace_synthetic.Time,
                                        trace_synthetic.displ_values(m, direction),
                                        kind='linear', fill_value='extrapolate')
                synthetic_extrapolated = interpolator(trace_truth.Time)
                misfit = synthetic_extrapolated[::-1] - trace_truth.displ_values(m, direction)[::-1]
            else:
                misfit = trace_synthetic.displ_values(m, direction)[::-1] - trace_truth.displ_values(m, direction)[::-1]
            if not os.path.exists("monitors_misfit"):
            # Si le dossier n'existe pas, créez-le
                os.makedirs("monitors_misfit")
            with (open (f"monitors_misfit/misfit_{m}_{direction}.txt", "w")) as f:
                for i in range(len(misfit)):
                    f.write(f"{trace_truth.Time[i]}, {misfit[i]}\n")
                    
            
            if dt_gt != dt_synth:
                plot_misfit(m, direction, trace_truth,
                synthetic_extrapolated, misfit)
            else:
                plot_misfit(m, direction, trace_truth,
                trace_synthetic.displ_values(m, direction), misfit)
            
            print('Done !')


if __name__ == "__main__":
    generate_misfit_files()