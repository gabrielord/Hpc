import numpy as np
import h5py
from parse_h5_traces import parse, components
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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

    for m in range(3):
        for direction in ['x','y','z']:
            print(f"Writing misfit file for monitor {m} and direction {direction}", end='\t')
            if dt_gt != dt_synth:
                interpolator = interp1d(trace_synthetic.Time,
                                        trace_synthetic.displ_values(m, direction),
                                        kind='linear', fill_value='extrapolate')
                synthetic_extrapolated = interpolator(trace_truth.Time)
                misfit = synthetic_extrapolated[::-1] - trace_truth.displ_values(m, direction)[::-1]
                print("The time delay between two successive dates for grand truth and synthetic traces are different.")
            else:
                misfit = trace_synthetic.displ_values(m, direction)[::-1] - trace_truth.displ_values(m, direction)[::-1]
            with (open (f"input_files/monitors_misfit/misfit_{m}_{direction}.txt", "w")) as f:
                for i in range(len(misfit)):
                    f.write(f"{trace_truth.Time[i]}, {misfit[i]}\n")
            print("Done!")

            plt.clf()
            plt.figure(figsize=(10, 10))
            plt.subplot(2, 1, 1)
            plt.title(f"Grand truth and synthetic traces for monitor {m} and direction {direction}")
            plt.plot(trace_truth.Time, trace_truth.displ_values(m, direction), "-g", label='Grand truth', alpha=0.6)
            if synthetic_extrapolated is not None:
                plt.plot(trace_truth.Time, synthetic_extrapolated, "-b", label='Synthetic (extrapolated)', alpha=0.6)
            else:
                plt.plot(trace_truth.Time, trace_synthetic.displ_values(m, direction), "-r", label='Synthetic', alpha=1)
            plt.legend()
            plt.margins(0)

            plt.subplot(2, 1, 2)
            plt.title(f"Misfit between grand truth and synthetic traces for monitor {m} and direction {direction}")
            plt.plot(trace_truth.Time, misfit[::-1], "-r", label='Misfit', alpha=0.8)
            plt.legend()

            plt.savefig(f"input_files/monitors_misfit/plot_{m}_{direction}.png")


if __name__ == "__main__":
    generate_misfit_files()