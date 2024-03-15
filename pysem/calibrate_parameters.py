import numpy as np
import h5py
from parse_h5_traces import parse, components
import matplotlib.pyplot as plt

def calibrate_parameters():
    nt = 1000
    dt_u = np.zeros((nt, 3))
    dt_psi = np.zeros((nt, 3))
    dpos_u = np.zeros((nt, 3))
    dpos_psi = np.zeros((nt, 3))

if __name__ == "__main__":
    calibrate_parameters()
    opt = {'wkd':'./input_files/traces2/',
           'fmt':'h5',
           'nam':['Source'],
           'var':['Displ'],
           'rdr':['x','y','z'],
           'mon':[0, 1, 2]}
    stream = parse(**opt)
    print("Database parsed!")
    trace = list(stream.values())[0]
    plt.figsize=(10, 5)
    plt.title("Displacement in x direction at source after\nsoliving the backward problem.")
    plt.plot(trace.Time, trace.displ_values(0, 'x'))
    plt.savefig("backward_pb.png")
