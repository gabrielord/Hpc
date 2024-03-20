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