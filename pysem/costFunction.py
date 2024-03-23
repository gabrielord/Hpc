import numpy as np
import h5py
from parse_h5_traces import parse, components
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from gradient import gradient

def costFunction(mu,lamb,R_lamb,R_mu, g_reg_lamb,g_reg_mu, g_lamb_mis, g_mu_mis):
    # define options
    opt = {'syd':'./input_files/traces_forward/',
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
    totalMisfit = 0
    for m in range(len(opt['mon'])):
        for direction in ['x','y','z']:
            if dt_gt != dt_synth:
                interpolator = interp1d(trace_synthetic.Time,
                                        trace_synthetic.displ_values(m, direction),
                                        kind='linear', fill_value='extrapolate')
                synthetic_extrapolated = interpolator(trace_truth.Time)
                misfit = np.linalg.norm((synthetic_extrapolated[::-1],trace_truth.displ_values(m, direction)[::-1]))**2
                
            else:
                misfit = np.linalg.norm((trace_synthetic.displ_values(m, direction)[::-1] - trace_truth.displ_values(m, direction)[::-1]))**2
    
    totalMisfit = sum(misfit) 
    grad_lamb, grad_mu = gradient(mu,lamb,R_lamb,R_mu, g_reg_lamb,g_reg_mu, g_lamb_mis, g_mu_mis)
    costFunction = 1/2*(totalMisfit + R_lamb*np.dot(grad_lamb,grad_lamb) + R_mu*np.dot(grad_mu,grad_mu))



