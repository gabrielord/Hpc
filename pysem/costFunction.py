import numpy as np
import h5py
from parse_h5_traces import parse, components
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import json

def costFunction(mu,lamb,R_lamb,R_mu):
    # define options
    opt = {'syd':'traces_forward/',
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
    grad_lamb = np.gradient(lamb)
    rad_lamb = np.array(grad_lamb)
    dot_p_lamb = grad_lamb[0] * grad_lamb[1]
    res_grad_lamb  = sum(dot_p_lamb.flatten())
    grad_mu = np.gradient(mu)
    grad_mu = np.array(grad_mu)
    dot_p_mu = grad_mu[0] * grad_mu[1]
    res_grad_mu  = sum(dot_p_mu.flatten())
    cost = 1/2*(totalMisfit + R_lamb*res_grad_lamb + R_mu*res_grad_mu)
    with open('cout.json') as json_file:
        json.dump(cost, json_file)
    return cost

