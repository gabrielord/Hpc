import numpy as np
import h5py
import argparse
from parse_h5_traces import parse, components
import matplotlib.pyplot as plt

if __name__=='__main__':
    parser = argparse.ArgumentParser(prefix_chars='@')
    parser.add_argument('@s','@@syd',type=str,default='./input_files/traces/',help='Input directory for synthetic traces')
    parser.add_argument('@t','@@trd',type=str,default='./Uobs/',help='Input directory for grand truth traces')
    parser.add_argument('@f','@@fmt',type=str,default='h5',help='Database format')
    parser.add_argument('@n','@@nam',type=str,nargs='+',default=['all'],help = 'Name of the set(s) of monitors')
    parser.add_argument('@v','@@var',type=str,nargs='+',default=['Displ'],help='Output variables') # ['Displ','Veloc','Accel']
    parser.add_argument('@r','@@rdr',type=str,nargs='+',default=['x','y','z'],help='Motion components')
    parser.add_argument('@m','@@mon',type=int,nargs='+',default=[0],help='Monitor number')
    opt = parser.parse_args().__dict__
    
    
    print("Parse {} database ({}*.{})\nVariables: {}-Comp: {}".format(opt['nam'],opt['syd'],opt['fmt'],
                                                              opt['var'],opt['rdr']))
    
    # parse database
    opt_truth = opt.copy()
    opt_truth['wkd'] = opt['trd']
    opt_synthetic = opt.copy()
    opt_synthetic['wkd'] = opt['syd']

    stream_truth = parse(**opt_truth)
    print("Grand truth database parsed!\n")
    stream_synthetic = parse(**opt_synthetic)
    print("Synthetic database parsed!\n")

    # compute misift
    trace_truth = list(stream_truth.values())[0]
    trace_synthetic = list(stream_synthetic.values())[0]

    for direction in ['x','y','z']:
        for m in range(3):
            misfit = np.linalg.norm(trace_truth.displ_values(m, direction) - trace_synthetic.displ_values(m, direction))
            print(f'Misfit for {direction} component of monitor {m}: {misfit}')
            
            plt.plot(trace_synthetic.Time,trace_synthetic.displ_values(m, 'z'), label="Synthetic")
            plt.plot(trace_truth.Time,trace_truth.displ_values(m, 'z'), label="Grand truth")
            plt.legend()
            plt.savefig(f'measures_{direction}_{m}.png')
            plt.clf()
