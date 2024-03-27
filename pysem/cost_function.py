import numpy as np
import argparse
from parse_h5_traces import parse
from scipy.interpolate import interp1d
import json

def integrate_on_volume(var, x_values, y_values, z_values):
    grad_var = np.gradient(var, x_values, y_values, z_values)
    dx, dy, dz = np.diff(x_values), np.diff(y_values), np.diff(z_values)
    # On ajoute le dernier élément pour que les dimensions soient les mêmes
    dx = np.append(dx, dx[-1])
    dy = np.append(dy, dy[-1])
    dz = np.append(dz, dz[-1])
    # On calcule l'intégrale sur le volume de la norme du gradient
    grad = 0
    for i, values in enumerate([dx, dy, dz]):
        shape = [1] * grad_var[i].ndim
        shape[i] = -1
        values = values.reshape(shape)
        grad += np.sum((grad_var[i]*values)**2)
    return grad
    

def costFunction(iter, mu, lamb, R_lamb, R_mu, x_values, y_values, z_values):
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
    
    totalMisfit = np.sum(misfit)

    res_grad_lamb = integrate_on_volume(lamb, x_values, y_values, z_values)
    res_grad_mu = integrate_on_volume(mu, x_values, y_values, z_values)

    cost = 1/2*(totalMisfit + R_lamb*res_grad_lamb + R_mu*res_grad_mu)

    with open(f'output_files/cost_{iter}.txt', 'w') as f:
        f.write(json.dumps({"cost":cost}))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(prefix_chars='@')
    parser.add_argument('@@iter',type=int,help="iteration number")

    # We get the data from the gradient_values file
    with open(f"output_files/gradient_values_{iter}.txt", "r") as f:
        data = json.loads(f.readline().strip())
        mu = data['mu']
        lamb = data['lamb']
        nodes = data['nodes']
    iter = parser.parse_args().__dict__['iter']

    # We process mu and lambda to get a 3D array
    x_values = np.sort(np.unique(nodes[:, 0]))
    y_values = np.sort(np.unique(nodes[:, 1]))
    z_values = np.sort(np.unique(nodes[:, 2]))

    lamb_3d = np.zeros((len(x_values), len(y_values), len(z_values)))
    mu_3d = np.zeros((len(x_values), len(y_values), len(z_values)))

    for index, node in enumerate(nodes):
        x, y, z = node
        i, j, k = np.where(x_values == x)[0][0], np.where(y_values == y)[0][0], np.where(z_values == z)[0][0]
        lamb_3d[i, j, k] = lamb[index]
        mu_3d[i, j, k] = mu[index]

    cost = costFunction(iter, mu_3d, lamb_3d, 1, 1, x_values, y_values, z_values)

    print("Fonction de coût calculée avec succès !")