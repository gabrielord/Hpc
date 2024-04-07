import numpy as np
from json import dumps, loads
import argparse

def recursion_two_loop(gradient_vec, s_stored, y_stored, m):
    q = gradient_vec
    a = np.zeros(m)
    rou = np.array([1/np.dot(y_stored[j, :], s_stored[j, :]) for j in range(m)])
    for i in range(m):
        a[m - 1 - i] = rou[m - 1 - i] * np.dot(s_stored[m - 1 - i, :], q)
        q = q - a[m - 1 - i]*y_stored[m - 1 - i, :]
    
    H_k0 = (np.dot(s_stored[m - 1], y_stored[m - 1])/np.dot(y_stored[m - 1], y_stored[m - 1]))
    r = H_k0 * q
    
    for i in range(m):
        beta = rou[i] * np.dot(y_stored[i, :], r)
        r = r + (a[i] - beta) * s_stored[i]
    return r

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prefix_chars='@')
    parser.add_argument('@@iter',type=int,help="iteration number")
    parser.add_argument('@@data_filename',type=str,help="data filename")
    opt = parser.parse_args().__dict__

    with open(opt['@@data_filename'], "r") as f:
        data = loads(f.read())

    s_mu = - recursion_two_loop(data["nabla_mu"], np.array(data["s_mu_stored"]), np.array(data["y_mu_stored"]), data["m"])
    s_lamb = - recursion_two_loop(data["nabla_lamb"], np.array(data["s_lamb_stored"]), np.array(data["y_lamb_stored"]), data["m"])

    with open(f"output_files/step_direction_{opt['@@iter']}.txt", "w") as f:
        dumps({"s_mu": s_mu.tolist(), "s_lambda": s_lamb.tolist()}, f)