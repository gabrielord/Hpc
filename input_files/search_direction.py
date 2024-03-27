import numpy as np

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