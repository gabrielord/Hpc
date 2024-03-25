
import numpy as np 
def gradient(mu,lamb,R_lamb,R_mu, g_reg_lamb,g_reg_mu, g_lamb_mis, g_mu_mis):

    dco_lamb = lamb.T * (R_lamb*g_reg_lamb + g_lamb_mis)
    dco_mu = mu.T * (R_mu*g_reg_mu + g_mu_mis)
    return dco_lamb, dco_mu


if __name__ == "__main__":
    import sys
    import numpy as np
    mu = np.array(sys.argv[1], dtype=float)
    lamb = np.array(sys.argv[2], dtype=float)
    R_lamb = np.array(sys.argv[3], dtype=float)
    R_mu = np.array(sys.argv[4], dtype=float)
    g_reg_lamb = np.array(sys.argv[5], dtype=float)
    g_reg_mu = np.array(sys.argv[6], dtype=float)
    g_lamb_mis = np.array(sys.argv[7], dtype=float)
    g_mu_mis = np.array(sys.argv[8], dtype=float)

    result = gradient(mu, lamb, R_lamb, R_mu, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis)

 
    with open("output_grad.txt", "w") as f:
        for array in result:
            f.write(' '.join(map(str, array)) + '\n')
