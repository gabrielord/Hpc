import numpy as np 

def gradient_lamb(mu,lamb, g_reg_lamb,g_reg_mu, g_lamb_mis, g_mu_mis):
    R_lamb = 3
    R_mu = 4
    dco_lamb = lamb.T @ (R_lamb * g_reg_lamb + g_lamb_mis)
    dco_mu = mu.T @ (R_mu * g_reg_mu + g_mu_mis)

    np.savetxt("output_grad.txt", (dco_lamb, dco_mu))
    
    return dco_lamb

def gradient_mu(mu,lamb, g_reg_lamb,g_reg_mu, g_lamb_mis, g_mu_mis):
    R_lamb = 3
    R_mu = 4
    dco_lamb = lamb.T @ (R_lamb * g_reg_lamb + g_lamb_mis)
    dco_mu = mu.T @ (R_mu * g_reg_mu + g_mu_mis)

    np.savetxt("output_grad.txt", (dco_lamb, dco_mu))
    
    return dco_mu