def gradient(mu,lamb,R_lamb,R_mu, g_reg_lamb,g_reg_mu, g_lamb_mis, g_mu_mis):

    dco_lamb = lamb.T * (R_lamb*g_reg_lamb + g_lamb_mis)
    dco_mu = mu.T * (R_mu*g_reg_mu + g_mu_mis)
    return dco_lamb, dco_mu

