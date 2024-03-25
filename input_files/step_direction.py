'''
    snp_forward = GetSnapshots(comm,size,rank)

    Mtilde = snp_forward.dset['Mass']/snp_forward.dset['Dens']

    lambda_k, mu_k, s_lambda_k, s_mu_k, g_tilde_lambda_k, g_tilde_mu_k = get_properties()
    alpha_lambda, alpha_mu = calculate_alpha(lambda_k, mu_k, alpha_lambda, s_lambda_k, alpha_mu, s_mu_k, g_tilde_lambda_k, g_tilde_mu_k)
    dlambda_misfit, dmu_misfit = gradient(snp_forward, Mtilde, g_tilde_lambda_k, g_tilde_mu_k)
    lambda_k = lambda_k + alpha_lambda * dlambda_misfit
    mu_k = mu_k + alpha_mu * dmu_misfit'''

