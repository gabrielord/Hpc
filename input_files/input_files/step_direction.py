'''
    snp_forward = GetSnapshots(comm,size,rank)

    Mtilde = snp_forward.dset['Mass']/snp_forward.dset['Dens']

    lambda_k, mu_k, s_lambda_k, s_mu_k, g_tilde_lambda_k, g_tilde_mu_k = get_properties()
    alpha_lambda, alpha_mu = calculate_alpha(lambda_k, mu_k, alpha_lambda, s_lambda_k, alpha_mu, s_mu_k, g_tilde_lambda_k, g_tilde_mu_k)
    dlambda_misfit, dmu_misfit = gradient(snp_forward, Mtilde, g_tilde_lambda_k, g_tilde_mu_k)
    lambda_k = lambda_k + alpha_lambda * dlambda_misfit
    mu_k = mu_k + alpha_mu * dmu_misfit'''

import os
import sys
import numpy as np
import subprocess
from math import floor
from datetime import datetime
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI
from ..pysem.gradient import gradient_mu, gradient_lamb
from full_script import get_current_time

from ..pysem.parse_h5_snapshots import GetSnapshots

def write_output(txt):
    """
    Function to write in the output file.
    """
    with open("output.txt", "a") as output_file:
        output_file.write(txt+"\n")

def step_length(s_lamb, g_lamb, s_mu, g_mu, comm, rank, size):
    #Init varibale
    N_iter = 5
    alpha_lamb = 1
    alpha_mu = 1
    c1= 10**(-4)
    xi = 0.5


    with open('cout.txt', 'r') as file:
        J_k_1 = file.readlines()[-1].strip()

    iter = 0 
    #On initialise J_k à None, pour, à la fin de la permière boucle, avoir J_k à la valeur calculé, et J_k_1 à la nouvelle valeur
    J_k = None

    #Backtracking
    while iter < N_iter and (J_k is None or J_k_1 > J_k + c1*(alpha_lamb*np.dot(s_lamb,g_lamb)+ alpha_mu*np.dot(s_mu,g_mu)) ):
        write_output("\n" + 2*f"***********************{'*'*floor(np.log10(max(1, iter)))}\n"+ \
                              f"***** ITERATION {iter} *****\n" + \
                            2*f"***********************{'*'*floor(np.log10(max(1, iter)))}\n")
        if J_k_1 is not None:
            write_output(f"Current cost function value, J_k_1 and J_k: {(J_k_1,J_k)}")
        alpha_lamb *= xi
        alpha_mu *= xi
    
        snp_forward, snp_backward = GetSnapshots(comm,size,rank)
        nodes = snp_forward.dset['Nodes']
        MPI.Finalize()

        ### We update the material parameters
        write_output(f"\n\n5. [{get_current_time()}] Updating the material parameters...")

        ###Change parametes in material file 
        write_output(f"\n\n5. [{get_current_time()}] Updating the material parameters...")

        ### We solve the forward problem
        write_output(f"\n\n1. [{get_current_time()}] Solving the forward problem...")

        # Launch the mesher
        write_output(f"\n==> [{get_current_time()}] Launching the mesher...")
        write_output("\t- Update the input.spec file...")
        os.system("cp input_forward.spec input.spec")
        write_output("\t- Launch the mesher...")
        os.system(f"mesher < mesh.input > output_files/forward_{iter}.mesher")
        write_output("\t- Move the mesh files...")
        os.system("mv mesh4spec.* ./sem/")

        # Launch the solver
        write_output(f"\n==> [{get_current_time()}] Launching the solver...")
        write_output("\t- Launch the solver...")
        cmd = f"mpirun -np 32 " + \
              f"-map-by ppr:1:core:PE=1 " + \
              f"sem3d.exe"
        f = open(f"output_files/forward_{iter}.solver", "w")
        subprocess.run(cmd, shell=True, stdout=f)
        f.close()
        write_output("\t- Move the output files...")
        os.system("mv traces forward/traces")
        os.system("mv res forward/res")
        os.system(f"mv stat.log output_files/stat_forward_{iter}.log")

        ### We determine the cost_function
        write_output(f"\n\n2. [{get_current_time()}] Determining the cost...")
        # Generate the misfit files
        write_output(f"\n==> [{get_current_time()}] Generating the cost files...")
        cmd = f"mpirun -np 1 " + \
              f"-map-by ppr:1:core:PE=1 " + \
              f"python3 ../pysem/costFonction.py"
        f = open(f"output_files/misfit_{iter}.output", "w")
        subprocess.run(cmd, shell=True, stdout=f)
        f.close()

        J_k = J_k_1
        with open('cout.txt', 'r') as file:
            J_k_1 = file.readlines()[-1].strip()
        iter +=10
    return alpha_lamb, alpha_mu


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

### calculate an initial nabla_mu and nabla_lamb

def L_bfgs(mu, lamb, y_mu_stored, y_lamb_stored, s_mu_stored, s_lamb_stored, nabla_mu, nabla_lamb, m, g_lamb, g_mu, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis, comm, rank, size):    
    '''
    Store the {y_i, s_i}
    '''
    grad_mu_old = nabla_mu[:]
    grad_lamb_old = nabla_lamb[:]

    nabla_mu = gradient_mu(mu, lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis) # initial gradient 
    nabla_lamb = gradient_lamb(mu, lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis)

    s_mu = - recursion_two_loop(nabla_mu, np.array(s_mu_stored), np.array(y_mu_stored), m)
    s_lamb = - recursion_two_loop(nabla_lamb, np.array(s_lamb_stored), np.array(y_lamb_stored), m)

    alpha_lamb, alpha_mu = step_length(s_lamb, g_lamb, s_mu, g_mu, comm, rank, size)

    # p = - nabla
    step_mu = - alpha_mu * nabla_mu
    step_lamb = - alpha_lamb * nabla_lamb

    s_mu_stored.append(step_mu)
    s_lamb_stored.append(step_lamb)
    
    # nabla_mu = gradient_mu(mu,lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis) # initial gradient 
    # nabla_lamb = gradient_lamb(mu,lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis)
    
    y_mu_stored.append(nabla_mu - grad_mu_old)
    y_lamb_stored.append(nabla_lamb - grad_lamb_old)

    mu = mu + step_mu
    lamb = lamb + step_lamb

    # s_mu = - recursion_two_loop(nabla_mu, np.array(s_mu_stored), np.array(y_mu_stored), m)
    # s_lamb = - recursion_two_loop(nabla_lamb, np.array(s_lamb_stored), np.array(y_lamb_stored), m)

    return mu,lamb
'''
    while np.linalg.norm(nabla_mu) > 1e-5 and np.linalg.norm(nabla_lamb) > 1e-5: # while gradient is positive
        if it > max_it: 
            print('Maximum iterations reached!')
            break

        if 0 < it and it < m :
            p_mu = - recursion_two_loop(nabla_mu, np.array(s_mu_stored), np.array(y_mu_stored), m_)
            p_lamb = - recursion_two_loop(nabla_mu, np.array(s_lamb_stored), np.array(y_lamb_stored), m_)
            
            s_mu_stored.append(alpha_mu * p_mu)
            s_lamb_stored.append(alpha_lamb * p_lamb)
            
            grad_mu_old = nabla_mu[:]
            grad_lamb_old = nabla_lamb[:]
            
            mu = mu + alpha_mu * p_mu
            lamb = lamb + alpha_lamb * p_lamb
            
            nabla_mu = gradient_mu(mu,lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis) 
            nabla_lamb = gradient_lamb(mu,lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis)

            y_mu_stored.append(nabla_mu - grad_mu_old)
            y_lamb_stored.append(nabla_lamb - grad_lamb_old)

            m_ = m_ + 1
            it = it + 1
            
        else:
            p_mu = - recursion_two_loop(nabla_mu, np.array(s_mu_stored), np.array(y_mu_stored), m)
            p_lamb = - recursion_two_loop(nabla_mu, np.array(s_lamb_stored), np.array(y_lamb_stored), m)
            
            # alpha_lamb, alpha_mu = step_length(s_lamb, g_lamb, s_mu, g_mu, comm, rank, size)

            #append the s_k+1 
            s_mu_stored.append(alpha_mu * p_mu)
            s_lamb_stored.append(alpha_lamb * p_lamb)

            #discard the s_(k-m)
            s_mu_stored.pop(0)
            s_lamb_stored.pop(0)
            
            grad_mu_old = nabla_mu[:]
            grad_lamb_old = nabla_lamb[:]

            mu = mu + alpha_mu * p_mu
            lamb = lamb + alpha_lamb * p_lamb
            
            nabla_mu = gradient_mu(mu,lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis) 
            nabla_lamb = gradient_lamb(mu,lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis)

            #append the y_k+1
            y_mu_stored.append(nabla_mu - grad_mu_old)
            y_lamb_stored.append(nabla_lamb - grad_lamb_old)

            #discard the y_k-m
            y_mu_stored.pop(0)
            y_lamb_stored.pop(0)
            
            it = it + 1
'''    
    

def main(mu, lamb, m, g_lamb, g_mu, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis, comm, rank, size):
    # Update les nodes pour changer le matériaux
    MPI.Init()
    comm = MPI.COMM_WORLD               # Get communicator
    size = MPI.COMM_WORLD.Get_size()    # Get size of communicator
    rank = MPI.COMM_WORLD.Get_rank()    # Get the current rank
    hostname = MPI.Get_processor_name() # Get the hostname

    y_mu_stored = []
    y_lamb_stored = []
    s_mu_stored = []
    s_lamb_stored = []

    nabla_mu_initial = gradient_mu(mu, lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis) # initial gradient 
    nabla_lamb_initial = gradient_lamb(mu, lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis)

    mu, lamb = L_bfgs(mu, lamb, y_mu_stored, y_lamb_stored, s_mu_stored, s_lamb_stored, nabla_mu_initial, nabla_lamb_initial, m, g_lamb, g_mu, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis, comm, rank, size)


if __name__=='__main__':
    main()