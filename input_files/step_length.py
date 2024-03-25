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
from ..pysem.gradient import gradient

from .parse_h5_snapshots import GetSnapshots

def write_output(txt):
    """
    Function to write in the output file.
    """
    with open("output.txt", "a") as output_file:
        output_file.write(txt+"\n")

def main():
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
        
        #Update les nodes pour changer le matériaux
        MPI.Init()
        comm = MPI.COMM_WORLD               # Get communicator
        size = MPI.COMM_WORLD.Get_size()    # Get size of communicator
        rank = MPI.COMM_WORLD.Get_rank()    # Get the current rank
        hostname = MPI.Get_processor_name() # Get the hostname
    
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

        

if __name__=='__main__':
    main()