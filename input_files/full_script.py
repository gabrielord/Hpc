"""
Main script to determine the Lamé parameters of the geological layer.
"""
import os
import sys
import numpy as np
import subprocess
from math import floor
from datetime import datetime
import json

### Constants
N_iter = 10
alpha_lamb = 1
alpha_mu = 1
c1 = 1e-4
xi = 0.5
tol = 1e-4



def write_output(txt):
    """
    Function to write in the output file.
    """
    with open("output.txt", "a") as output_file:
        output_file.write(txt+"\n")




def get_current_time():
    """
    Return the current date with format DD/MM/YYYY à HHhMMminSS
    """
    current_datetime = datetime.now()
    return current_datetime.strftime("%d/%m/%Y à %Hh%Mmin%S")




def create_useful_folders():
    """
    Create the folders that are missing.
    """
    folders = ["output_files", "monitors_misfit", "sem", "stats"]
    for folder in folders:
        os.makedirs(folder, exist_ok=True)
    
    sem3d_folders = ["traces", "res"]
    for folder in sem3d_folders:
        os.makedirs(f"forward/{folder}", exist_ok=True)
        os.makedirs(f"backward/{folder}", exist_ok=True)




def copy_material_files():
    """
    Copy the material files in the material folder to keep a version
    of the initial state of the layer.
    """
    os.makedirs("material", exist_ok=True)
    # On copie tous les fichiers de la forme example_*.h5 et example_*.xml
    # au format *.h5 et *.xml dans le dossier material
    os.system("""for file in material/initial_*; do cp "$file" "$(echo $file | sed 's/initial_//')"; done""")




def get_output(file_path):
    """
    Get the output of a function stored in a file with json dumps.
    """
    with open(file_path, "r") as f:
        ligne = f.readline().strip()
        return json.loads(ligne)




def main():
    """Main function of the script."""
    ### Initialization
    create_useful_folders()
    copy_material_files()

    ### Variables
    iter = 0
    J = None

    while iter < N_iter and (J is None or J > tol):
        write_output("\n" + 2*f"***********************{'*'*floor(np.log10(max(1, iter)))}\n"+ \
                              f"***** ITERATION {iter} *****\n" + \
                            2*f"***********************{'*'*floor(np.log10(max(1, iter)))}\n")
        if J is not None:
            write_output(f"Current cost function value: {J}")
        else:
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
            with open(f"output_files/forward_{iter}.solver", "w") as f:
                subprocess.run(cmd, shell=True, stdout=f)
            
            write_output("\t- Move the output files...")
            os.system("rm -r forward/*")
            os.system("mv traces forward")
            os.system("mv res forward")
            os.system(f"mv stat.log output_files/stat_forward_{iter}.log")



        ### We determine the misfit
        write_output(f"\n\n2. [{get_current_time()}] Determining the misfit...")
        # Generate the misfit files
        write_output(f"\n==> [{get_current_time()}] Generating the misfit files...")
        cmd = f"mpirun -np 1 " + \
              f"-map-by ppr:1:core:PE=1 " + \
              f"python3 ../pysem/compute_misfit.py"
        with open(f"output_files/misfit_{iter}.output", "w") as f:
            subprocess.run(cmd, shell=True, stdout=f)


        ### We solve the adjoint (backward) problem
        write_output(f"\n\n3. [{get_current_time()}] Solving the adjoint problem...")
        # Launch the mesher
        write_output(f"\n==> [{get_current_time()}] Launching the mesher...")
        write_output("\t- Update the input.spec file...")
        os.system("cp input_backward.spec input.spec")
        write_output("\t- Launch the mesher...")
        os.system(f"mesher < mesh.input > output_files/backward_{iter}.mesher")
        write_output("\t- Move the mesh files...")
        os.system("mv mesh4spec.* ./sem/")

        # Launch the solver
        write_output(f"\n==> [{get_current_time()}] Launching the solver...")
        write_output("\t- Launch the solver...")
        cmd = f"mpirun -np 32 " + \
                f"-map-by ppr:1:core:PE=1 " + \
                f"sem3d.exe"
        with open(f"output_files/backward_{iter}.solver", "w") as f:
            subprocess.run(cmd, shell=True, stdout=f)
            
        write_output("\t- Move the output files...")
        os.system("rm -r backward/*")
        os.system("mv traces backward")
        os.system("mv res backward")
        os.system(f"mv stat.log output_files/stat_backward_{iter}.log")

        ### We determine the gradient
        write_output(f"\n\n4. [{get_current_time()}] Determining the gradient...")
        # Calculate the gradient
        write_output(f"\n==> [{get_current_time()}] Calculating the gradient...")
        cmd = f"mpirun -np 1 " + \
                f"python3 ../pysem/parse_h5_snapshots.py @@iter {iter}"
        with open(f"output_files/snapshot_processing_{iter}.output", "w") as f:
            subprocess.run(cmd, shell=True, stdout=f)
        
        # Get the output values
        write_output(f"\n==> [{get_current_time()}] Getting the resulting output values...")
        data = get_output(f"output_files/gradient_values_{iter}.txt")
        Lambda = data["lambda"]
        Mu = data["mu"]
        g_lambda = data["g_lambda"]
        g_mu = data["g_mu"]
        M = data["M"]
        Nodes = data["nodes"]

        if J is None:
            # Compute the cost function
            write_output(f"\n==> [{get_current_time()}] Computing the cost function...")
            cmd = f"mpirun -np 1 " + \
                    f"python3 ../pysem/cost_function.py @@iter {iter}"
            with open(f"output_files/cost_computation_{iter}.output", "w") as f:
                subprocess.run(cmd, shell=True, stdout=f)
            J = get_output(f"output_files/cost_{iter}.txt")["cost"]


        # Calculate the search direction
        write_output(f"\n==> [{get_current_time()}] Calculating the search direction...")
        m = 10
        if iter == 0:
            y_mu_stored = []
            y_lamb_stored = []
            s_mu_stored = []
            s_lamb_stored = []
            nabla_mu = data["g_mu"] # initial gradient 
            nabla_lamb = data["g_lambda"]
        
        grad_mu_old = nabla_mu[:]
        grad_lamb_old = nabla_lamb[:]

        with open("output_files/data_for_search_direction.txt", "w") as f:
            f.write(json.dumps({"nabla_mu": list(nabla_mu), "nabla_lamb": list(nabla_lamb), "s_mu": list(s_mu_stored), "s_lamb": list(s_lamb_stored), "y_mu": list(y_mu_stored), "y_lamb": list(y_lamb_stored)}))

        cmd = f"mpirun -np 1 " + \
                f"python3 ../pysem/step_direction.py @@iter {iter} @@data_filename output_files/data_for_search_direction.txt"
        with open(f"output_files/step_direction_{iter}.output", "w") as f:
            subprocess.run(cmd, shell=True, stdout=f)
        
        ### We read the output values of s_mu and s_lambda
        with open(f"output_files/step_direction_{iter}.txt", "r") as f:
            data = json.loads(f.read())
            s_mu = np.array(data["s_mu"])
            s_lambda = np.array(data["s_lambda"])

        # Calculate the step lengths
        write_output(f"\n==> [{get_current_time()}] Calculating the step lengths...")
        newJ = None
        alpha_lamb, alpha_mu = 1, 1
        while newJ is None or newJ >= J + c1*(alpha_lamb*np.dot(s_lambda, g_lambda) + alpha_mu*np.dot(s_mu, g_mu)):
            if newJ is not None:
                alpha_lamb /= xi
                alpha_mu /= xi

            ### We update the material parameters
            write_output(f"\n\t--> [{get_current_time()}] Updating the material parameters...")
            cmd = f"mpirun -np 1 " + \
                  f"python3 ../pysem/generate_h5_materials.py @@iter {iter} @@alpha_lamb {alpha_lamb} @@alpha_mu {alpha_mu}"
            with open(f"output_files/materials_{iter}.output", "w") as f:
                subprocess.run(cmd, shell=True, stdout=f)
            
            ## We solve the forward problem
            write_output(f"\n\n1. [{get_current_time()}] Solving the forward problem...")
            # Launch the mesher
            write_output(f"\n==> [{get_current_time()}] Launching the mesher...")
            write_output("\t- Update the input.spec file...")
            os.system("cp input_forward.spec input.spec")
            write_output("\t- Launch the mesher...")
            os.system(f"mesher < mesh.input > output_files/forward_{iter+1}.mesher")
            write_output("\t- Move the mesh files...")
            os.system("mv mesh4spec.* ./sem/")
            
            # Launch the solver
            write_output(f"\n==> [{get_current_time()}] Launching the solver...")
            write_output("\t- Launch the solver...")
            cmd = f"mpirun -np 32 " + \
                f"-map-by ppr:1:core:PE=1 " + \
                f"sem3d.exe"
            with open(f"output_files/forward_{iter}.solver", "w") as f:
                subprocess.run(cmd, shell=True, stdout=f)
            write_output("\t- Move the output files...")
            os.system("rm -r forward/*")
            os.system("mv traces forward/traces")
            os.system("mv res forward/res")
            os.system(f"mv stat.log output_files/stat_forward_{iter+1}.log")

            ## We determinate the cost_function
            write_output(f"\n\n2. [{get_current_time()}] Determining the cost...")
            # Compute the cost function
            write_output(f"\n==> [{get_current_time()}] Computing the cost function...")
            cmd = f"mpirun -np 1 " + \
                  f"python3 ../pysem/cost_function.py @@iter {iter+1}"
            # We get the output values
            write_output(f"\n==> [{get_current_time()}] Getting the resulting output values...")
            newJ = get_output(f"output_files/cost_{iter+1}.txt")
        
        ### We update the variables for the next iteration of L-BFGS
        step_mu = - alpha_mu * nabla_mu
        step_lamb = - alpha_lamb * nabla_lamb

        s_mu_stored.append(step_mu)
        s_lamb_stored.append(step_lamb)

        ### We increment the iteration counter
        iter += 1
        J = newJ




if __name__ == "__main__":
    ### Call the main function
    main()