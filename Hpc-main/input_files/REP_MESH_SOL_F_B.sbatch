#!/bin/bash
#SBATCH --job-name=FULL_PROJECT_SEISM
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --output=output.txt
#SBATCH --error=error.txt 
#SBATCH --partition=cpu_prod

# Définition des variables N_iter à 10 (nombre d'itérations), alpha_lamb, alpha_mu, c1, xi, et tol
export N_iter=10
export alpha_lamb=1
export alpha_mu=1
export c1=1e-4
export xi=0.5
export tol=1e-4

# Premier sbatch
sbatch "MESHER_FORWARD.sbatch"
sbatch --dependency=afterok:$J "SOLVER_FORWARD.sbatch"
J=$(sbatch --parsable --dependency=afterok:$J "MISFIT.sbatch")

# Initialisation de l'itération
iter=0

# Boucle while tant que l'itération est inférieure à N_iter et que J est supérieur à tol
while (( $iter < $N_iter )) && (( $(bc <<< "$J > $tol") ))
do
    sbatch --dependency=afterok:$J "MESHER_FORWARD.sbatch"
    sbatch --dependency=afterok:$J "SOLVER_FORWARD.sbatch"
    J=$(sbatch --parsable --dependency=afterok:$J "MISFIT.sbatch")
    sbatch --dependency=afterok:$J "MESHER_BACKWARD.sbatch"
    sbatch --dependency=afterok:$J "SOLVER_BACKWARD.sbatch"
    read g_lamb_k g_mu_k < <(sbatch --parsable --dependency=afterok:$J "SNAPSHOTS.sbatch")
    s_lamb_k=$g_lamb_k
    s_mu_k=$g_mu_k
    J_k=$J
    sbatch --dependency=afterok:$J "UPDATE_LAMB_MU_BACKTRACK.sbatch"
    sbatch --dependency=afterok:$J "MESHER_FORWARD.sbatch"
    sbatch --dependency=afterok:$J "SOLVER_FORWARD.sbatch"
    J_k_1=$(sbatch --parsable --dependency=afterok:$J "MISFIT.sbatch")
    while (( $(bc <<< "$J_k_1 >= $J_k + $c1*($alpha_lamb*$g_lamb_k*$s_lamb_k+$alpha_mu*$g_mu_k*$s_mu_k)") ))
    do
        J_k=$J_k_1
        sbatch --dependency=afterok:$J "UPDATE_LAMB_MU_BACKTRACK.sbatch"
        sbatch --dependency=afterok:$J "MESHER_FORWARD.sbatch"
        sbatch --dependency=afterok:$J "SOLVER_FORWARD.sbatch"
        J_k_1=$(sbatch --parsable --dependency=afterok:$J "MISFIT.sbatch")
    done
    ((iter++))
done