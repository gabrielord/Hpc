#!/bin/bash
#SBATCH --job-name=FULL_PROJECT_SEISM
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --output=output.txt
#SBATCH --error=error.txt 
#SBATCH --partition=cpu_prod

# Définition des variables N_iter à 10 (nombre d'itérations), alpha_lamb, alpha_mu, c1, xi, et tol
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
    sbatch "MESHER_FORWARD.sbatch"
    wait $!
    sbatch "SOLVER_FORWARD.sbatch"
    wait $!
    sbatch "MISFIT.sbatch"
    J =$(cat misfit.txt)
    sbatch "MESHER_BACKWARD.sbatch"
    wait $!
    sbatch "SOLVER_BACKWARD.sbatch"
    wait $!
    sbatch "SNAPSHOTS.sbatch"
    wait $!
    (g_lamb_k,g_lamb_mu) = read g_lamb_mu.txt
    wait $!
    s_lamb_k=$g_lamb_k
    s_mu_k=$g_mu_k
    J_k=$J
    sbatch "UPDATE_LAMB_MU_BACK
    wait $!
    sbatch "MESHER_FORWARD.sbatch"
    wait $!
    sbatch "SOLVER_FORWARD.sbatch"
    wait $!
    sbatch "MISFIT.sbatch"
    J_k_1 = $(cat misfit.txt)
    while [ "$J_k_1 -lt $J_k + $c1*($alpha_lamb*$g_lamb_k*$s_lamb_k+$alpha_mu*$g_mu_k*$s_mu_k)" ]
    do
        J_k=$J_k_1
        sbatch "UPDATE_LAMB_MU_BACKTRACK.sbatch"
        wait $!
        sbatch "MESHER_FORWARD.sbatch"
        wait $!
        sbatch "SOLVER_FORWARD.sbatch"
        wait $!
        J_k_1=sbatch "MISFIT.sbatch"
        J_k_1 = $(cat misfit.txt)
    done
    ((iter++))
done