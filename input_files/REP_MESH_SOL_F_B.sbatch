#!/bin/bash

# Définition des variables N_iter à 10 (nombre d'itérations), alpha_lamb, alpha_mu, c1, xi, et tol
export N_iter=10
export alpha_lamb=1
export alpha_mu=1
export c1=1e-4
export xi=0.5
export tol=1e-4

sbatch "MESHER_FORWARD.sbatch"
sbatch "SOLVER_FORWARD.sbatch"
J=$(sbatch "MISFIT.sbatch")

# Boucle while pour lancer plusieurs jobs sbatch
while (( $(bc <<< "$J > $tol") ))
do
    sbatch "MESHER_FORWARD.sbatch"
    sbatch "SOLVER_FORWARD.sbatch"
    J=$(sbatch "MISFIT.sbatch")
    sbatch "MESHER_BACKWARD.sbatch"
    sbatch "SOLVER_BACKWARD.sbatch"
    read g_lamb_k g_mu_k < <(sbatch "SNAPSHOTS.sbatch")
    s_lamb_k=$g_lamb_k
    s_mu_k=$g_mu_k
    J_k=$J
    sbatch "UPDATE_LAMB_MU_BACKTRACK.sbatch"
    sbatch "MESHER_FORWARD.sbatch"
    sbatch "SOLVER_FORWARD.sbatch"
    J_k_1=$(sbatch "MISFIT.sbatch")
    while (( $(bc <<< "$J_k_1 >= $J_k + $c1*($alpha_lamb*$g_lamb_k*$s_lamb_k+$alpha_mu*$g_mu_k*$s_mu_k)") ))
    do
        J_k=$J_k_1
        sbatch "UPDATE_LAMB_MU_BACKTRACK.sbatch"
        sbatch "MESHER_FORWARD.sbatch"
        sbatch "SOLVER_FORWARD.sbatch"
        J_k_1=$(sbatch "MISFIT.sbatch")
    done
done
