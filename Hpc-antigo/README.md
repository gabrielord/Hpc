# SEM3D

The following SEM3D test case SEM3D_ST7 contains the input file to create a three-layered geological domain, with the top layer corresponding to a heterogeneous soil layer (whose profile must be discovered by applying the adjoint problem) and the two deepest homogeneous layers (bedrock).

List of files:
	
	+ `mat.dat` : three layers, no PML. The top layer (top surface: 0 m; bottom surface: -300 m) is the object of the inversion problem, since its velocity profile (therefore, Lambda and Mu parameters) are unknown to the students.
	
	+ `mater.in` : contains the standard homogeneous characteristics for each layer. It is necessary to run the mesher, and produce `material.input`
	
	+ `material.spec` : to be used to run the SEM3D propagation solver on heterogeneous domain (only the top layer; corresponding to material 0 in `material.input`) 

	+ `stations.txt` : it contains the list of monitor's coordinates (x y z)

	+ `gaussian_stf.txt` : source time function according to Gaussian function and stored as two column file (time-displacement)

	+ `input.spec` : it contains all the informations necessary to run the wave propagation solver and record the traces (referring to `stations.txt` and creating a set of monitor called `Uobs` in `traces` folder) and save the snapshots (in `res` folder).

	+ `MESHER.sbatch` : slurm batch file to run the mesher 
	
	+ `SOLVER.sbatch` : slurm batch file to run the wave propagation solver

# PYSEM

## Preliminaries

	1. Connect to ruche01:

	```
	ssh -X your_login@ruche01.centralesupelec.fr
	```

	2. Connect to interactive mode (optional if you don't use a batch file)

	```
	srun --nodes=1 --time=00:30:00 -p cpu_short --pty /bin/bash
	```

	3. Create anaconda3 environment (if you don't have one already: this one is called invsem0X-st7)

	```
	module load anaconda3/2020.02/gcc-9.2.0
	export thisuser=$(whoami)
	export hmd="/gpfs/users"
	export wkd="/gpfs/workdir"
	
	mkdir -p ${wkd}/${thisuser}/.conda
	conda create -n ${thisuser}-st7 -y
	mkdir -p ${wkd}/${thisuser}/.conda
	source activate ${thisuser}-st7
	ln -s ${wkd}/${thisuser}/.conda ~/.conda

	conda install numpy matplotlib h5py scipy openmpi 
	conda install -c conda-forge mpi4py
	```

	4. Activate the anaconda3 environment (called ${thisuser}-st7)

	```
	source activate ${thisuser}-st7
	```

	5. Locate the `pysem` folder (referenced hereafter as `/path/to/pysem`) and the traces folder (referenced hereafter as `/path/to/traces`)

## Generate h5 material files for `material.spec`

Script to create hdf5 material files for SEM3D, corresponding to Lambda and Mu parameters following a linear gradient variation along z direction for a cube of [-100.0,100.0]x[-100.0,100.0]x[-100.0,100.0] with step [20x20x20]

```
python3 generate_h5_materials.py @@tag linear_gradient @@pfx example @@dir z @@xlim -100.0 100.0 @@ylim -100.0 100.0 @zlim -100.0 100.0 @@step 20 20 20
```

Be careful! The limits of your domain must be larger than the mesh domain considered.

## Parse h5 traces files (called `Uobs` in `input.spec`):

Parse hdf5 Displ (displacement traces) for all directions x y z, monitor 0 and 1 (in according to the list in `stations.txt`) and plot them in the current folder (that can be found by typing `pwd`)
        
```
python3 /path/to/pysem parse_h5_traces.py @@wkd /path/to/traces/ @@fmt h5 @@nam Uobs @@var Displ @@rdr x y z @@mon 0 1 2 @@plt
```