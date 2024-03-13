# -*- mode: perl -*-
run_name = "Test_name";

# duration of the run
sim_time = 2.0;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim=3;
fmax=0.1;
ngll=3;

snapshots {
    save_snap = true;
    # seconds of real simulation
    select material = 0;
    snap_interval = 0.5; 
};

# Mpml attenuation
mpml_atn_param = 0.02;

# Monitor structure
save_traces   = true;
traces_format = hdf5;
#traces_format = text;

capteurs "Uobs" {
    type   = points;
    file   = "stations.txt";
    period = 1; 
};
capteurs "Mesh1DLambdaMu" {
    type   = points;
    file   = "stations_z.txt";
    period = 1; 
};

# Fichier protection reprise
prorep=false;
prorep_iter=200000;
restart_iter=0;

# introduce a source for monitor 0 for direction x
source {
    # coordinates of the sources (x,y,z)
    coords = 0. 0. 0.;
    # Type (1.Impulse, 2.moment, 3.fluidpulse)
    type = impulse;
    # Direction 0.x,1.y ou 2.z (only for Impulse)
    dir = 1 0 0 ;
    # Source Time Function (STF)
    func = file;
    # Read STF from external file
    time_file = "monitors_misfit/monitor_0_x.txt";
    # Scaling amplitude
    amplitude = 1.;
};
# introduce a source for monitor 0 for direction y
source {
    # coordinates of the sources (x,y,z)
    coords = 0. 0. 0.;
    # Type (1.Impulse, 2.moment, 3.fluidpulse)
    type = impulse;
    # Direction 0.x,1.y ou 2.z (only for Impulse)
    dir = 0 1 0 ;
    # Source Time Function (STF)
    func = file;
    # Read STF from external file
    time_file = "monitors_misfit/monitor_0_y.txt";
    # Scaling amplitude
    amplitude = 1.;
};
# introduce a source for monitor 0 for direction z
source {
    # coordinates of the sources (x,y,z)
    coords = 0. 0. 0.;
    # Type (1.Impulse, 2.moment, 3.fluidpulse)
    type = impulse;
    # Direction 0.x,1.y ou 2.z (only for Impulse)
    dir = 0 0 1 ;
    # Source Time Function (STF)
    func = file;
    # Read STF from external file
    time_file = "monitors_misfit/monitor_0_z.txt";
    # Scaling amplitude
    amplitude = 1.;
};

# introduce a source for monitor 1 for direction x
source {
    # coordinates of the sources (x,y,z)
    coords = 0. 0. -1.;
    # Type (1.Impulse, 2.moment, 3.fluidpulse)
    type = impulse;
    # Direction 0.x,1.y ou 2.z (only for Impulse)
    dir = 1 0 0 ;
    # Source Time Function (STF)
    func = file;
    # Read STF from external file
    time_file = "monitors_misfit/monitor_1_x.txt";
    # Scaling amplitude
    amplitude = 1.;
};
# introduce a source for monitor 1 for direction y
source {
    # coordinates of the sources (x,y,z)
    coords = 0. 0. -1.;
    # Type (1.Impulse, 2.moment, 3.fluidpulse)
    type = impulse;
    # Direction 0.x,1.y ou 2.z (only for Impulse)
    dir = 0 1 0 ;
    # Source Time Function (STF)
    func = file;
    # Read STF from external file
    time_file = "monitors_misfit/monitor_1_y.txt";
    # Scaling amplitude
    amplitude = 1.;
};
# introduce a source for monitor 1 for direction z
source {
    # coordinates of the sources (x,y,z)
    coords = 0. 0. -1.;
    # Type (1.Impulse, 2.moment, 3.fluidpulse)
    type = impulse;
    # Direction 0.x,1.y ou 2.z (only for Impulse)
    dir = 0 0 1 ;
    # Source Time Function (STF)
    func = file;
    # Read STF from external file
    time_file = "monitors_misfit/monitor_1_z.txt";
    # Scaling amplitude
    amplitude = 1.;
};

# introduce a source for monitor 2 for direction x
source {
    # coordinates of the sources (x,y,z)
    coords = 0. 0. -2.;
    # Type (1.Impulse, 2.moment, 3.fluidpulse)
    type = impulse;
    # Direction 0.x,1.y ou 2.z (only for Impulse)
    dir = 1 0 0 ;
    # Source Time Function (STF)
    func = file;
    # Read STF from external file
    time_file = "monitors_misfit/monitor_2_x.txt";
    # Scaling amplitude
    amplitude = 1.;
};

# introduce a source for monitor 2 for direction y
source {
    # coordinates of the sources (x,y,z)
    coords = 0. 0. -2.;
    # Type (1.Impulse, 2.moment, 3.fluidpulse)
    type = impulse;
    # Direction 0.x,1.y ou 2.z (only for Impulse)
    dir = 0 1 0 ;
    # Source Time Function (STF)
    func = file;
    # Read STF from external file
    time_file = "monitors_misfit/monitor_2_y.txt";
    # Scaling amplitude
    amplitude = 1.;
};

# introduce a source for monitor 2 for direction z
source {
    # coordinates of the sources (x,y,z)
    coords = 0. 0. -2.;
    # Type (1.Impulse, 2.moment, 3.fluidpulse)
    type = impulse;
    # Direction 0.x,1.y ou 2.z (only for Impulse)
    dir = 0 0 1 ;
    # Source Time Function (STF)
    func = file;
    # Read STF from external file
    time_file = "monitors_misfit/monitor_2_z.txt";
    # Scaling amplitude
    amplitude = 1.;
};

time_scheme {
    accel_scheme = false;  # Acceleration scheme for Newmark
    veloc_scheme = true;   # Velocity scheme for Newmark
    alpha = 0.5;           # alpha (Newmark parameter)
    beta =  0.5;           # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
    courant=0.2;
};


pml_infos {
pml_type = CPML;
cpml_kappa0 = 1.0;
cpml_kappa1 = 0.0;
cpml_rc = 0.00001;
cpml_integration = Order2;
};

out_variables {
    enP  = 0;   # P-wave energy (scalar field)
    enS  = 0;   # S-wave energy (scalar field)
    evol = 1;   # volumetric strain (scalar field)
    pre  = 0;   # pressure (scalar field)
    dis  = 1;   # displacement (vector field)
    vel  = 0;   # velocity (vector field)
    acc  = 0;   # acceleration (vector field)
    edev = 1;   # deviatoric strain (tensor field)
    sdev = 0;   # deviatoric stress (tensor field)
    edevpl = 0; # deviatoric strain (tensor field)
    gradla = 1; # gradient Lambda
    gradmu = 1; # gradient Mu
};

