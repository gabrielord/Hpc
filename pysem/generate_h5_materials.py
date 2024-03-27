# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Script to create h5 material files for SEM3D

    Ex.1 generate linear gradient along z direction for a cube of [-100.0,100.0]x[-100.0,100.0]x[-100.0,100.0] 
    with step [20x20x20]
        
        python3 generate_h5_materials.py @@tag linear_gradient @@pfx example @@dir z @@xlim -100.0 100.0 @@ylim -100.0 100.0 @zlim -100.0 100.0
            @@step 20 20 20
"""
# Required modules
import argparse
import h5py
import numpy as np
import json
from scipy.interpolate import RegularGridInterpolator

# General informations
__author__ = "Filippo Gatti"
__copyright__ = "Copyright 2020, CentraleSupélec (MSSMat UMR CNRS 8579)"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"

dirs_dict = {'x':0,'y':1,'z':2}
trnsp = (2,1,0)

def heterogeneous_mat(la, mu, nodes):
    """
    Function to create the new material for lambda and mu by interpolating the values
    from la and mu on nodes to the grid points.
    """
    ### Definition of constants
    step = 30
    x_min, x_max = -1200, 1200
    y_min, y_max = -1200, 1200
    z_min, z_max = -300, 0
    nx, ny, nz = (x_max - x_min)//step+1, (y_max - y_min)//step+1, (z_max - z_min)//step+1
    lims = {'xmin':x_min, 'xmax':x_max, 'ymin':y_min, 'ymax':y_max, 'zmin':z_min, 'zmax':z_max, 'nx':nx, 'ny':ny, 'nz':nz}
    pfx = ''
    prop = ['la','mu']


    ### Handle SEM3D data
    # Get SEM3D nodes coordinates
    sem_x, sem_y, sem_z = zip(*nodes)
    sem_unique_x, sem_unique_y, sem_unique_z = np.unique(sem_x), np.unique(sem_y), np.unique(sem_z)
    sem_grid = np.meshgrid(sem_unique_x, sem_unique_y, sem_unique_z)
    # Reshape the SEM3D data
    la = la.reshape(sem_grid[0].shape)
    mu = mu.reshape(sem_grid[0].shape)

    ### Handle the grid coordinates
    # Create the grid
    grd = grid(lims)
    # Create list of points from the material grd
    points_grd = np.array([np.ravel(i) for i in grd]).T

    ### Calculate the new material
    # Interpolate the values of lambda and mu on the grid points
    mat_la = RegularGridInterpolator((sem_unique_x, sem_unique_y, sem_unique_z), la)(points_grd)
    mat_mu = RegularGridInterpolator((sem_unique_x, sem_unique_y, sem_unique_z), mu)(points_grd)
    # Reshape the interpolated values to the shape of the grid
    mat_la = mat_la.reshape(grd[0].shape).transpose(*trnsp)
    mat_mu = mat_mu.reshape(grd[0].shape).transpose(*trnsp)

    ### Generate the material
    # Definition of variable
    mat = {'la':mat_la,'mu':mat_mu}
    # Create the material
    write_h5(pfx,prop,mat,lims)
    write_xdmf(pfx,prop,mat,lims)

    print("Material files generated successfully!")


def base_smooth_heterogeneous(d,grd,nu=0.3):
    z = grd[d]
    la = np.zeros_like(grd[d])
    mu = np.zeros_like(grd[d])
    ds = np.full_like(grd[d],2000.).transpose(*trnsp)
    
    la = 20.*(80.0+0.45*np.abs(z)+\
        35.0*np.exp(-(np.abs(z)-22.5)**2/150.0))*1.e6 # Mpa
    # la = 20.*(100.0+0.45*np.abs(z)+\
#         50.0*np.exp(-(np.abs(z)-100.0)**2/1000.0))*1.e6 # Mpa
    mu = 0.5*(1.-2.*nu)*la/nu
    la = la.transpose(*trnsp)
    mu = mu.transpose(*trnsp)
    vp = np.sqrt((la+2.*mu)/ds)
    vs = np.sqrt(mu/ds)
    return {'la':la,'mu':mu,'ds':ds,'vp':vp,'vs':vs}

def linear_gradient(d,grd,nu=0.3):
    z = grd[d]
    la = np.zeros_like(grd[d])
    mu = np.zeros_like(grd[d])
    ds = np.full_like(grd[d],2000.).transpose(*trnsp)
    
    la = (100.0+0.45*np.abs(z))*1.e6 # Mpa
    mu = 0.5*(1.-2.*nu)*la/nu
    la = la.transpose(*trnsp)
    mu = mu.transpose(*trnsp)
    vp = np.sqrt((la+2.*mu)/ds)
    vs = np.sqrt(mu/ds)
    return {'la':la,'mu':mu,'ds':ds,'vp':vp,'vs':vs}

func = {'base_smooth_heterogeneous':base_smooth_heterogeneous,
        'linear_gradient':linear_gradient}

def grid(lims):
    xv = np.linspace(lims['xmin'],lims['xmax'],lims['nx'],dtype=np.float64)
    yv = np.linspace(lims['ymin'],lims['ymax'],lims['ny'],dtype=np.float64)
    zv = np.linspace(lims['zmin'],lims['zmax'],lims['nz'],dtype=np.float64)
    xg,yg,zg = np.meshgrid(xv,yv,zv,indexing='xy')
    return (xg,yg,zg)

def gen_mat(model,dirs,prop,grd,nu=0.3):
    d = dirs_dict[dirs]
    mats = func[model.lower()](d,grd,nu)
    return dict(tuple([(v,mats[v]) for v in prop]))

def write_h5(pfx,prop,mat,lims,xdmf=True):
    for v in prop:
        with h5py.File("{}_{}.h5".format(pfx,v),"w") as fid:
            samples = fid.create_dataset("samples",mat[v].shape,data=mat[v])
            # chunks=tuple([l//10 for l in mat[v].shape]))
            xMinGlob = np.array([lims['xmin'],lims['ymin'],lims['zmin']])
            xMaxGlob = np.array([lims['xmax'],lims['ymax'],lims['zmax']])
            fid.attrs['xMinGlob'] = xMinGlob
            fid.attrs['xMaxGlob'] = xMaxGlob
            fid.close()

def write_xdmf(pfx,prop,mat,lims):
    
    xMinGlob = np.array([lims['xmin'],lims['ymin'],lims['zmin']])
    xMaxGlob = np.array([lims['xmax'],lims['ymax'],lims['zmax']])
    dxV = (xMaxGlob-xMinGlob)/np.array([lims['nx']-1,lims['ny']-1,lims['nz']-1])
                               
    for v in prop:
        szs = mat[v].shape
        with open("{}_{}.xmf".format(pfx,v),"w") as fid:
            fnm = "{}_{}".format(pfx,v)
            fid.write('''<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n'''+
                      '''<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n'''+
                      '''<Domain>\n''')
            fid.write('   <DataItem Name="{}" Format="HDF" DataType="Float" Precision="8" Dimensions="{} {} {}">\n'.format(v,*szs))
            fid.write('        %s.h5:/samples\n'%fnm)
            fid.write('   </DataItem>\n')
            fid.write('  <Grid GridType="Collection" CollectionType="Spatial">\n')
            fid.write('   <Grid Name="Group1">\n')
            fid.write('     <Topology TopologyType="3DCoRectMesh" Dimensions="%u %u %u"/>\n'%szs)
            fid.write('     <Geometry GeometryType="ORIGIN_DXDYDZ">\n')
            fid.write('   <DataItem Name="origin" Format="XML" DataType="Float" Precision="8" Dimensions="3">\n')
            fid.write('         %30.10f\n'%lims['zmin'])
            fid.write('         %30.10f\n'%lims['ymin'])
            fid.write('         %30.10f\n'%lims['xmin'])
            fid.write('   </DataItem>\n')
            fid.write('   <DataItem Name="step" Format="XML" DataType="Float" Precision="8" Dimensions="3">\n')
            fid.write('         %30.10f\n'%dxV[2])
            fid.write('         %30.10f\n'%dxV[1])
            fid.write('         %30.10f\n'%dxV[0])
            fid.write('   </DataItem>\n')
            fid.write('     </Geometry>\n')
            fid.write('     <Attribute Name="%s" Center="Node" AttributeType="Scalar">\n'%v)
            fid.write('       <DataItem Reference="XML">\n')
            fid.write('         /Xdmf/Domain/DataItem[@Name="%s"]\n'%v)
            fid.write('       </DataItem>\n')
            fid.write('     </Attribute>\n')
            fid.write('   </Grid>\n')
            fid.write('  </Grid>\n')
            fid.write(' </Domain>\n')
            fid.write('</Xdmf>\n')
            fid.close()

if __name__=='__main__':
    step = 30
    x_min, x_max = -1200, 1200
    y_min, y_max = -1200, 1200
    z_min, z_max = -300, 0
    
    parser = argparse.ArgumentParser(prefix_chars='@')
    parser.add_argument('@@prop',type=str,nargs='+',default= ['la','mu'],help="list of properties to be generated")
    parser.add_argument('@@iter',type=int,help="iteration number")
    opt = parser.parse_args().__dict__

    with open(f"output_files/gradient_values_{iter}.txt", "r") as f:
        data = json.loads(f.readline().strip())
    
    heterogeneous_mat(data["lambda"],data["mu"],data["Nodes"])