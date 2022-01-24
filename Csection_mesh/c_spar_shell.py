import numpy as np
import os
import subprocess
import time
import psutil
#import pdb

from utils.surfaces import *
from utils.parameters import *
from utils.geometry import *

from write_shell_parameters import write_parameters
from write_shell_inp import write_inp_parameters

import node_element_extraction

class c_spar_shell_model():
    
    def __init__(self, gaussian_thickness=1, corner_thickness=1, output_node_indices=None, output_node_coords_2d=None):
        # gaussian_thickness  0: thickness has no randomness, 1: thickness is a gaussian random field
        self.gaussian_thickness = gaussian_thickness
        
        # corner_thickness  0: no thickness reduction at corner,  1: with thickness reduction at corner
        self.corner_thickness = corner_thickness
        
        # node indices and coords(2d) for model output
        if output_node_indices is None or output_node_coords_2d is None or self.gaussian_thickness == 1:
            self.write_geometry()
            self.get_mesh_gmsh()
            # transformed nodes in 2D plane (half of the total nodes)
            nodes_trans_index, nodes_trans = node_element_extraction.get_nodes_trans()
            
        if output_node_indices is None:
            self.output_node_indices = nodes_trans_index
        else:
            self.output_node_indices = output_node_indices
        
        if output_node_coords_2d is None:
            self.output_node_coords_2d = nodes_trans
        else:
            self.output_node_coords_2d = output_node_coords_2d
        
        if self.gaussian_thickness == 1:
            self.white_noise_size = len(nodes_trans_index)
            
        # random number generator
        self.rng = np.random.default_rng()
        

    def write_geometry(self, Ylength=7.0, ThicknessVar=2.0, Sigma=0.1, Length=10.0):
        '''
        write geometry parameters to 'parameters.txt' 
        '''
        write_parameters(Ylength=Ylength, ThicknessVar=ThicknessVar, Sigma=Sigma, Length=Length, is_GaussianThickness=self.gaussian_thickness, is_CornerThickness=self.corner_thickness)
    
    def get_mesh_gmsh(self):
        '''
        use gmsh to get the mesh files '.msh'
        '''
        filename_param = 'parameters'
        self.param = parameters()
        self.param.init(filename_param)

        geo = geometry()
        geo.init(self.param)
        geo.delta_p, geo.delta_l, geo.delta_s, geo.delta_v = front_surface(self.param, geo)
        self.param.write_input_orientation('input.txt')
        geo.dump(self.param)
        
        command = []
        command.append('gmsh')
        command.append(self.param.name+'.geo')
        command.append('-3')
        command.append('-format')
        command.append('msh2')
        command.append('-o')
        command.append(self.param.name+'.msh')
        
        out_gmsh = subprocess.run(command, capture_output=True, check=True)
    
    def get_mesh_ori_inp(self):
        '''
        modify the '.msh' file and get the inp file (mesh, oritation) for Abaqus 
        '''
        out_gridMod = subprocess.run(os.path.join('.','gridMod'), capture_output=True, check=True)
    
    def white_noise_to_file(self, white_noise=None):
        '''
        write white noise into file which will be use by 'gridMod'
        '''
        if white_noise is None:
            white_noise = self.rng.standard_normal(self.white_noise_size)
        
        np.savetxt('white_noise.txt', white_noise)
    
    def model_run(self, Ylength=7.0, ThicknessVar=2.0, Sigma=0.1, Length=10.0, white_noise=None):
        '''
        run c-spar shell model by Abaqus and extract the results
        '''
        # write geometry parameters to 'parameters.txt'
        self.write_geometry(Ylength=Ylength, ThicknessVar=ThicknessVar, Sigma=Sigma, Length=Length)
        
        # get the mesh file with gmsh
        self.get_mesh_gmsh()
        
        # write white noise to file if 'gaussian_thickness' is 1
        if self.gaussian_thickness == 1:
            self.white_noise_to_file(white_noise=white_noise)
        
        # modify mesh file and get the inp file (mesh, oritation) for Abaqus
        self.get_mesh_ori_inp()
        
        # write inp file (shell, material, boundary condition, step, ...) for Abaqus 
        write_inp_parameters(Ylength=Ylength)

        # change directory to folder 'Abaqus'
        os.chdir('Abaqus')
        
        # run Abaqus
        command_abq = []
        command_abq.append('abaqus')
        command_abq.append('job=sinan')
        command_abq.append('input='+self.param.name+'.inp')
        command_abq.append('cpus=1')
        command_abq.append('double=BOTH')
        # command_abq.append('int')
        out_abq = subprocess.run(command_abq, capture_output=True, check=True)
        
        
        # wait until the "odb" file exists and big enough 
        while os.path.isfile('sinan.odb') == False:
            time.sleep(0.1)
        
        start = time.time()
        while os.path.getsize('sinan.odb') < 2e6:
            time.sleep(0.5)
            if time.time() - start > 20: # time is too long, break
                break
        # print(os.path.getsize('sinan.odb'))
        # out_abq_t = subprocess.run('abaqus terminate job=sinan', shell=True, capture_output=True)
        
        
        # # get the strain and stress into a file
        # command_abq_1 = []
        # command_abq_1.append('abaqus')
        # command_abq_1.append('python')
        # command_abq_1.append('GetCellFields_for_Sinan.py')
        # command_abq_1.append('sinan.odb')
        # out_abq_1 = subprocess.run(command_abq_1, capture_output=True, check=True)
        
        # get the displacement into a file
        command_abq_2 = []
        command_abq_2.append('abaqus')
        command_abq_2.append('python')
        command_abq_2.append('Getdisplacement.py')
        command_abq_2.append('sinan.odb')
        
        out_abq_2 = subprocess.run(command_abq_2, capture_output=True, check=False)
        
        # read the file to get the output of interest.
        # if out_abq_1.returncode == 0:
        #     out_strain = np.loadtxt('sinan_SE.txt', skiprows=1, usecols=(5,6,7,8))
        #     out_stress = np.loadtxt('sinan_SE.txt', skiprows=1, usecols=(1,2,3,4))
        
        if out_abq_2.returncode == 0: # get displacement correctly
            out_disp = np.loadtxt('sinan_U.txt', skiprows=1, usecols=(1,2,3))
        
        # remove the Abaqus files
        out_rm = subprocess.run('rm -r sinan*', shell=True, capture_output=True)
        
        # kill process 'standard' 'pre' if they are still running
        for proc in psutil.process_iter():
            if proc.name() in ['standard', 'pre']:
                try:
                    proc.kill()
                except:
                    print('Process killing does not work')
        
        # change the directory back to the main directory
        os.chdir('..')
            
        if out_abq_2.returncode == 0:
            self.out_disp = out_disp
            return out_disp[np.array(self.output_node_indices)-1,-1]
        else:
            self.out_disp = None
            return np.ones(len(self.output_node_indices))*1e5
        
        
    # def c_spar_shell_model_derivative():

if __name__ == '__main__':
    
    # c_spar_shell_model()  # get initial '.msh' file for node extraction
    # node_indices, node_coords = node_element_extraction.get_output_nodes()
    model_c_spar = c_spar_shell_model(gaussian_thickness=1, corner_thickness=1, output_node_indices=None, output_node_coords_2d=None)
    
    start = time.time()
    for i in range(1):
        out_disp = model_c_spar.model_run(Ylength=7.0, ThicknessVar=2.0, Sigma=1.0, Length=10.0)
        print(time.time()-start)
        print(out_disp)
    
    # # generate data
    # data_disp_select_z = out_disp + 0.01*np.random.randn(30, out_disp.size)
    # np.savetxt(os.path.join('.','data','disp_select_z_full.csv'), out_disp)
    # np.savetxt(os.path.join('.','data','data_disp_select_z_full.csv'), data_disp_select_z)
