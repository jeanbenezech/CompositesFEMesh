import numpy as np
import os
from math import sqrt, pi, asin
from utils.parameters import parameters
from scipy.spatial.distance import pdist, cdist, squareform


def node_trans(point, param):
    '''
    Transform the nodes in flange and corner parts (C-spar is unfolded to be a 2D plane)
    '''
    
    local_xmax = param.Height;
    local_ymin = param.R;
    local_ymax = param.X+local_ymin;
    interior_radius = param.R;
    
    init = point
    ref = [[], [], []]
    moved = init.copy()      # nodes in the flange and corner parts will be transformed
    dist_from_bottom_surf = 0.0
    
    if init[1]>=local_ymax:
        if(init[0]<=local_xmax):
            dist_from_bottom_surf = (init[1] - local_ymax) - interior_radius;
        else:
            dist_from_bottom_surf = sqrt((init[1]-local_ymax)**2+(init[0]-local_xmax)**2)-interior_radius;
            
        ref[0] = local_xmax + interior_radius + dist_from_bottom_surf;
        ref[1] = param.X + interior_radius;
        moved[0] = ref[0];
        
        local_radius = interior_radius + dist_from_bottom_surf;

        if init[0]<=local_xmax:
            theta = pi/2.0;
            # moved[1] = ref[1] + param.R*theta + abs(init[0]-ref[0])-local_radius;
            moved[1] = ref[1] + param.R*theta + local_xmax-init[0];
        else:
            d = sqrt((init[1]-ref[1])**2+(init[0]-ref[0])**2);
            theta = 2*asin(d/2.0/local_radius);
            moved[1] = ref[1] + theta*param.R;
        
    elif init[1]<=local_ymin:
        if init[0]<=local_xmax:
            dist_from_bottom_surf = -init[1] + local_ymin - interior_radius;
        else:
            dist_from_bottom_surf = sqrt((init[1]-local_ymin)**2+(init[0]-local_xmax)**2)-interior_radius;
        
        ref[0] = local_xmax + interior_radius + dist_from_bottom_surf;
        ref[1] = interior_radius;
        moved[0] = ref[0];
        local_radius = interior_radius + dist_from_bottom_surf;

        if init[0]<=local_xmax:
            theta = pi/2.0;
            moved[1] = ref[1] - (param.R*theta + local_xmax-init[0]);
        else:
            d = sqrt((init[1]-ref[1])**2+(init[0]-ref[0])**2);
            theta = 2*asin(d/2.0/local_radius);
            moved[1] = ref[1] - theta*param.R;
        
    else:
        dist_from_bottom_surf = init[0] - local_xmax - interior_radius;
    
    if dist_from_bottom_surf>0.01:
        tmp = [[], []]  # 2D point in a plane
        tmp[0]=moved[1];
        tmp[1]=moved[2];
        return tmp
    else:
        return None
    
def get_nodes_trans():
    '''
    Get the transformed nodes in 2D plane
    '''
    # number of nodes: 2132
    num_nodes = 2132
    nodes = np.loadtxt(os.path.join('.','C-spar_continuum_shell.msh'), skiprows=5, max_rows=num_nodes, delimiter=' ')
    
    nodes_trans_index = [] # index of the transformed nodes
    nodes_trans = [] # transformed 2D nodes (y,z)
    
    param = parameters()
    param.init('parameters')
    for i in range(num_nodes):
        tmp = node_trans(nodes[i,1:], param)
        if tmp is not None:
            nodes_trans_index.append(i+1)
            nodes_trans.append(tmp)
    
    nodes_trans_index = np.array(nodes_trans_index)
    nodes_trans = np.array(nodes_trans)
    
    return nodes_trans_index, nodes_trans
    
def get_output_nodes():
    
    nodes_trans_index, nodes_trans = get_nodes_trans()
    
    #### target position of y and z in transformed 2D plane ####
    # z_targets = np.linspace(25, 475, 10)
    # y_targets = np.array((40, 100))
    
    # locations including the flange, corner and top surface
    z_targets = np.linspace(25, 475, 4)
    y_targets = np.array((-33, 1, 75, 149, 183))
    
    # zy_targets = []
    # zy_targets.append([75,50])
    # zy_targets.append([75,175])
    # zy_targets.append([75,225])
    # zy_targets.append([75,325])
    # zy_targets.append([75,450])
    
    node_indices = []
    node_coords = []
    
    # for zy_target in zy_targets:
        # z_target = zy_target[1]
        # y_target = zy_target[0]
    for z_target in z_targets:
        for y_target in y_targets:
            yz_target = np.array([[y_target, z_target]])
            dist = cdist(yz_target, nodes_trans)
            arg_target = dist.argmin()
            node_indices.append(nodes_trans_index[arg_target])
            node_coords.append(nodes_trans[arg_target].tolist())
    
    return node_indices, node_coords
    
    
    # # number of elements: 1800
    # # ELEMENT, type=SC8R, ELSET=Hexahedron
    # # each element consists of 8 nodes
    # elements = np.loadtxt('.\Abaqus\Sample_mesh.inp', skiprows=3782,               max_rows=1800, delimiter=',', comments='*')
    
def cov_mat_eval(node_coords, Sigma, Length):
    '''
    Calculate the covariance matrix of all the measurements
    '''
    node_coords = np.array(node_coords)
    dist2 = squareform(pdist(node_coords, 'sqeuclidean'))
    cov_mat = Sigma**2 * np.exp( - dist2 / (Length**2) );
    
    return cov_mat
    
if __name__ == '__main__':
    node_indices, node_coords = get_output_nodes()
    print(node_indices)
    print(node_coords)
    
    cov_mat = cov_mat_eval(node_coords, 1, 10)
