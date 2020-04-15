import os
import numpy as np

from bmtk.builder.networks import NetworkBuilder

# helper functions


def generate_positions(N, x0=0.0, x1=300.0, y0=0.0, y1=100.0):
    X = np.random.uniform(x0, x1, N)
    Y = np.random.uniform(y0, y1, N)
    return np.column_stack((X, Y))

# returns the stack in 3D space
def generate_positions_3D(N, x0=0.0, x1=300.0, z0=0.0, z1=300, y0=0.0, y1=300.0):
    X = np.random.uniform(x0, x1, N)
    Y = np.random.uniform(y0, y1, N)
    Z = np.random.uniform(z0, z1, N)
    result=np.column_stack((X, Y))
    return np.column_stack((result, Z))

# this function generates the positions of cells on the arch in 3D space
# phi and rho - are two input arrays that determine the polar coordinates of cells
# Y coordinates are equal to 0, so all cells are on the same plane
def generate_arc_positions_3D(phi, rho):
    factor=0.5
    X = np.multiply(rho,np.cos(phi*factor))
    Z = np.multiply(rho,np.sin(phi*factor))
    Y = np.zeros(len(phi))
    result=np.column_stack((X, Y))
    return np.column_stack((result, Z))


def generate_rotation_angle_1D(rho, angle):
  # rho is needed only for lenght
    n_cells = len(rho)
    population_angles = np.ones(int(n_cells))*angle

    return population_angles



# This function randomly picks up the spike trains from the external input and drives the
# GC and BS populations
def select_random_srcs(sources_lst, trg_cell, nconns=500, nsyns_min=2, nsyns_max=5):
  syns_selected = np.zeros(len(sources_lst), dtype=np.uint)
  select_syns = np.random.choice(len(sources_lst), nconns, replace=False)
  # np.random.randint(nsyns_min,nsyns_max,len(sources_lst),dtype='int') 
  #syns_selected[select_syns] = [nsyns]*len(sources_lst)
  syns_selected[select_syns] = np.random.randint(nsyns_min,nsyns_max,len(sources_lst),dtype='int')
  return syns_selected



def connection_phi_layer(phi_1,phi,n):
# this function takes the phi as an argument then determines if cells are connected or not
# function returns the list of coordinates of connected cells
# connections are in the same layer
# phi_1 - cell polar coordinate, phi - all phi coordinates available, n - number of connections
# returns the list of connected cells to the phi_1
  # make another array for shifted distances
  phi_shifted=np.array(phi)
  for i in np.arange(0,len(phi)):
    if phi[i]>np.pi:
      phi_shifted[i]=phi[i]-2*np.pi
    else:
      phi_shifted[i]=phi[i]
  # compute the miminum distance
  phi_dist = abs(phi_shifted - phi_1)
  #min_element_idx=np.argmin(phi_dist)
  min_element_idx=phi_dist.argsort()[:2*n]
  # take all elements, including itself
  min_element_idx=min_element_idx[:]
  # get the indexes of minimal elements
  phi_neighbours=phi[min_element_idx]
  # convert the np.array to list
  phi_neighbours=phi_neighbours.tolist()  

  return phi_neighbours


def connection_phi_same_layer(src, trgs, nedges=10, nsyns_min=2, nsyns_max=5, divergence=2):
# this function takes the phi as an argument then determines if cells are connected or not
# function returns the list of coordinates of connected cells
# connections are in the same layer
# phi_1 - cell polar coordinate, phi - all phi coordinates available, n - number of connections
# returns the list of connected cells to the phi_1
  # make another array for shifted distances
  phi_1 = src['phi']
  phi = np.array([node['phi'] for node in trgs], dtype=np.float64)
  phi_shifted=np.array(phi)
  #print phi_shifted
  for i in np.arange(0,len(phi)):
    if phi[i]>np.pi:
      phi_shifted[i]=phi[i]-2*np.pi
    else:
      phi_shifted[i]=phi[i]
  # compute the miminum distance
  phi_dist = abs(phi_shifted - phi_1)
  #print phi_dist
  #min_element_idx=np.argmin(phi_dist)
  min_element_idx=phi_dist.argsort()[:2*nedges+1]
  # remove the first element, i.e. itself
  min_element_idx=min_element_idx[1:]
  #print len(min_element_idx)
  #print '\n'
  # randomly set the number of edges according to the divergence pattern
#  print 'nedges: ' + str(len(min_element_idx))
#  print 'divergence: ' + str(divergence)
#  print
  assert(len(min_element_idx) >= divergence)
  min_element_idx=np.random.choice(min_element_idx, divergence, replace=False)
  #print min_element_idx
  # get the indexes of minimal elements
  phi_neighbours = np.zeros(len(phi), dtype=np.uint)
#  phi_neighbours[min_element_idx] = nsyns  # # set up the random number of synapses
  phi_neighbours[min_element_idx] = np.random.randint(nsyns_min,nsyns_max,len(min_element_idx),dtype='int') 
#  print phi_neighbours
  return phi_neighbours



def connection_phi_other_layer(src, trgs, nedges=10, nsyns_min=2, nsyns_max=5, divergence=2):
# this function takes the phi as an argument then determines if cells are connected or not
# function returns the list of coordinates of connected cells
# connections are in the same layer
# phi_1 - cell polar coordinate, phi - all phi coordinates available, n - number of connections
# returns the list of connected cells to the phi_1
  # make another array for shifted distances
  phi_1 = src['phi']
  phi = np.array([node['phi'] for node in trgs], dtype=np.float64)
  phi_shifted=np.array(phi)
  #print phi_shifted
  for i in np.arange(0,len(phi)):
    if phi[i]>np.pi:
      phi_shifted[i]=phi[i]-2*np.pi
    else:
      phi_shifted[i]=phi[i]
  # compute the miminum distance
  phi_dist = abs(phi_shifted - phi_1)
  #print phi_dist
  #min_element_idx=np.argmin(phi_dist)
  min_element_idx=phi_dist.argsort()[:2*nedges]
  # in the other layer it could connect to itself  
  min_element_idx=min_element_idx[:]
  #print len(min_element_idx)
  #print '\n'
  # randomly set the number of edges according to the divergence pattern
  min_element_idx=np.random.choice(min_element_idx, divergence, replace=False)
  #print len(min_element_idx)
  # get the indexes of minimal elements
  phi_neighbours = np.zeros(len(phi), dtype=np.uint)
#  phi_neighbours[min_element_idx] = nsyns  # always n synapses
  phi_neighbours[min_element_idx] = np.random.randint(nsyns_min,nsyns_max,len(min_element_idx),dtype='int')  # always n synapses
#  print phi_neighbours
  return phi_neighbours



def connection_phi_focal_input(src, trgs, nedges=10, nsyns_min=2, nsyns_max=5, phi_1=0):
# this function takes the phi as an argument then determines if cells are connected or not
# function returns the list of coordinates of connected cells
# connections are in the same layer
# phi_1 - cell polar coordinate, phi - all phi coordinates available, n - number of connections
# returns the list of connected cells to the phi_1
  # make another array for shifted distances
  #phi_1 = src['phi']
 
  phi = np.array([node['phi'] for node in trgs], dtype=np.float64)
  phi_shifted=np.array(phi)
  #print phi_shifted
  for i in np.arange(0,len(phi)):
    if phi[i]>np.pi:
      phi_shifted[i]=phi[i]-2*np.pi
    else:
      phi_shifted[i]=phi[i]
  # compute the miminum distance
  phi_dist = abs(phi_shifted - phi_1)
  #print phi_dist
  # get only one connection per external cell
  min_element_idx=phi_dist.argsort()[:nedges]
  min_element_idx=min_element_idx[:]

  # connect only to one element
  #min_element_idx=min_element_idx[1]
  # get the indexes of minimal elements
  phi_neighbours = np.zeros(len(phi), dtype=np.uint)
  # phi_neighbours[min_element_idx] = nsyns  # always n synapses
  phi_neighbours[min_element_idx] = np.random.randint(nsyns_min,nsyns_max,len(min_element_idx),dtype='int')  # always n synapses
  # return all phi neighbours

  # show connected cells, test code
  #print phi_neighbours

  return phi_neighbours


# Step 1: Create GC network coupled synaptically on the ring. Large ring for GC, small ring for BC

# set up the number of cells in different populations

# GC populations
N_GC_1=42
N_GC_2=42
N_GC_3=42
N_GC_4=42
N_GC_5=42
N_GC_6=42
N_GC_7=42
N_GC_8=42
N_GC_9=42
N_GC_10=42
N_GC_11=42
N_GC_12=38
N_GC_all=N_GC_1+N_GC_2+N_GC_3+N_GC_4+N_GC_5+N_GC_6+N_GC_7+N_GC_8+N_GC_9+N_GC_10+N_GC_11+N_GC_12

# BS populations
N_BC_1=3
N_BC_2=3
N_BC_all=N_BC_1+N_BC_2

# generate the polar coordinates for GC
phi_GC = np.linspace(2*np.pi/N_GC_all, 2*np.pi, N_GC_all)
phi_GC = np.random.permutation(phi_GC) # permute the numbers to assign cells randomly
rho_GC=np.ones(len(phi_GC))*800

# generate the polar coordinates for BC
phi_BC = np.linspace(2*np.pi/N_BC_all, 2*np.pi, N_BC_all)
phi_BC = np.random.permutation(phi_BC) # permute the numbers to assign cells randomly
rho_BC=np.ones(len(phi_BC))*750

# set up the index for phi
phi_index=0


net = NetworkBuilder("DG")
net.add_nodes(N=N_GC_1,  # specifiy the number of cells belong to said group.
              pop_name='GC', location='Granule_cell_layer', ei='e',  # pop_name, location, and ei are optional parameters that help's identifies properties of the cells. The modeler can choose whatever key-value pairs as they deem appropiate.
              positions=generate_arc_positions_3D(phi_GC[phi_index:N_GC_1],rho_GC[phi_index:N_GC_1]),  # The following properties we are passing in lists     # of size N. Doing so will uniquely assign different
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_GC[phi_index:N_GC_1],3.130742385893708),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_GC[phi_index:N_GC_1],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_GC[phi_index:N_GC_1],-0.524013338073954),
              phi=phi_GC[phi_index:N_GC_1],
              rho=rho_GC[phi_index:N_GC_1],
              model_type='biophysical',  # The type of cell we are using
              model_template='ctdb:Biophys1.hoc',  # Tells the simulator that when building cells models use a hoc_template specially created for parsing Allen Cell-types file models. Value would be different if we were using NeuronML or different model files
              model_processing='aibs_allactive',  # further instructions for how to processes a cell model. In this case aibs_perisomatic is a built-in directive to cut the axon in a specific way
              dynamics_params='1706014110201.json',  # Name of file (downloaded from Allen Cell-Types) used to set model parameters and channels
              morphology_file='1706014110201.swc')  # Name of morphology file downloaded

phi_index=phi_index+N_GC_1




net.add_nodes(N=N_GC_2, pop_name='GC', location='Granule_cell_layer', ei='e',
              positions=generate_arc_positions_3D(phi_GC[phi_index:phi_index+N_GC_2],rho_GC[phi_index:phi_index+N_GC_2]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_2],-3.0679289479122605),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_2],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_2],-0.7220925381616236),              
              phi=phi_GC[phi_index:phi_index+N_GC_2],
              rho=rho_GC[phi_index:phi_index+N_GC_2],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='1706014110202.json',
              morphology_file='1706014110202.swc')
phi_index=phi_index+N_GC_2



net.add_nodes(N=N_GC_3, pop_name='GC', location='Granule_cell_layer', ei='e',
              positions=generate_arc_positions_3D(phi_GC[phi_index:phi_index+N_GC_3],rho_GC[phi_index:phi_index+N_GC_3]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_3],-3.002976341092172),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_3],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_3],1.2146231788969273),                                       
              phi=phi_GC[phi_index:phi_index+N_GC_3],
              rho=rho_GC[phi_index:phi_index+N_GC_3],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='1706014110502.json',
              morphology_file='1706014110502.swc')
phi_index=phi_index+N_GC_3


net.add_nodes(N=N_GC_4, pop_name='GC', location='Granule_cell_layer', ei='e',
              positions=generate_arc_positions_3D(phi_GC[phi_index:phi_index+N_GC_4],rho_GC[phi_index:phi_index+N_GC_4]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_4],-3.100357711390048),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_4],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_4],0.2392386920845879),                                      
              phi=phi_GC[phi_index:phi_index+N_GC_4],
              rho=rho_GC[phi_index:phi_index+N_GC_4],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='1706014110703.json',
              morphology_file='1706014110703.swc')
phi_index=phi_index+N_GC_4

net.add_nodes(N=N_GC_5, pop_name='GC', location='Granule_cell_layer', ei='e',
              positions=generate_arc_positions_3D(phi_GC[phi_index:phi_index+N_GC_5],rho_GC[phi_index:phi_index+N_GC_5]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_5],0.34531062164929244),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_5],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_5],-1.223480622451367),                                
              phi=phi_GC[phi_index:phi_index+N_GC_5],
              rho=rho_GC[phi_index:phi_index+N_GC_5],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='1706014110704.json',
              morphology_file='1706014110704.swc')
phi_index=phi_index+N_GC_5

net.add_nodes(N=N_GC_6, pop_name='GC', location='Granule_cell_layer', ei='e',
              positions=generate_arc_positions_3D(phi_GC[phi_index:phi_index+N_GC_6],rho_GC[phi_index:phi_index+N_GC_6]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_6],0.16694361153334078),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_6],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_6],0.8339207965166622),                                  
              phi=phi_GC[phi_index:phi_index+N_GC_6],
              rho=rho_GC[phi_index:phi_index+N_GC_6],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='1706015210802.json',
              morphology_file='1706015210802.swc')
phi_index=phi_index+N_GC_6

net.add_nodes(N=N_GC_7, pop_name='GC', location='Granule_cell_layer', ei='e',
              positions=generate_arc_positions_3D(phi_GC[phi_index:phi_index+N_GC_7],rho_GC[phi_index:phi_index+N_GC_7]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_7],0.10355188205115783),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_7],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_7],0.42402951286132445),                               
              phi=phi_GC[phi_index:phi_index+N_GC_7],
              rho=rho_GC[phi_index:phi_index+N_GC_7],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='1706015210901.json',
              morphology_file='1706015210901.swc')
phi_index=phi_index+N_GC_7

net.add_nodes(N=N_GC_8, pop_name='GC', location='Granule_cell_layer', ei='e',
              positions=generate_arc_positions_3D(phi_GC[phi_index:phi_index+N_GC_8],rho_GC[phi_index:phi_index+N_GC_8]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_8],0.07658838304843384),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_8],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_8],0.2913391930286345),                               
              phi=phi_GC[phi_index:phi_index+N_GC_8],
              rho=rho_GC[phi_index:phi_index+N_GC_8],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='1706015210902.json',
              morphology_file='1706015210902.swc')
phi_index=phi_index+N_GC_8

net.add_nodes(N=N_GC_9, pop_name='GC', location='Granule_cell_layer', ei='e',
              positions=generate_arc_positions_3D(phi_GC[phi_index:phi_index+N_GC_9],rho_GC[phi_index:phi_index+N_GC_9]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_9],-2.9809784679136024),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_9],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_9],0.263187572746322),                                          
              phi=phi_GC[phi_index:phi_index+N_GC_9],
              rho=rho_GC[phi_index:phi_index+N_GC_9],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='1706015210904.json',
              morphology_file='1706015210904.swc')
phi_index=phi_index+N_GC_9

net.add_nodes(N=N_GC_10, pop_name='GC', location='Granule_cell_layer', ei='e',
              positions=generate_arc_positions_3D(phi_GC[phi_index:phi_index+N_GC_10],rho_GC[phi_index:phi_index+N_GC_10]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_10],-3.0783195291389207),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_10],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_10],0.37643107877632764),                            
              phi=phi_GC[phi_index:phi_index+N_GC_10],
              rho=rho_GC[phi_index:phi_index+N_GC_10],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='1706015210906.json',
              morphology_file='1706015210906.swc')
phi_index=phi_index+N_GC_10

net.add_nodes(N=N_GC_11, pop_name='GC', location='Granule_cell_layer', ei='e',
              positions=generate_arc_positions_3D(phi_GC[phi_index:phi_index+N_GC_11],rho_GC[phi_index:phi_index+N_GC_11]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_11],-0.2844161166418587),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_11],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_11],0.6271830350104461),                                    
              phi=phi_GC[phi_index:phi_index+N_GC_11],
              rho=rho_GC[phi_index:phi_index+N_GC_11],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='1706015210907.json',
              morphology_file='1706015210907.swc')              
phi_index=phi_index+N_GC_11

net.add_nodes(N=N_GC_12, pop_name='GC', location='Granule_cell_layer', ei='e',
              positions=generate_arc_positions_3D(phi_GC[phi_index:phi_index+N_GC_12],rho_GC[phi_index:phi_index+N_GC_12]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_12],-0.013492500548193464),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_12],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_GC[phi_index:phi_index+N_GC_12],0.07790578339305411),                                     
              phi=phi_GC[phi_index:phi_index+N_GC_12],
              rho=rho_GC[phi_index:phi_index+N_GC_12],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='1706015211001.json',
              morphology_file='1706015211001.swc')

phi_index=0

# Note that in the previous cells we set the tuning_angle, but for BC_1 and BC_2 such parameter is absent (as it is not
# applicable for inhibitory cells). The BMTK builder allows heterogeneous cell properties as dictated by the model
net.add_nodes(N=N_BC_1, pop_name='BC', location='Basket_cell_layer', ei='i',
              positions=generate_arc_positions_3D(phi_BC[phi_index:phi_index+N_BC_1],rho_BC[phi_index:phi_index+N_BC_1]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_BC[phi_index:phi_index+N_BC_1],-0.41682952238355436),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_BC[phi_index:phi_index+N_BC_1],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_BC[phi_index:phi_index+N_BC_1],-1.2877151425739988),                            
              phi=phi_BC[phi_index:phi_index+N_BC_1],
              rho=rho_BC[phi_index:phi_index+N_BC_1],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='529807751.json',
              morphology_file='529807751.swc')
phi_index=phi_index+N_BC_1


net.add_nodes(N=N_BC_2, pop_name='BC', location='Basket_cell_layer', ei='i',
              positions=generate_arc_positions_3D(phi_BC[phi_index:phi_index+N_BC_2],rho_BC[phi_index:phi_index+N_BC_2]),
              rotation_angle_xaxis=generate_rotation_angle_1D(rho_BC[phi_index:phi_index+N_BC_2],-3.056112588143585),
              rotation_angle_yaxis=generate_rotation_angle_1D(rho_BC[phi_index:phi_index+N_BC_2],0),
              rotation_angle_zaxis=generate_rotation_angle_1D(rho_BC[phi_index:phi_index+N_BC_2],-0.339536404497591),                                          
              phi=phi_BC[phi_index:phi_index+N_BC_2],
              rho=rho_BC[phi_index:phi_index+N_BC_2],
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_allactive',
              dynamics_params='541536216.json',
              morphology_file='541536216.swc')



# INTERNAL SYNAPTIC CONNECTIONS IN THE NETWORK

# Step 2: We want to connect our network. Just like how we have node-types concept we group our connections into
# "edge-types" that share rules and properties


# EXCITATORY POPULATION, CONNECTIONS


# gEE reccurent excitation
net.add_edges(source={'ei': 'e', 'model_type': 'biophysical'}, target={'ei': 'e', 'model_type': 'biophysical'},
              connection_rule=connection_phi_same_layer,
              connection_params={'nedges': 100, 'nsyns_min' : 2, 'nsyns_max': 5, 'divergence': 50},
              iterator='one_to_all',
              syn_weight=0.00013,
              weight_function='wmax',
              distance_range=[0.0, 150],
              target_sections=['apical', 'basal'],
              delay=0.8,              
              model_template='exp2syn',
              dynamics_params='AMPA_ExcToExc_GC_GC.json')

# additional reccurent connections
net.add_edges(source={'ei': 'e', 'model_type': 'biophysical'}, target={'ei': 'e', 'model_type': 'biophysical'},
              connection_rule=connection_phi_same_layer,
              connection_params={'nedges': 100, 'nsyns_min' : 2, 'nsyns_max': 5, 'divergence': 50},
              iterator='one_to_all',
              syn_weight=0.0000,
              weight_function='wmax',
              distance_range=[0.0, 150],
              target_sections=['apical', 'basal'],
              delay=0.8,              
              model_template='exp2syn',
              dynamics_params='AMPA_ExcToExc_GC_GC.json')


# gEI Feedforward inhibition
net.add_edges(source={'ei': 'e', 'model_type': 'biophysical'}, target={'ei': 'i', 'model_type': 'biophysical'},
              connection_rule=connection_phi_other_layer,
              connection_params={'nedges': 3, 'nsyns_min' : 2, 'nsyns_max': 5, 'divergence': 1},
              iterator='one_to_all',
              syn_weight=0.0047,
              weight_function='wmax',
              distance_range=[0.0, 150],
              target_sections=['apical', 'basal'],
              delay=0.8,              
              model_template='exp2syn',
              dynamics_params='AMPA_ExcToInh_GC_BC.json')


# INHIBITORY POPULATION, CONNECTIONS

# g_II RECCURENT Inhibition
net.add_edges(source={'ei': 'i', 'model_type': 'biophysical'}, target={'ei': 'i', 'model_type': 'biophysical'},
              connection_rule=connection_phi_same_layer,
              connection_params={'nedges': 3, 'nsyns_min' : 2, 'nsyns_max': 5, 'divergence': 2},
              iterator='one_to_all',
              syn_weight=0.0076,  # synaptic weight
              target_sections=['basal'],  # Gives the simulator the target sections and
              distance_range=[0.0, 150],           # distances (from soma) when creating connections
              delay=0.8,              
              model_template='exp2syn',
              dynamics_params='GABA_InhToInh_BC_BC.json')

# g_IE Feedforward inhibition
net.add_edges(source={'ei': 'i', 'model_type': 'biophysical'}, target={'ei': 'e', 'model_type': 'biophysical'},
              connection_rule=connection_phi_other_layer,
              connection_params={'nedges': 140, 'nsyns_min' : 2, 'nsyns_max': 5, 'divergence': 100},
              iterator='one_to_all',
              syn_weight=0.0016,
              weight_function='wmax',
              distance_range=[0.0, 50.0],
              target_sections=['somatic'],
              delay=0.85,
              model_template='exp2syn',
              dynamics_params='GABA_InhToExc_BC_GC.json')


net.build()
net.save(output_dir='network')




# BASELINE EXTERNAL DRIVE TO THE GC population
GC_external_input = NetworkBuilder("GC_external_input")
GC_external_input.add_nodes(N=9000, pop_name='perforant_path', ei='e',
              positions=generate_positions(9000),
              model_type='virtual')
# drive to the GC population
GC_external_input.add_edges(source=GC_external_input.nodes(), target=net.nodes(ei='e'),
             connection_rule=select_random_srcs,
             connection_params={'nconns': 5, 'nsyns_min': 5, 'nsyns_max': 15},
             iterator='all_to_one',
             syn_weight=0.0003,
             weight_function='wmax',
             distance_range=[150.0, 1e20],
             target_sections=['basal','apical'],
             delay=3.0,
             model_template='exp2syn',
			       dynamics_params='AMPA_ExcToExc_perforant_path.json')
GC_external_input.build()
GC_external_input.save(output_dir='network')



# ADDITIONAL DRIVE TO THE GC population
GC_additional_input = NetworkBuilder("GC_additional_input")
GC_additional_input.add_nodes(N=9000, pop_name='perforant_path_additional', ei='e',
              positions=generate_positions(9000),
              model_type='virtual')
# additional drive to the GC population (double the synapse number)
GC_additional_input.add_edges(source=GC_additional_input.nodes(), target=net.nodes(ei='e'),
             connection_rule=select_random_srcs,
             connection_params={'nconns': 5, 'nsyns_min': 5, 'nsyns_max': 15},
             iterator='all_to_one',
             syn_weight=0.00,
             weight_function='wmax',
             distance_range=[150.0, 1e20],
             target_sections=['basal','apical'],
             delay=3.0,
             model_template='exp2syn',
             dynamics_params='AMPA_ExcToExc_perforant_path.json')
GC_additional_input.build()
GC_additional_input.save(output_dir='network')


# FOCAL EXTERNAL DRIVE TO THE GC population
GC_focal_input = NetworkBuilder("GC_focal_input")
GC_focal_input.add_nodes(N=200, pop_name='perforant_path', ei='e',
              positions=generate_positions(200),
              model_type='virtual')
# focal drive to the GC population
GC_focal_input.add_edges(source=GC_focal_input.nodes(), target=net.nodes(ei='e', model_type='biophysical'),
             connection_rule=connection_phi_focal_input,
             connection_params={'nedges': 100, 'nsyns_min': 5, 'nsyns_max': 15, 'phi_1': np.pi},
             iterator='one_to_all',
             syn_weight=0.001,
             weight_function='wmax',
             distance_range=[150.0, 1e20],
             target_sections=['basal','apical'],
             delay=3.0,
             model_template='exp2syn',
             dynamics_params='AMPA_ExcToExc_perforant_path.json')
GC_focal_input.build()
GC_focal_input.save(output_dir='network')



# BASELINE EXTERNAL DRIVE TO THE BC population
BC_external_input = NetworkBuilder("BC_external_input")
BC_external_input.add_nodes(N=9000, pop_name='perforant_path', ei='e',
              positions=generate_positions(9000),
              model_type='virtual')
# drive to the BC population
BC_external_input.add_edges(source=BC_external_input.nodes(), target=net.nodes(ei='i'),
             connection_rule=select_random_srcs,
             connection_params={'nconns': 5, 'nsyns_min': 5, 'nsyns_max': 15},
             iterator='all_to_one',
             syn_weight=0.0005,
             weight_function='wmax',
             distance_range=[150.0, 1e20],
             target_sections=['basal'],
             delay=3.0,
             model_template='exp2syn',
             dynamics_params='AMPA_ExcToInh_perforant_path.json')
BC_external_input.build()
BC_external_input.save(output_dir='network')


# FOCAL EXTERNAL DRIVE TO THE BC population
BC_focal_input = NetworkBuilder("BC_focal_input")
BC_focal_input.add_nodes(N=50, pop_name='perforant_path', ei='e',
              positions=generate_positions(50),
              model_type='virtual')
# focal drive to the BC population
# focal drive to the BC population
BC_focal_input.add_edges(source=BC_focal_input.nodes(), target=net.nodes(ei='i', model_type='biophysical'),
             connection_rule=connection_phi_focal_input,
             connection_params={'nedges': 2, 'nsyns_min': 5, 'nsyns_max': 15, 'phi_1': np.pi},
             iterator='one_to_all',
             syn_weight=0.001,
             weight_function='wmax',
             distance_range=[150.0, 1e20],
             target_sections=['basal','apical'],
             delay=3.0,
             model_template='exp2syn',
             dynamics_params='AMPA_ExcToExc_perforant_path.json')
BC_focal_input.build()
BC_focal_input.save(output_dir='network')


