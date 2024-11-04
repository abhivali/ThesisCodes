"""
Created on Wed May 11 10:19:45 2022

@author: valisama

Applying the resistance network method 
no contacts only bonds
"""

from time import perf_counter
from pathlib import Path
import pickle

from edempy import Deck
import numpy as np

from scipy.spatial import Delaunay
from scipy.sparse import coo_matrix, csc_matrix, csr_matrix, hstack, vstack, eye
from scipy.sparse.linalg import spsolve, lsqr
from scipy.sparse.csgraph import reverse_cuthill_mckee

import matplotlib.pyplot as plt
from sys import getsizeof

"""
Saving the conductivity data in the form of .vtk file
/!\ remove the pyfedic module and its dependent code snips if other saving methods are used
"""
from pyfedic.mesh import Mesh
from pyfedic.cells import T4
from pyfedic.io import write_mesh

# ---------------------------------------------------------------------------------
#%% Function definitions
# ---------------------------------------------------------------------------------
def neighbour_list(current_node,contact_loc):
    """
    function to extract the neighbours based on the contact between particles
    returns -> list of neighbours
    """
    neighbours = []
    for i in np.arange(len(contact_loc)):
        if contact_loc[i,0] == current_node:
            neighbours = np.append(neighbours,contact_loc[i,1])
        if contact_loc[i,1] == current_node:
            neighbours = np.append(neighbours,contact_loc[i,0])
    return neighbours

def DeckExtract(D,T):
    """
    Function to extract data from the edempy deck object sort it and return as a numpy array
    ----------------------------------
    Parameters
    ----------------------------------
    D : edempy deck object 
    T : float
    ----------------------------------
    Returns
    ----------------------------------
    dataSorted : numpy array [Id Posx Posy Posz]
    p_radSorted : numpy array of particle radii
    """
    pids = []
    ppos = []
    p_radUS = []

    ####
    # data for bond formation timestep:
    ####
    for i in range(D.timestep[T].numTypes):
        #condition to check if particles of this type are created
        if D.timestep[T].particle[i].getNumParticles() > 0: 
            p = D.timestep[T].particle[i].getIds()
            pos = D.timestep[T].particle[i].getPositions()
            rad = D.timestep[T].particle[i].getSphereRadii()
            
            pids.extend(p) # list of particle ids
            ppos.extend(pos) # list of particle positions
            p_radUS.extend(rad) # list of particle radii
            
    ppos = np.array(ppos)
    dataUS = np.column_stack((pids,ppos)) #combined data in the array format and is unsorted

    dataSorted = dataUS[dataUS[:,0].argsort()] # Sorted data
    p_radSorted = np.array(p_radUS)[dataUS[:,0].argsort()] # sorted data
    
    return dataSorted, p_radSorted

def ExtractConnectivity(D, T):
    """
    ----------------------------------
    Parameters
    ----------------------------------
    D : edempy deck object 
    T : float. Timestep
    T_steps : list of all the required timesteps in the
    NumParticles : int. number of particles
    ----------------------------------
    Returns
    ----------------------------------
    Returns Interactions as [particle A, particle B]
    IntactBonds : numpy array (num_bonds x 2)
    """
    BondStatus = D.timestep[T].contact.surfSurf.customProperties['BondStatus'].getData() #get bond status
    cont_ids = D.timestep[T].contact.surfSurf.getIds() #pairs of contact ids from deck
    cont_loc = cont_ids - 1 # convert ids to row numbers
    
    # filtering the intact bonds
    #----------------------------------------
    IntactBonds = cont_loc[BondStatus == 1]
    return IntactBonds

def PboundsBondLengthCorrection(Pdata, p_radius, Intactbondlist, X_bounds = [-25.5,25.5], Y_bounds = [-25.5,25.5], CtoC = False):
    """
    This function estimates estimates the bond length of bonded particles of an assembly with 
    periodic boundary condition along X and Y.
    ----------------------------------
    Parameters
    ----------------------------------
    Pdata : numpy array of particle positions [X Y Z], dimention (num_particles x 3)
    p_radius : numpy array of particle radii, dimention (num_particles x 1)
    Intactbondlist : numpy array of bond pairs [A, B], dimention (num_bonds x 2)
    
    X_bounds : list of float, [X_min ,X_max] optional, limits of boundary along x. The default is [-25.5,25.5].
    Y_bounds : list of float, [Y_min ,Y_max] optional, limits of boundary along y. The default is [-25.5,25.5].
    
    CtoC : Boolean, True or False optional, if True estimates the bond as centre to centre distance. The default is False.
    ----------------------------------
    Returns
    ----------------------------------
    bond_length : numpy array of bond length
        dimention (num_bonds x 1)
    """
    Xbl, = np.diff(X_bounds)
    Ybl, = np.diff(Y_bounds)
    
    bond_diff = Pdata[Intactbondlist[:, 0], :] - Pdata[Intactbondlist[:, 1], :]
    
    # identify the indices to correct and modify the value
    bond_diff[ bond_diff[:,0] > Xbl*0.5  , 0] -= Xbl
    bond_diff[ bond_diff[:,0] < -Xbl*0.5 , 0] += Xbl
    bond_diff[ bond_diff[:,1] > Ybl*0.5  , 1] -= Ybl
    bond_diff[ bond_diff[:,1] < -Ybl*0.5 , 1] += Ybl
    
    if CtoC == True:
        bond_length = np.sqrt(np.sum(bond_diff ** 2, axis=1))
    else:
        bond_length = np.abs(np.sqrt(np.sum(bond_diff ** 2, axis=1)) - (p_radius[Intactbondlist[:, 0]] + p_radius[Intactbondlist[:, 1]]))
    
    return bond_length

def Geometry_NodeInteraction(D,T, NumParticles):
    """
    function to identify the particles in contact with the geometry
    # step 1: identifying particles in contact with the start and end nodes (virtual seperator and electrode nodes)
    # step 2: establishing connection between particles and the virtual nodes based on their interaction with geometries
    ----------------------------------
    Parameters
    ----------------------------------
    D : edempy deck object
    T : float, required timestep.
    NumParticles: int, number of particles
    ----------------------------------
    Returns
    ----------------------------------
    VirtualElectrodeNodes : numpy array, particles in contact with the electrode plate.
    VirtualNode_overlaps : numpy array, stack of overlaps of particle and geometry.
    """
    cont_geom_ids = D.timestep[T].contact.surfGeom.getIds() # get surface to geometry contacts
    GeomOverlap = D.timestep[T].contact.surfGeom.getNormalOverlap() # particle geometry overlaps
    ##identifying the start nodes and end nodes through plate contacts
    TopGeomID = deck.timestep[T].geometry['TopPlate'].getTriangleIds()
    BottomGeomID = deck.timestep[T].geometry['BottomPlate'].getTriangleIds()
    
    bottom_loc = np.isin(cont_geom_ids[:, 1], BottomGeomID)
    top_loc = np.isin(cont_geom_ids[:, 1], TopGeomID)

    start_nod = cont_geom_ids[bottom_loc, 0] - 1  # ids to locations
    end_nod = cont_geom_ids[top_loc, 0] - 1  # ids to locations

    # Extract start and end node overlaps
    start_nod_overlap = GeomOverlap[bottom_loc]
    end_nod_overlap = GeomOverlap[top_loc]
    
    VirtualNode_overlaps = np.vstack((start_nod_overlap.reshape(-1,1),end_nod_overlap.reshape(-1,1)))
    
    # virtual nodes later used in resistance network method
    #------------------------------------------------------
    # creating a virtual node(electrode) that is connected to the start nodes
    # assign the row number of the virtual electrode to "maximum row number + 1"
    nodeElec = NumParticles
    nodeElec_cont = np.column_stack((np.full(len(start_nod), nodeElec), start_nod))
    
    # creating a virtual node(seperator) that is connected to the end nodes
    # assign the row number of the virtual seperator to "maximum row number + 2"
    nodeSep = NumParticles + 1
    nodeSep_cont = np.column_stack((np.full(len(end_nod), nodeSep), end_nod))

    #combining all the contacts virtual
    VirtualElectrodeNodes = (np.vstack((nodeElec_cont,nodeSep_cont))).astype(int)
    
    return VirtualElectrodeNodes, VirtualNode_overlaps

def GMat_generation(p_radius, 
                    Gfac_bond, Gfac_P2G,
                    p_bonded,
                    p_floaters, 
                    Virtual_nodes, 
                    P_Type, g_bond = 375, beta = 0.1):
    
    """
    Generates the conductivity matrix based on the current state of particle interactions.
    ----------------------------------
    Parameters
    ----------------------------------
    p_radius : [numpy array] (num_particles x 1) Array of particle radii.
    
    Gfac_bond : [numpy array] (num_bonds x 1) Conductivity factors associated with bonds.
    Gfac_P2G : [numpy array] (num_contacts x 1) Conductivity factors between particle contacts and the geometry.
    
    p_bonded : [numpy array] (num_bonds x 2) Array of bonded particle pairs (indices refer to particle IDs - 1).
    p_floaters : [numpy array] (num_floaters x 1) Array of particles not currently bonded or in contact with other particles.
    Virtual_nodes : [numpy array] (num_virtual_nodes x 2) Array of virtual nodes and their associated particles (indices refer to node and particle IDs - 1).
    
    P_Type : [numpy boolean array] (num_particles x 1) Boolean array indicating particle type: True for Graphite, False for SiOx.

    g_bond : [float, optional] The electrical conductivity of the conductive binder domain (CBD) bonds. 
        Default is 375 S/m.
    beta : [float, optional] Parameter used to estimate the conductivity between particle-geometry interactions. 
        Default is 0.1.
    ----------------------------------
    Returns
    ----------------------------------
    G_coo_All : [scipy.sparse.coo_matrix] Sparse matrix in COO format representing the full conductivity matrix including floaters.
    G_mat : [scipy.sparse.coo_matrix] Sparse matrix in COO format representing the conductivity matrix excluding floaters.
    """
    MatrixSize = len(p_radius) + 2
    G_virtual = g_bond * beta # beta

    # =============================================================================
    # G factor for contacts and bonds
    # =============================================================================
    G_data_bnd = (-g_bond*Gfac_bond).flatten()
    G_data_virtual = (-G_virtual * Gfac_P2G).flatten() 
    
    # =============================================================================
    # matrix generation
    # =============================================================================
    G_coo_bnd = coo_matrix((G_data_bnd , ((p_bonded[:,0]).flatten(),(p_bonded[:,1]).flatten())), shape=(MatrixSize,MatrixSize))
    G_coo_virtual = coo_matrix((G_data_virtual , ((Virtual_nodes[:,0]).flatten(),(Virtual_nodes[:,1]).flatten())), shape=(MatrixSize,MatrixSize))
    
    G_coo2 = G_coo_bnd + G_coo_virtual

    G_non_diag = G_coo2 + G_coo2.transpose()
    
    G_sum = -(G_non_diag.sum(axis = 0))
    G_sum = np.ravel(G_sum)
    G_diag = coo_matrix((G_sum, (np.arange(0,MatrixSize),np.arange(0,MatrixSize))))
    
    G_coo_All = G_diag + G_non_diag
    
    # =============================================================================
    # removing the floaters to avoid issues with matrix invertion
    # =============================================================================
    mask = ~np.isin(np.arange(G_coo_All.shape[0]), p_floaters)
    G_coo_no_floater = G_coo_All.tocsc()[mask, :].tocsr()[:, mask]
    
    return G_coo_All, G_coo_no_floater.tocoo()

def sparse_insert_zero(matrix, indices):
    for index in sorted(indices):
        
        index = int(index) # ensure that its an integer
        num_rows, num_cols = matrix.shape
        new_row_data = np.zeros(num_cols)
        new_col_data = np.zeros(num_rows+1)
        
        new_coo_row = coo_matrix((new_row_data, ([0]*num_cols, range(num_cols))), shape=(1, num_cols))
        new_coo_col = coo_matrix((new_col_data, (range(num_rows+1), [0]*(num_rows+1))), shape=(num_rows+1, 1))

        if index == 0:
            matrix = vstack([new_coo_row, matrix])
            matrix = hstack([new_coo_col, matrix])
        elif index == num_rows:
            matrix = vstack([matrix, new_coo_row])
            matrix = hstack([matrix, new_coo_col])
        else:
            upper_part = matrix.tocsc()[:index, :]
            lower_part = matrix.tocsc()[index:, :]
            matrix = vstack([upper_part, new_coo_row, lower_part])

            left_part = matrix.tocsr()[:, :index]
            right_part = matrix.tocsr()[:, index:]
            matrix = hstack([left_part, new_coo_col, right_part])
            
    return matrix


def solve_conductivity_coo(G_coo_mat, U_0_1 = 1, regularise = False, regularization_value = 1e-8):
    """
    Solves the conductivity matrix using the method described by Oleg, involves refining the conductivity matrix.
    ----------------------------------
    Parameters
    ----------------------------------
    G_coo_mat : [scipy.sparse.coo_matrix] The conductivity matrix in COO (Coordinate) format, representing the conductance between nodes in the network.
    U_0_1 : [float, optional] The potential difference between the electrodes. 
        The default value is 1.
    
    regularization : add a small value to the diagonal to avoid singular matrix
    use_lsqr : if True, uses scipy.sparse.linalg.lsqr instead of spsolve
    ----------------------------------
    Returns
    ----------------------------------
    R_eff : [float] The effective resistance of the network calculated based on the conductivity matrix.
    i_eff : [float] The effective current through the network, derived from the potential difference and the conductivity matrix.
    phi_vecf : [numpy array] The voltage vector representing the potential at each node in the network after solving the system.
    """    
    # Set initial and end node potentials
    phi_start = 0
    phi_end = phi_start - U_0_1

    # Create the current vector with a random initial value
    curr_vec = np.zeros([G_coo_mat.shape[0], 1])
    i_eff = 56  # Arbitrary current for shrinking operation
    curr_vec[-2, 0] = i_eff
    curr_vec[-1, 0] = -i_eff
    
    # Update the current vector based on the potential difference U_0_1
    curr_vec -= (G_coo_mat.tocsc()[:, -1].toarray().reshape(-1, 1)) * U_0_1

    # Create a copy of the matrix and shrink it by merging the last two rows/columns
    Shrink_G_mat = G_coo_mat.copy().tocsc()
    print(f'-- Matrix size: {Shrink_G_mat.shape[0]} x {Shrink_G_mat.shape[1]}')

    Shrink_G_mat[:, -2] += Shrink_G_mat[:, -1]
    Shrink_G_mat[-2, :] += Shrink_G_mat[-1, :]
    Shrink_G_mat = Shrink_G_mat[:-1, :-1]

    # Check for zero elements on the diagonal
    if np.any(Shrink_G_mat.diagonal() == 0):
        print("-- /!\\ Warning: The matrix has zero elements on the diagonal.")

    
    if regularise:
        # Apply regularization to treat the zero diagonals (not needed if the clusters are properly treated)
        Shrink_G_mat += regularization_value * eye(Shrink_G_mat.shape[0])
        print(f"-- Regularization applied with value {regularization_value}.")

    # Store nodal values
    G_Electrodenode = Shrink_G_mat[-1, :].toarray().reshape(-1, 1)

    # Adjust the current vector after shrinking
    curr_vec[-2, 0] += curr_vec[-1, 0]
    Shrink_curr_vec = curr_vec[:-1]

    # Reorder the matrix and vector using Reverse Cuthill-McKee for better numerical stability
    order = reverse_cuthill_mckee(Shrink_G_mat)
    Shrink_G_mat = Shrink_G_mat[order, :][:, order]
    Shrink_curr_vec = Shrink_curr_vec[order]

    # Solve the system 
    phi_vec = spsolve(Shrink_G_mat, Shrink_curr_vec).reshape(-1, 1)

    # Reapply the original node order
    phi_vec = phi_vec[np.argsort(order)]

    # Apply boundary conditions to obtain the final potential vector
    phi_offset = phi_end - phi_vec[-1, 0]
    phi_vecf = phi_vec + phi_offset
    phi_vecf = np.vstack((phi_vecf, 0)).reshape(-1, 1)

    print('-- Matrix solved for voltage vector.')
    
    return G_Electrodenode, phi_vecf

def estimate_EffectiveResistances(G_Electrodenode, phi_vecf, Floaters, VirtualNodes, U_0_1 = 1):
    """
    Estimates the effective resistance and current of a network given a refined conductivity matrix and potential vector.
    ----------------------------------
    Parameters
    ----------------------------------
    G_Electrodenode : [scipy.sparse.csc_matrix]
        The shrunk conductance matrix in Compressed Sparse Column (CSC) format, representing the conductance between nodes 
        in the network, with the last row corresponding to the electrode node.
    phi_vecf : [numpy array] (num_nodes x 1)
        The voltage vector representing the potential at each node in the network.
    Floaters : [numpy array] (num_floaters x 1)
        Array of particle indices that are floaters, meaning they have no interactions with other particles.
    VirtualNodes : [numpy array] (num_virtual_nodes x 2)
        Array of virtual contact pairs, where each row represents a virtual node and the node it connects to.
    U_0_1 : [float, optional] The potential difference between the electrodes. 
        The default value is 1.
    ----------------------------------
    Returns
    ----------------------------------
    R_eff : [float] The effective resistance of the network calculated from the potential difference and the effective current.
    i_eff : [float] The effective current through the network, calculated based on the conductance matrix and the potential differences.
    phi_vecf : [numpy array] (num_nodes x 1) The updated voltage vector, including any inserted values for floater nodes.
    """
    
    if len(Floaters) > 0:
        for v in sorted(Floaters):
            phi_vecf = np.insert(phi_vecf,int(v),np.nan)
            G_Electrodenode = np.insert(G_Electrodenode,int(v),np.nan)
    
    phi_vecf = phi_vecf.reshape(-1,1)
    G_Electrodenode = G_Electrodenode.reshape(-1,1)
    
    # effective current i_eff entering the start node
    print('-- Estimating effective resistance and current')
    i_elecNode = [0]
    
    ElectrodeNode = VirtualNodes[0,0] # get the node number of the electrode
    
    # particles that are connected to the electrode node which is virtual node
    k = VirtualNodes[VirtualNodes[:,0] == ElectrodeNode , 1]
    
    for co in k:
        i_elecNode = np.vstack((i_elecNode,(G_Electrodenode[int(co),0] * (phi_vecf[ElectrodeNode,0] - phi_vecf[int(co),0]))))
    
    i_eff = np.nansum(i_elecNode)
    R_eff = U_0_1/i_eff
    return R_eff, i_eff, phi_vecf

def CorrectLengthPeriodic(Pos1, Pos2, X_bounds , Y_bounds):
    """
    ----------------------------------
    Parameters
    ----------------------------------
    Pos1 : [X,Y,Z] (nx3)
    Pos2 :  [X,Y,Z] (nx3)
    X_bounds : list [Xmin, Xmax]
    Y_bounds : list [Ymin, Ymax]
    ----------------------------------
    Returns     Distance : array (nx1)
    ----------------------------------
    """
    Xbl, = np.diff(X_bounds)
    Ybl, = np.diff(Y_bounds)
    diff = Pos1 - Pos2
    diff[ diff[:,0] > Xbl*0.5  , 0] -= Xbl
    diff[ diff[:,0] < -Xbl*0.5 , 0] += Xbl
    diff[ diff[:,1] > Ybl*0.5  , 1] -= Ybl
    diff[ diff[:,1] < -Ybl*0.5 , 1] += Ybl
    
    Distance = np.sum(diff**2, axis=1)**0.5
    
    return Distance

def Find_Neighbours(pairs, num_spheres):
    """
    ----------------------------------
    Parameters
    ----------------------------------
    pairs = numpy array of pairs (nx2) 
    example : [[0,1],[0,2],[1,2],[2,3]]
    
    num_spheres = integer 
    Number of sphere ids to check
    ----------------------------------
    Returns
    ----------------------------------
    neighbour_list = pairs converted to dictionary of neighbours
    {"0": ["1", "2"],
     "1": ["0","2"],
     "2": ["0","1","3"],
     "3": ["2"]  
     "4": Nan # 4 not present in the pairs
     "5": Nan # 4 not present in the pairs }
    """
    neighbour_list = {i: [] for i in range(num_spheres+1)}

    for pair in pairs:
        neighbour_list[pair[0]].append(pair[1])
        neighbour_list[pair[1]].append(pair[0])

    for key in neighbour_list:
        neighbour_list[key] = sorted(set(neighbour_list[key]))

    return neighbour_list


def find_connected_points(points):
    """
    ----------------------------------
    Parameters
    ----------------------------------
    points : dictionary of neighbours
         {"0": ["1", "2"],
          "1": ["0","2"],
          "2": ["0","1","3"],
          "3": ["2"]  
          "4": Nan # 4 not present in the pairs
          "5": Nan # 4 not present in the pairs }
    ----------------------------------
    Returns
    ----------------------------------
    components : dictionary
        dictionary of connected clusters.

    """
    # if the arrays are too long the search recursion will stop
    if len(points) > getrecursionlimit():
        setrecursionlimit(len(points) + 2) 
    
    def dfs(node, component):
        component.add(node)
        visited.add(node)
        for neighbor in points[node]:
            if neighbor not in visited:
                dfs(neighbor, component)

    visited = set()
    components = []

    for point in points:
        if point not in visited:
            component = set()
            dfs(point, component)
            components.append(component)
    return components


def Exclude_PeriodicBonding(data_L, data_R, X_bounds = [-25.5,25.5], Y_bounds = [-25.5,25.5]):
    """
    particles of data_i are bonded to particles of data_f
    the input arrays are the positions of these particles
    Bonded particles : | A1 B1 | 
                       | A2 B2 |
                       | :  :  |
                       | An Bn |
    ----------------------------------
    Parameters
    ----------------------------------      
    data_L : Position array of particles A  [X Y Z]
    data_R : Position array of particles B  [X Y Z]
    
    X_bounds : list [X_min, X_max], optional
        boundary limits across x. The default is [-25.5,25.5].
    Y_bounds : list [Y_min, Y_max], optional
        boundary limits across y. The default is [-25.5,25.5].
    ----------------------------------
    Returns
    ----------------------------------
    boolean array of Points to exclude
    True : Exclude
    False: do not exclude
    """
    Xbl, = np.diff(X_bounds)
    Ybl, = np.diff(Y_bounds)

    # estimate distance
    delx, dely, delz = (data_R - data_L).T 
    
    # check for the periodic condition
    outxp = delx > Xbl/2
    outxn = delx < -Xbl/2
    
    outyp = dely > Ybl/2
    outyn = dely < -Ybl/2
    
    out = outxp | outxn | outyp | outyn

    return out


#%% 
# =============================================================================
# ############################ MAIN PROGRAM ###################################
# =============================================================================
# =============================================================================
# ##### program settings
# =============================================================================
# /!\ program suitable for bonding v2 condition ensure before starting the simulation

file_name =  r"xxxxxFileLocationxxxxx/SimulationFileName.dem"
                                 
RunForParametricTest = True
VtkSave = True
Binning = False
StepsManualInput = False
ElecNodesFromEDEM = True



# EDEM deck file opening and saving as pickle for optimised future use
# -----------------------------------------------------------------------------
print("________Extracting Deck ________")
file_path = Path(file_name)
pickle_Save = file_path.parent
pickleFileName = file_path.parts[-1][0:-4] +'.pickle'

try:
    saveDeck = False
    with open( pickle_Save / pickleFileName , 'rb') as file:
        deck = pickle.load(file)
    print("________Pickle file available")
    
    #  check if the file is in proper location
    if(deck.__dict__['deckname'] != file_name):
        print('Old pickle file needs update')
        saveDeck = True

except FileNotFoundError:
    print("________Pickle file not found")
    saveDeck = True

if saveDeck == True:

    print("________Extracting deck")
    deck = Deck(file_name) # opening the deck
    print("________Extracting deck : done")
    print("________Creating a Pickle file")

    with open( pickle_Save / pickleFileName, 'wb') as file:
        pickle.dump(deck, file) # saving as pickle file as it is easier to load
# -----------------------------------------------------------------------------


XBounds = [deck.creatorData[0].domain.getXMin(), deck.creatorData[0].domain.getXMax()]
YBounds = [deck.creatorData[0].domain.getYMin(), deck.creatorData[0].domain.getYMax()]


#Binning dimensions (optional)
# -----------------------------------------------------------------------------
if Binning:
    XBounds = [-30,30]
    YBounds = [-30,30]
    ZBounds = [-50,60]
    XboundLen, = np.diff(XBounds)
    YboundLen, = np.diff(YBounds)


#Simulation Timesteps
# -----------------------------------------------------------------------------
t_b = 5 #bond formation timestep
t_start = 40 #start step for swelling
step_size = 25
t_final = len(deck.timestepKeys)-1
t_steps = np.arange(t_start,t_final,step_size).tolist()
if (t_steps[-1] != t_final): t_steps.append(t_final)
t_steps = np.array(t_steps)

if RunForParametricTest == True:
    #Modify t_steps to have max an min timesteps
    SimSwellTimes = np.array([65000,50000,30000,13000,10000,65000,65000,65000,52000,36000,18000,40000,10000]) # Swell times used for each simulation
    sim = int(file_path.parts[-2][3::]) # Extract the simulation number
    
    # use a general step size of 5 and combine the tsteps with required t_steps if not present
    t_gen = np.arange(40,t_final,5).tolist()
    if (t_gen[-1] != t_final): t_gen.append(t_final)
    t_gen = np.array(t_gen)
    
    t_rel =  ((t_gen - 40) / (SimSwellTimes[sim-1]/200)) % 2
    t_max = t_gen[t_rel == 1]
    t_str = t_gen[t_rel == 0]
    
    # t_steps = np.union1d( np.union1d(t_steps,t_max) , t_str)
    t_steps = t_str.copy()
    
if StepsManualInput == True:
    t_steps = np.array([40,300,305])
    
# domain length
# -----------------------------------------------------------------------------
DomainXlen = (deck.timestep[0].domainMax[0] - deck.timestep[0].domainMin[0])

# constants
# -----------------------------------------------------------------------------
# /!\ The values of conductivity do not matter as the code aims only for the normalised impact
G_bond = 375 # electronic conductivity of CBD = 0.375e3[S/m] reference: DOI:10.1149/2.0981813jes
BondDiskScale = 0.2

gamma = 0.2# Contact radius for overlap correction
model_num = 3


# Save voltages in a visualisation file
# -----------------------------------------------------------------------------
if VtkSave:
    # establish the saving folder
    save_folder = Path(file_name).parent
    Sub_folder = 'Conductivity_Files_test_onlyBinder'
    (save_folder / Sub_folder).mkdir(parents=True, exist_ok=True)


print(f"________Processing data for the bond formation time step {t_b} ________")
data0, p_rad0 = DeckExtract(deck, t_b)
# =============================================================================
# ##### file processing and extracting basic data from the deck 
# =============================================================================
dict_data = {}
dict_rad = {}
dict_IntactCont = {}
dict_IntactBnd = {}
dict_BrokenBndCont = {}
dict_G_coo_combined = {}
dict_VoltageVector = {}
dict_R_eff = {}
dict_i_eff = {}
dict_G_eff = {}
dict_k_eff = {}
dict_Virtualconts = {}
dict_GeomPressure ={}

# search for the propper bin data locations
# -----------------------------------------------------------------------------
if Binning:
    xSearch = np.logical_and( data0[:, 1] > XBounds[0] , data0[:, 1] < XBounds[1])
    ySearch = np.logical_and( data0[:, 2] > YBounds[0] , data0[:, 2] < YBounds[1])
    zSearch = np.logical_and( data0[:, 3] > ZBounds[0] , data0[:, 3] < ZBounds[1])
    
    BinSearch = np.logical_and(xSearch, ySearch, zSearch)
    BinIDs = data0[BinSearch, 0]
    BinLocs = BinIDs-1
    
    data0 = data0[BinSearch] # replace with new data
    p_rad0 = p_rad0[BinSearch] # replace




# Identification of particle type
# -----------------------------------------------------------------------------
# Graphite = np.asarray(["A","B","C","D","E","F","G"])
Graphite_Radii = np.asarray([1.99, 2.8, 3.6, 4.8, 6, 6.5, 7.5])
# Silicon = np.asarray(["H","I","J","K","L","M","N"])
Silicon_Radii =  np.asarray([0.61, 0.81, 1.02, 1.33, 1.63, 1.84, 2.04])
PType = np.array([False] * len(p_rad0))
for i in Graphite_Radii:
    PType[p_rad0 == i] = True # identify graphite particles and set type to True


# Save initial mesh for future use to save data (optional)
# -----------------------------------------------------------------------------
if VtkSave:
    print("________Generating initial mesh for particles")
    nodes = data0[:,1:4]
    cells = Delaunay(nodes).simplices
    m = Mesh.new(nodes, {T4: cells})
    

for t in t_steps:
    print(f"________Processing data for the time step: {t} ________")
    # extracting values for the required time step
    data, p_rad = DeckExtract(deck, t) 
    thick_z = abs(np.max(data[:,3]) - np.min(data[:,3]))
    
    
    
    print("________Extracting Particle connectivity data")
    # Extracting particle connectivity data
    # IntactCont, IntactBnd, BrokenBndCont, ContOverlap, ContOverlappos = ExtractConnectivity(deck, t)
    IntactBnd = ExtractConnectivity(deck, t)
    
    
    
    if t == t_steps[0] or t == t_steps[-1]:
        posi1 = data[IntactBnd[:,0],1:4]
        posi2 = data[IntactBnd[:,1],1:4]

        if t == t_steps[0]:
            fig, axs = plt.subplots(1, 3, subplot_kw={'polar': True}, figsize=(20, 4))
            
            plotBondpolor(posi1, posi2, fig, axs, c='b')
        else:
            plotBondpolor(posi1, posi2, fig, axs, c='r', superpose=True)
    
    
    #Select the interactions with good ids that are binned
    # -----------------------------------------------------------------------------
    if Binning:
        # Initiate Arrays
        BondBool0 = np.array([False] * len(IntactBnd))
        BondBool1 = np.array([False] * len(IntactBnd))
        # search for particles using their Location in Data
        # as the ids are converted to locations in the interaction extraction
        for Loc in BinLocs:
            BondBool0[ IntactBnd[:,0] == Loc] = True
            BondBool1[ IntactBnd[:,1] == Loc] = True
        # Accept the interactions if both the particles are inside the Bin and replace the array
        IntactBnd = IntactBnd[np.logical_and(BondBool0, BondBool1)]
    
    
    
    # Identifying particles in contact with the geometries
    # -----------------------------------------------------------------------------
    if(t == t_steps[0]):
        if ElecNodesFromEDEM:
            VnodeStack, VnodesOverlaps = Geometry_NodeInteraction(deck,t, NumParticles = data.shape[0])
            if Binning:
                # initialise boolean array
                VnodeBool1= np.array([False] * len(VnodeStack))
                # particle ids are only in the first column of the stack
                for Loc in BinLocs:
                    VnodeBool1[ VnodeStack[:,1] == Loc] = True
                VnodeStack = VnodeStack[VnodeBool1]
                VnodesOverlaps = VnodesOverlaps[VnodeBool1]
        else:
            # /!\ considering that the positions of particles are within -40 and +40 during this step
            # if not change the values rspectively 
            # needs to be updated to make it compatible with binning
            Z_min = np.min(data[:,3])
            Z_max = np.max(data[:,3])
            
            CollectorNodesB = (data[:, 3] - Z_min) <= (p_rad * 1.25)
            CollectorNodesT = (Z_max - data[:, 3])  <= (p_rad * 1.25)
            
            VnodeB = CollectorNodesB.nonzero()[0]
            VnodeT = CollectorNodesT.nonzero()[0]
            
            OverlapB = ((p_rad * 1.1) - (data[:, 3] - Z_min))[VnodeB]
            OverlapT = ((p_rad * 1.1) - (Z_max - data[:, 3]))[VnodeT]
            
            StartStack = np.column_stack((np.full(len(VnodeB), data.shape[0]) , VnodeB))
            EndStack = np.column_stack((np.full(len(VnodeT), data.shape[0] + 1 ) , VnodeT))
            
            VnodeStack = (np.vstack((StartStack, EndStack))).astype(int)
            VnodesOverlaps = (np.hstack((OverlapB, OverlapT))).astype(int)
            
            
    # identifying the floaters (Including inactive clusters that are bonded):
    # -----------------------------------------------------------------------------
    all_interactions = np.vstack((IntactBnd, VnodeStack))
    
    # Remove inactive clusters: identify the ones that dont have the electrode nodes 
    neighbours_b = Find_Neighbours(IntactBnd, data.shape[0])
    bonded_clusters_b = find_connected_points(neighbours_b)

    V1 = set(VnodeStack[VnodeStack[:,0] == VnodeStack[0,0], 1])
    V2 = set(VnodeStack[VnodeStack[:,0] == VnodeStack[-1,0],1])
    
    active_Clusters = []
    for cluster in bonded_clusters_b:
        temp1 = cluster.intersection(V1)
        temp2 = cluster.intersection(V2)
        if len(temp1)>0 and len(temp2)>0:
            active_Clusters = np.append(active_Clusters, np.array(list(cluster)))
    active_Clusters = np.append(active_Clusters, [data.shape[0], data.shape[0]+1])
    
    if Binning:
        floater = np.setdiff1d(BinLocs, active_Clusters)
    else :  
        floater = np.setdiff1d(np.arange(len(data)), active_Clusters)
    
    if Binning:
        # create new data set with all the data within the bin
        dataBin = data[BinSearch,:]
        IntactBnd_bin = np.copy(IntactBnd)
        # BrokenBndCont_bin = np.copy(BrokenBndCont)
        VnodeStack_bin = np.copy(VnodeStack)
        floater_bin = np.copy(floater)
        
        for i, loc in enumerate(BinLocs):
            IntactBnd_bin[IntactBnd == loc] = i
            # BrokenBndCont_bin[BrokenBndCont == loc] = i
            VnodeStack_bin[VnodeStack == loc] = i
            floater_bin[floater == loc] = i
        VnodeStack_bin[VnodeStack == len(data)] = len(dataBin)
        VnodeStack_bin[VnodeStack == len(data)+1] = len(dataBin)+1
        
        data = np.copy(dataBin)
        p_rad = p_rad[BinSearch]
        IntactBnd = np.copy(IntactBnd_bin)
        VnodeStack = np.copy(VnodeStack_bin)
        floater = np.copy(floater_bin)
    
    
    
    print("________Applying RNM")
    # =============================================================================
    # # Applying resitance network method
    # =============================================================================   
    # # maintaining constant volume of bonds
    # virtual bond area will be calculated during the conductivity matrix generation
    # it represents the change of area of the bond during the change in length to 
    # keep the volume constant
    # implemented to avoid over estimation and under estimation of the conductivity matrix
    if(t == t_steps[0]):
        InitialList = np.copy(IntactBnd)
        BondLength_i = PboundsBondLengthCorrection(data0[:,1:], p_rad0, IntactBnd)
        BondLength_t = np.copy(BondLength_i)
        BondDiskArea_i = np.power(((p_rad0[IntactBnd]).min(axis=1)) * BondDiskScale, 2) * np.pi
        BondDiskArea_t = np.copy(BondDiskArea_i)
    else:
        BondLength_t = np.zeros(len(IntactBnd))
        BondDiskArea_t = np.zeros(len(IntactBnd))
        
        # Create a dictionary for quick lookup of InitialList indices
        bond_map = {tuple(bond): idx for idx, bond in enumerate(InitialList)}
        
        # Iterate through IntactBnd and check both the original and reversed bond pairs 
        for idx, bond in enumerate(IntactBnd):
            bond_tuple = tuple(bond)
            reverse_bond_tuple = tuple(bond[::-1])
            
            if bond_tuple in bond_map:
                match_idx = bond_map[bond_tuple]
            elif reverse_bond_tuple in bond_map:
                match_idx = bond_map[reverse_bond_tuple]
            else:
                continue
            BondLength_t[idx] = BondLength_i[match_idx]
            BondDiskArea_t[idx] = BondDiskArea_i[match_idx]
    
    
    # conductivity matrix generation:
    # -----------------------------------------------------------------------------
    print('________Generating G_mat')
    cpu0 = perf_counter()
    
    #### Estimate the conductivity factor 
    # -----------------------------------------------------------------------------
    # Bond conductivity
    # -----------------------------------------------------------------------------
    # bond_Gfac = np.divide(BondDiskArea_t,BondLength_t)
    bond_Gfac = np.ones(len(BondDiskArea_t))
    
    
    # conductivity between particles and the plate particle to geometry (P2G)
    # -----------------------------------------------------------------------------
    if (t == t_steps[0]) and (model_num == 2 or model_num == 3):
        VnodesOverlaps = VnodesOverlaps.flatten() + gamma*(p_rad[VnodeStack[:,1]])
        VparticlesRadii = p_rad[VnodeStack[:,1]]
        P2G_Gfac = np.pi * (((1+gamma)*VparticlesRadii * 2) - VnodesOverlaps.flatten()) # A/L for the contact of particle with geometry => A/L = pi * (2R - h)
    
    
    # modifying the contact overlaps to consider the contact radius for Gmat estimation
    # -----------------------------------------------------------------------------
    G_coo_combined, G_mat_no_Floaters = GMat_generation(p_rad, 
                                                        bond_Gfac, 
                                                        P2G_Gfac,
                                                        IntactBnd, floater, VnodeStack,
                                                        PType, g_bond = G_bond, beta = 10)
    
    print(f'-- Matrix size {getsizeof(G_coo_combined)} bytes')
    
    cpu1 = perf_counter()
    print('-- time elapsed',"{:.4e}".format(cpu1-cpu0),'s')
    
    
    
    # solving matrix for voltage vectors:
    # -----------------------------------------------------------------------------
    print('Solving G_mat')
    cpu0 = perf_counter()
    
    G_mat_nodes, VoltageVector_shrinked = solve_conductivity_coo(G_mat_no_Floaters, regularise = True)
    R_eff, i_eff, VoltageVector = estimate_EffectiveResistances(G_mat_nodes, VoltageVector_shrinked, floater, VnodeStack)
    
    G_eff = 1/R_eff # effective conductance
    if Binning:
        k_eff = G_eff * (thick_z/(XboundLen*YboundLen)) # effective conductivity of the domain
    else:
        k_eff = G_eff * (thick_z/(DomainXlen*DomainXlen)) # effective conductivity of the domain
    
    cpu1 = perf_counter()
    print('-- time elapsed',"{:.4e}".format(cpu1-cpu0),'s')
    
    
    
    
    
    
    if RunForParametricTest:
        t_rel = ((t - 40) / (SimSwellTimes[sim-1]/200)) % 2
    else:
        t_rel = 0
        
    if VtkSave:
        if t_rel == 0 or t_rel == 1:
            VoltageVector = np.absolute(VoltageVector)
            VoltageVector[VoltageVector < 1e-6] = 0
            write_mesh(save_folder / Sub_folder / f"K{t:04}.vtk", m, point_values = {'Radius': p_rad0.flatten() , 'Type' : PType.flatten(), 'Potential': VoltageVector.flatten()})
    
    
    dict_data[t] = data
    dict_rad[t] = p_rad
    dict_IntactBnd[t] = IntactBnd
    dict_G_coo_combined[t] = G_coo_combined
    dict_VoltageVector[t] = VoltageVector
    dict_R_eff[t] = R_eff
    dict_i_eff[t] = i_eff
    dict_G_eff[t] = G_eff
    dict_k_eff[t] = k_eff
    dict_Virtualconts[t] = VnodeStack
    dict_GeomPressure[t] = np.sum(deck.timestep[t].geometry['TopPlate'].getPressure())/2
         
    if StepsManualInput == True:
         np.save(save_folder / Sub_folder / f"data_{t}", data)
         np.save(save_folder / Sub_folder / f"rad_{t}", p_rad)
         np.save(save_folder / Sub_folder / f"IntactBnd_{t}", IntactBnd)
         np.save(save_folder / Sub_folder / f"R_eff_{t}", R_eff)
         np.save(save_folder / Sub_folder / f"k_eff_{t}", k_eff)
         np.save(save_folder / Sub_folder / f"VoltageVector_{t}", VoltageVector)
    
    
    
if RunForParametricTest == True:
     np.save(save_folder / Sub_folder / "dict_data", dict_data)
     np.save(save_folder / Sub_folder / "dict_rad", dict_rad)
     np.save(save_folder / Sub_folder / "dict_R_eff", dict_R_eff)
     np.save(save_folder / Sub_folder / "dict_k_eff", dict_k_eff)
     np.save(save_folder / Sub_folder / "dict_VoltageVector", dict_VoltageVector)
     np.save(save_folder / Sub_folder / "dict_IntactBnd", dict_IntactBnd)


    
