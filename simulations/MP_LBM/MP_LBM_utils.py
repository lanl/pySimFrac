from subprocess import run
import os
from os import mkdir
from shutil import copy, move
from argparse import Namespace

import pandas as pd
import numpy as np

import pickle
from copy import copy as cp

from skimage import measure
from scipy.ndimage import distance_transform_edt as edist




def read_permeability():
    """
    Reads the permeability value from the data file 'frac.txt'. 

    The function retrieves the permeability value, which is expected to be found after the first 
    equal sign ('=') and before the substring '\nName' on the fourth line from the end of the file.

    It uses the pandas library to read and handle the data file. If the permeability value cannot be 
    converted to a float or if the file 'frac.txt' cannot be found, exceptions will be raised.
    
    Parameters
    -----------
        None
        
    Returns
    --------
        perm : float
            The permeability value parsed from 'frac.txt'.

    Notes
    ------
        None
    
    Raises
    -------
        ValueError
            If the permeability value is not a valid float.
        FileNotFoundError
            If 'frac.txt' cannot be found.
    """
    # Read the data file using pandas
    file = pd.read_csv('frac.txt', header=None)

    # Find the index of the first equal sign in the last fourth line
    ind1 = str(file.iloc[-4]).find('=')+2
    # Find the index of the substring '\nName' in the last fourth line
    ind2 = str(file.iloc[-4]).find('\nName')

    # Extract the permeability value from the string between the two indices and convert it to float
    perm = float(str(file.iloc[-4])[ind1:ind2])
    return perm





def check_lbm_install(sim_mode='single-phase'):
    """
    Checks whether the MPLBM software has been installed. 

    The function verifies the presence of the 'permeability' or 'ShanChen' file in the 
    corresponding directory based on the simulation_mode input. 
    If the file is not present, it indicates that the LBM software is not installed and raises a 
    NotImplementedError. If the file is present, it prints a message to the console to indicate that the 
    LBM installation was found.

    Parameters
    -----------
    simulation_mode: str, optional
        Simulation mode ('single-phase' or 'two-phase'). Defaults to 'single-phase'.

    Returns
    --------
        None
        
    Raises
    -------
        NotImplementedError
            If the required file cannot be found in the corresponding directory, suggesting that 
            the LBM software has not been installed.

    Notes
    ------
        The function makes use of the 'os' library to interact with the operating system.

    Examples
    --------
        >>> check_lbm_install('single-phase')
        LBM installation was found
    """
    if sim_mode == 'single-phase':
        file_path = 'lbm_light/src/1-phase_LBM/permeability'
    elif sim_mode == 'two-phase':
        file_path = 'MPLBM-UT/src/2-phase_LBM/ShanChen'
    else:
        raise ValueError(f'Invalid simulation_mode: {sim_mode}')

    if not os.path.isfile(file_path):
        raise NotImplementedError(f'\nLBM for {sim_mode} mode has not been installed\n')
    else:
        print(f'\nLBM for {sim_mode} mode installation was found\n')

def create_folder(folder_name):
    """
    Creates a new folder with a given name.

    The function attempts to create a new directory with the specified name. If a directory with the 
    same name already exists, the function will not create a new directory and will not raise any errors.

    Parameters
    -----------
        folder_name : str
            The name of the folder to be created. 

    Returns
    -------
        None
        
    Raises
    -------
        FileExistsError
            If a directory with the specified name already exists. This error is caught internally and 
            does not interrupt the program execution.

    Notes
    ------
        The function uses the 'os' library function 'mkdir' to create the new directory.

    """
    try:
        mkdir(folder_name)
    except FileExistsError:
        pass
   
    
        
    
def replace_word(file, old_word, new_word):
    """
    Replaces occurrences of a word in a file with a new word.

    The function reads the contents of the specified file, replaces all occurrences of 'old_word' 
    with 'new_word', and overwrites the file with the new contents. If the file does not exist or 
    cannot be read, a FileNotFoundError will be raised. If the file cannot be written to, an 
    OSError may be raised.

    Parameters
    ----------
        file : str
            The path of the file in which the word will be replaced.
        old_word : str
            The word to be replaced.
        new_word : str
            The word to replace 'old_word' with.

    Returns
    --------
        None
        
    Raises
    ------
        FileNotFoundError
            If the specified file does not exist or cannot be read.
        OSError
            If the file cannot be written to.

    Notes
    -----
        None 

    """
    
    f = open(file,'r')
    filedata = f.read()
    f.close()
    
    newdata = filedata.replace(old_word, new_word)
    
    f = open(file,'w')
    f.write(newdata)
    f.close()

def create_single_phase_input_file(input_file_name, geom_name, domain_size, periodic, io_folders, sim_settings, save_vtks):
    """
    Creates an XML input file for a single phase flow simulation based on given parameters.

    The function uses the provided parameters to write the settings for a single phase flow simulation 
    into an XML file. If a file with the specified name already exists, it will be overwritten.

    Parameters
    ----------
        input_file_name : str
            The name of the input file to create. If it already exists, it will be overwritten.
        geom_name : str
            The name of the .dat file that contains the entire geometry. Do not include the ".dat" extension.
        domain_size : list
            The size (in voxels) of the simulation domain, provided as a list of the form [x, y, z].
        periodic : list
            Determines whether the simulation should be periodic in the x, y, and z directions. This should be a list of three strings, either "True" or "False".
        io_folders : list
            A list containing the paths to the input and output folders, in that order.
        sim_settings : list
            A list containing the number of geometries/simulations, pressure, maximum number of iterations, and convergence, in that order.
        save_vtks : str
            Determines whether to save VTK files for the medium and velocity. This should be a string, either "True" or "False".

    Returns
    --------
        None
        
    Notes
    -----
        The function uses Python's built-in 'open' function to create and write to the XML file.

    Examples
    --------
        >>> create_single_phase_input_file('input.xml', 'geometry', [100, 100, 100], ['True', 'False', 'True'], ['input/', 'output/'], [5, 0.02, 2000, 0.0001], 'True')
        # This will create 'input.xml' with the specified settings.
    
    """
    # Parse geometry inputs
    x_size, y_size, z_size = domain_size
    periodic_x, periodic_y, periodic_z = periodic

    # Parse I/O inputs
    input_folder, output_folder = io_folders

    # Parse simulation inputs
    num_geoms_or_sims, pressure, max_iter, convergence = sim_settings

    # Open the input file, creating it if necessary
    with open(input_file_name, 'w') as file:
        file.write('<?xml version="1.0" ?>\n\n')  # Write XML header

        # Write the geometry settings to the file
        file.write('<geometry>\n')
        file.write(f'\t<file_geom> {geom_name} </file_geom>\n')
        file.write(f'\t<size> <x> {x_size} </x> <y> {y_size} </y> <z> {z_size} </z> </size>\n')
        file.write(f'\t<per> <x> {periodic_x} </x> <y> {periodic_y} </y> <z> {periodic_z} </z> </per>\n')
        file.write('</geometry>\n\n')

        # Write the I/O settings to the file
        file.write('<folder>\n')
        file.write(f'\t<in_f> {input_folder} </in_f>\n')
        file.write(f'\t<out_f> {output_folder} </out_f>\n')
        file.write('</folder>\n\n')

        # Write the simulation settings to the file
        file.write('<simulations>\n')
        file.write(f'\t<num> {num_geoms_or_sims} </num>\n')
        file.write(f'\t<press> {pressure} </press>\n')
        file.write(f'\t<iter> {max_iter} </iter>\n')
        file.write(f'\t<conv> {convergence} </conv>\n')
        file.write(f'\t<vtk_out> {save_vtks} </vtk_out>\n')
        file.write('</simulations>')

def create_two_phase_input_file(input_file_name, geom_name, domain_size, periodic, io_folders, sim_settings):
    # Generate input xml file
    x_size, y_size, z_size = domain_size
    periodic_x, periodic_y, periodic_z = periodic

    # Parse I/O inputs
    input_folder, output_folder = io_folders

    # Parse simulation inputs
    sim_max_iter, sim_convergence, convergence_relperm,max_iterations_relperm, fluid_1_init, fluid_2_init, Gc, omega_f1, omega_f2,G_ads_f1_s1, G_ads_f1_s2, G_ads_f1_s3, G_ads_f1_s4, convergence_iter, convergence_iter,save_sim, save_iter, gif_iter, vtk_iter, rho_f2_vtk, rho_f2_vtk, print_geom, print_stl,save_vtks, pressure_bc, minimum_radius, num_pressure_steps, force_f1, force_f2, rho_f1, rho_f2 = sim_settings

    input_xml_file = '2_phase_sim_input.xml'
    
    restart_sim = 'False'
    
    if pressure_bc == True:
        minimum_radius = minimum_radius
        num_pc_steps = num_pressure_steps
    else:
        minimum_radius = 1
        num_pc_steps = 0
    
    load_fluid_type = 'geom'
    if load_fluid_type == 'geom':
        load_fluid_from_geom = True
    else:
        load_fluid_from_geom = False

    with open(input_file_name, 'w') as file:    
        file.write('<?xml version="1.0" ?>\n\n')  # Write xml header
        
        # Restart sim?
        file.write(f'<load_savedstated> {restart_sim} </load_savedstated>\n\n')
        
        # Write geometry section
        file.write('<geometry>\n')
        # Geometry name
        file.write(f'\t<file_geom> input/{geom_name}.dat </file_geom>\n')
        # Geometry size
        file.write(f'\t<size> <x> {x_size} </x> <y> {y_size} </y> <z> {z_size} </z> </size>\n')
        # Periodicity
        file.write(f'\t<per>\n')
        file.write(f'\t\t<fluid1> <x> {periodic_x} </x> <y> {periodic_y} </y> <z> {periodic_z} </z> </fluid1>\n')
        file.write(f'\t\t<fluid2> <x> {periodic_x} </x> <y> {periodic_y} </y> <z> {periodic_z} </z> </fluid2>\n')
        file.write(f'\t</per>\n')
        file.write('</geometry>\n\n')
        
        # Write initial position of fluids
        file.write(f'<init>\n')
        file.write(f'\t<fluid_from_geom> {load_fluid_from_geom} </fluid_from_geom>\n')
        file.write(f'\t<fluid1>\n')
        file.write(f'\t\t <x1> {fluid_1_init[0]} </x1> <y1> {fluid_1_init[1]} </y1> <z1> {fluid_1_init[2]} </z1>\n')
        file.write(f'\t\t <x2> {fluid_1_init[3]} </x2> <y2> {fluid_1_init[4]} </y2> <z2> {fluid_1_init[5]} </z2>\n')
        file.write(f'\t</fluid1>\n')
        file.write(f'\t<fluid2>\n')
        file.write(f'\t\t <x1> {fluid_2_init[0]} </x1> <y1> {fluid_2_init[1]} </y1> <z1> {fluid_2_init[2]} </z1>\n')
        file.write(f'\t\t <x2> {fluid_2_init[3]} </x2> <y2> {fluid_2_init[4]} </y2> <z2> {fluid_2_init[5]} </z2>\n')
        file.write(f'\t</fluid2>\n')
        file.write('</init>\n\n')
        
        # Write fluid data
        file.write('<fluids>\n')
        
        file.write(f'\t<Gc> {Gc} </Gc>\n')
        file.write(f'\t<omega_f1> {omega_f1} </omega_f1>\n')
        file.write(f'\t<omega_f2> {omega_f2} </omega_f2>\n')
        file.write(f'\t<force_f1> {force_f1} </force_f1>\n')
        file.write(f'\t<force_f2> {force_f2} </force_f2>\n')
        file.write(f'\t<G_ads_f1_s1> {G_ads_f1_s1} </G_ads_f1_s1>\n')
        file.write(f'\t<G_ads_f1_s2> {G_ads_f1_s2} </G_ads_f1_s2>\n')
        file.write(f'\t<G_ads_f1_s3> {G_ads_f1_s3} </G_ads_f1_s3>\n')
        file.write(f'\t<G_ads_f1_s4> {G_ads_f1_s4} </G_ads_f1_s4>\n')
        
        file.write(f'\t<rho_f1> {rho_f1} </rho_f1>\n')
        file.write(f'\t<rho_f2> {rho_f2} </rho_f2>\n')
        
        file.write(f'\t<pressure_bc> {pressure_bc} </pressure_bc>\n')
        file.write(f'\t<rho_f1_i> {rho_f1} </rho_f1_i>\n')
        file.write(f'\t<rho_f2_i> {rho_f2} </rho_f2_i>\n')
        file.write(f'\t<num_pc_steps> {num_pc_steps} </num_pc_steps>\n')
        file.write(f'\t<min_radius> {minimum_radius} </min_radius>\n')
        file.write(f'\t<rho_d> 0.06 </rho_d>\n')
        
        file.write('</fluids>\n\n')
        
        # Write output section
        file.write('<output>\n')
        
        file.write(f'\t<out_folder> {output_folder} </out_folder>\n')
        file.write(f'\t<save_it> {save_iter} </save_it>\n')
        file.write(f'\t<save_sim> {save_sim} </save_sim>\n')
        file.write(f'\t<convergence> {sim_convergence} </convergence>\n')
        file.write(f'\t<it_max> {sim_max_iter} </it_max>\n')
        file.write(f'\t<it_conv> {convergence_iter} </it_conv>\n')
        file.write(f'\t<it_gif> {gif_iter} </it_gif>\n')
        file.write(f'\t<rho_vtk> {rho_f2_vtk} </rho_vtk>\n')
        file.write(f'\t<it_vtk> {vtk_iter} </it_vtk>\n')
        file.write(f'\t<print_geom> {print_geom} </print_geom>\n')
        file.write(f'\t<print_stl> {print_stl} </print_stl>\n')
        
        file.write('</output>')

    #return


    
class SimulationInitializer:
    def __init__(self, init_type):
        self.init_type = init_type

    def initialize_simulation_matrix(self, rock, wetting_saturation_ratio=0.5, args=None):
        if self.init_type == 'random':
            return self.random_simulation_matrix(rock, wetting_saturation_ratio, args)
        elif self.init_type == 'geometric':
            return self.geometric_simulation_matrix(rock, wetting_saturation_ratio, args)
        else:
            raise ValueError("Invalid initialization type. Choose 'random' or 'geometric'.")

    def random_simulation_matrix(self, rock, wetting_saturation_ratio, args=None):
        erock = edist(rock)
        
        # Ensure all the BCs have bounce back nodes
        erock[0, :, :] = 1
        erock[:, 0, :] = 1
        erock[:, :, 0] = 1
        erock[-1, :, :] = 1
        erock[:, -1, :] = 1
        erock[:, :, -1] = 1
        
        # Re open the pores
        erock[rock == 0] = 0
        
        # Get the final matrix [0,1,2]
        erock[(erock > 0) & (erock < 2)] = 1
        erock[erock > 1] = 2
        
        # Get indices of all pore spaces
        pore_indices = np.argwhere(erock==0)
        
        # Determine number of each type of fluid to insert
        num_pores = len(pore_indices)
        num_wetting = int(wetting_saturation_ratio * num_pores)
        num_nonwetting = num_pores - num_wetting

        # Shuffle pore indices and separate into wetting and non-wetting
        np.random.shuffle(pore_indices)
        wetting_indices = pore_indices[:num_wetting]
        nonwetting_indices = pore_indices[num_wetting:]

        # Assign new values to wetting phase, boundary, inside solid, non-wetting phase
        erock = erock.astype(np.int16)
        erock[erock == 1] = 2609  # boundary
        erock[erock == 2] = 2610  # grains
        for index in wetting_indices:
            erock[tuple(index)] = 2608  # wetting fluid
        for index in nonwetting_indices:
            erock[tuple(index)] = 2611  # non-wetting fluid
        erock = erock.astype(np.int16)
        
        # Determine the geometry name based on the 'print_size' flag
        if args.print_size:
            size = erock.shape
            geom_name = f'{args.name}_{size[0]}_{size[1]}_{size[2]}'
        else:
            geom_name = args.name
        
        erock.flatten().tofile(f'{args.loc}/input/{geom_name}.dat')
        return erock

    def geometric_simulation_matrix(self, rock, wetting_saturation_ratio, args=None):
        erock = edist(rock)
        
        # Ensure all the BCs have bounce back nodes
        erock[0, :, :] = 1
        erock[:, 0, :] = 1
        erock[:, :, 0] = 1
        erock[-1, :, :] = 1
        erock[:, -1, :] = 1
        erock[:, :, -1] = 1
        
        # Re open the pores
        erock[rock == 0] = 0
        
        # Get the final matrix [0,1,2]
        erock[(erock > 0) & (erock < 2)] = 1
        erock[erock > 1] = 2
        
        pore_voxel_num = np.count_nonzero(erock==0)
        wetting_voxel_num = int(pore_voxel_num*wetting_saturation_ratio)
        non_wetting_voxel_num = pore_voxel_num - wetting_voxel_num 
        
        # Fill non-wetting ratio from the farthest points from the grains 
        inverse_edist = edist(np.where(rock==1, 0, 1))
        sorted_indices = np.unravel_index(np.argsort(-inverse_edist, axis=None), inverse_edist.shape)
        
        # Assuming 'arr_to_change' is the array you want to change
        for i in range(non_wetting_voxel_num):
            erock[sorted_indices[0][i], sorted_indices[1][i], sorted_indices[2][i]] = 3
        
        # Assign new values to wetting phase, boundary, inside solid, non-wetting phase
        erock = erock.astype(np.int16)
        erock[erock == 0] = 2608  # pore space / wetting fluid
        erock[erock == 1] = 2609  # boundary
        erock[erock == 2] = 2610  # grains
        erock[erock == 3] = 2611  # non-wetting fluid
        erock = erock.astype(np.int16)
        
        # Determine the geometry name based on the 'print_size' flag
        if args.print_size:
            size = erock.shape
            geom_name = f'{args.name}_{size[0]}_{size[1]}_{size[2]}'
        else:
            geom_name = args.name
        
        erock.flatten().tofile(f'{args.loc}/input/{geom_name}.dat')
        return erock

def create_geom_edist(rock, args):
    """
    Modifies a given rock matrix, calculates its Euclidean distance, creates a geometry file based on the 
    modified matrix and returns it.

    This function primarily modifies an input rock matrix according to several configuration options provided 
    in the 'args' parameter, calculates the Euclidean distance for the modified matrix, and writes it to a .dat 
    file. Additionally, the function ensures that all the boundary conditions have bounding box nodes, and it 
    pads the matrix if the 'num_slices' argument is provided.

    Parameters
    ----------
        rock : np.array
            The input rock matrix to be modified.
        args : argparse.Namespace
            An object containing various configuration options. This includes the 'swapXZ', 'scale_2', 'add_mesh', 
            'num_slices', 'print_size' and 'name' flags, and 'loc' directory path.

    Returns
    -------
        erock : np.array
            The modified rock matrix after the Euclidean distance calculation.

    Raises
    ------
        NotImplementedError
            If the 'scale_2' or 'add_mesh' features are attempted to be used as they are not yet implemented.

    Notes
    -----
        The function uses numpy's 'transpose', 'pad' and 'tofile' functions, as well as a custom 'edist' function 
        for Euclidean distance calculation.

    Examples
    --------
        >>> import numpy as np
        >>> import argparse
        >>> rock = np.array([[[0, 1, 0], [1, 0, 1], [0, 1, 0]], [[1, 0, 1], [0, 1, 0], [1, 0, 1]], [[0, 1, 0], [1, 0, 1], [0, 1, 0]]])
        >>> args = argparse.Namespace(swapXZ=False, scale_2=False, add_mesh=False, num_slices=2, print_size=True, name='geom', loc='./')
        >>> mod_rock = create_geom_edist(rock, args)
    """
    # If the 'swapXZ' flag is set, transpose the rock matrix accordingly
    if args.swapXZ:
        rock = rock.transpose([2, 1, 0])

    # If the 'scale_2' flag is set, raise an error as this feature is not yet implemented
    if args.scale_2:
        raise NotImplementedError('Feature not yet implemented')

    # Calculate the Euclidean distance for the rock matrix
    erock = edist(rock)
    
    # Ensure all the BCs have BB nodes
    erock[0,:,:] = erock[:,0,:] = erock[:,:,0] = 1
    erock[-1,:,:] = erock[:,-1,:] = erock[:,:,-1] = 1
    
    # Reopen the pores
    erock[rock==0] = 0
    
    # Get the final matrix with values [0,1,2]
    erock[(erock > 0) * (erock < 2)] = 1
    erock[erock > 1] = 2
    
    # If the 'add_mesh' flag is set, raise an error as this feature is not yet implemented
    if args.add_mesh:
        raise NotImplementedError('Feature not yet implemented')
    
    # If the 'num_slices' argument is provided, pad the 'erock' array accordingly
    if args.num_slices:
        erock = np.pad(erock, [(args.num_slices, args.num_slices), (0, 0), (0, 0)])
    
    # Determine the geometry name based on the 'print_size' flag
    if args.print_size:
        size = erock.shape
        geom_name = f'{args.name}_{size[0]}_{size[1]}_{size[2]}'
    else:
        geom_name = args.name
    
    # Modify the 'erock' array's data type and values for final output
    erock = erock.astype(np.int16)
    erock[erock==0] = 2608
    erock[erock==1] = 2609
    erock[erock==2] = 2610
    
    # Write the 'erock' array to a .dat file
    erock.flatten().tofile(f'{args.loc}/input/{geom_name}.dat')
    
    return erock




def erase_regions(rock):
    """
    Identifies and erases isolated regions within a given rock matrix.

    The function employs the `label` function from the skimage.measure module to detect connected components 
    in the input rock matrix. Isolated regions are identified as those components which are not connected to 
    the main body of the matrix (where the corresponding element in the labeled matrix is greater than 1). 
    These regions are then removed (set to 0) from the rock matrix.

    Parameters
    ----------
        rock : np.array
            The input rock matrix to be processed.

    Returns
    -------
        rock : np.array
            The modified rock matrix with isolated regions removed.

    Notes
    -----
        The function uses skimage.measure's 'label' function for connected components detection. The background 
        is set to 1, and connectivity is set to 1 to consider only orthogonally adjacent points. 

    Examples
    --------
    >>> import numpy as np
    >>> rock = np.array([[[0, 1, 0], [1, 0, 1], [0, 1, 0]], [[1, 0, 1], [0, 1, 0], [1, 0, 1]], [[0, 1, 0], [1, 0, 1], [0, 1, 0]]])
    >>> mod_rock = erase_regions(rock)
    """
    # Use the measure.label function from the skimage.measure module to identify connected components in the rock matrix
    # Here, the background is set to 1, and connectivity is set to 1 to consider only orthogonally adjacent points
    blobs_labels = measure.label(rock, background=1, connectivity=1)
    
    # Erase all isolated regions within the rock matrix
    # This is done by setting all elements in 'rock' that are part of an isolated region (where the corresponding element in 'blobs_labels' is greater than 1) to 0
    rock[blobs_labels > 1] = 1
    
    return rock




class write_MPLBM():
    """
    This class encapsulates methods to set up, manage and execute an MPLBM (MultiPhase Lattice Boltzmann Method) 
    simulation. It uses provided or default settings to create the necessary files and directories, write geometry,
    write inputs, and execute the simulation. 
    """

    def __init__(self, frac_obj, buffer_layers=2, cpus=4, num_hrs=14, allocation=None, sim_mode='single-phase', fluid_ini_conf='geometric', wetting_f_sat=0.5):
        """
        Initializes the write_MPLBM class
        
        Parameters
        ----------
            frac_obj : object
                An object containing a fracture object.
            buffer_layers : int, optional
                The number of empty slices at the beginning and end of the domain for pressure boundary conditions. 
                Defaults to 2.
            cpus : int, optional
                The number of CPUs to be used in the simulation. Defaults to 4.
            num_hrs : int, optional
                The number of hours the simulation should run. Defaults to 14.
            allocation : str, optional
                A specific allocation for the simulation. Defaults to None.

        Returns
        --------
            None
        
        Notes
        -----
            The class initializes by defining the source location of LBM (Lattice Boltzmann Method), creating required
            directories, copying the surface fraction object, and pre-processing the 3D fracture to generate an efficient 
            geometry for simulation. It also writes a shell script to run the simulation.

            The class includes methods to create new folders, write geometry, write various input files, copy necessary 
            data, edit inputs, and save the configuration as a pickle file for future reference. 
        """ 
        # Location of LBM source file
        if sim_mode == 'single-phase':
            lbm_loc = '../../lbm_light/src/1-phase_LBM/permeability' 
        elif sim_mode == 'two-phase':
            lbm_loc = '../../MPLBM-UT/src/2-phase_LBM/ShanChen'
                    
        
        # Define geometry namespace
        geom = Namespace()
        name = 'frac'
        geom.name = name
        geom.print_size = True
        geom.add_mesh = False
        geom.num_slices = buffer_layers
        geom.swapXZ = True              
        geom.scale_2 = False                  
        
        # Create required directories
        create_folder('sims')
        self.folder_num = self.create_new_folder()
        self.folder_path = f'sims/frac_{self.folder_num}'
        create_folder(f'{self.folder_path}/input')
        
        # Assign folder path to geom.loc
        geom.loc = self.folder_path
        
        # Copy the surface fraction object and remove 'frac_3D' to reduce pickle size
        self.my_frac = cp(frac_obj)
        del self.my_frac.frac_3D
        
        self.buffer_layers = buffer_layers
        self.cpus = cpus
        self.num_hrs = num_hrs
        self.allocation = allocation
    
        # Preprocess the 3D fracture and generate an efficient geometry for simulation
        frac_3D = 1 - np.transpose(frac_obj.frac_3D, (2, 0, 1))
        frac_3D = erase_regions(frac_3D)


        if sim_mode == 'single-phase':
            frac_3D = create_geom_edist(frac_3D, geom)
            # Create input file
            frac_size = frac_3D.shape
            input_file_name = f"{self.folder_path}/{name}.xml"
            geom_name = f'{geom.name}_{frac_size[0]}_{frac_size[1]}_{frac_size[2]}'
            domain_size = [frac_size[0], frac_size[1], frac_size[2]]
            periodic = ["false", "false", "false"]
            
            create_folder(f'{self.folder_path}/output')
            
            io_folders = ['input/', 'output/']
            num_sims = 1
            sim_pressure = 0.0005
            sim_max_iter = 1000000
            sim_convergence = 0.0001
            sim_settings = [num_sims, sim_pressure, sim_max_iter, sim_convergence]
            save_vtks = "true"
            
            create_single_phase_input_file(input_file_name, geom_name, domain_size, periodic, io_folders, sim_settings, save_vtks)

        elif sim_mode == 'two-phase':
            sim_init = SimulationInitializer(fluid_ini_conf)  # or 'geometric'
            frac_3D = sim_init.initialize_simulation_matrix(frac_3D, wetting_f_sat, geom)
            
            frac_size = frac_3D.shape
            input_file_name = f"{self.folder_path}/{name}.xml"
            geom_name = f'{geom.name}_{frac_size[0]}_{frac_size[1]}_{frac_size[2]}'
            domain_size = [frac_size[0], frac_size[1], frac_size[2]]
            periodic = ["True", "True", "True"]
            create_folder(f'{self.folder_path}/output')
            
            io_folders = ['input/', 'output/']


            # sim settings
            sim_pressure = 0.0005
            sim_max_iter = 500000
            sim_convergence = 1e-4
            convergence_relperm = 1e-6
            max_iterations_relperm = 10**5
            
            fluid_1_init = [1, 2, 1, 75, 1, 75]
            fluid_2_init = [3, 150, 1, 75, 1, 75]
            Gc = 0.9
            omega_f1 = 1
            omega_f2 = 1
            G_ads_f1_s1 =-0.4
            G_ads_f1_s2 = 0
            G_ads_f1_s3 = 0
            G_ads_f1_s4 = 0

            convergence_iter = 1000

            save_sim = True
            save_iter = 20000
            gif_iter = 2000
            vtk_iter = 2000
            rho_f2_vtk = False
            print_geom = True
            print_stl = False
            save_vtks =True

            pressure_bc = False
            minimum_radius = 3
            num_pressure_steps = 1
            force_f1 = 1e-4
            force_f2 = 1e-4
            rho_f1 = 2
            rho_f2 = 2   

            sim_settings = [sim_max_iter, sim_convergence, convergence_relperm,
                            max_iterations_relperm, fluid_1_init, fluid_2_init, Gc, omega_f1, omega_f2,
                            G_ads_f1_s1, G_ads_f1_s2, G_ads_f1_s3, G_ads_f1_s4, convergence_iter, convergence_iter,
                            save_sim, save_iter, gif_iter, vtk_iter, rho_f2_vtk, rho_f2_vtk, print_geom, print_stl,
                            save_vtks, pressure_bc, minimum_radius, num_pressure_steps, force_f1, force_f2,
                            rho_f1, rho_f2]
            
            create_two_phase_input_file(input_file_name, geom_name, domain_size, periodic, io_folders, sim_settings)
                 
        else:
            raise ValueError(f'Invalid simulation_mode: {sim_mode}')

        # Write shell script to run the simulation
        np.savetxt(f'{self.folder_path}/run_{name}.sh', [f'mpirun -np {cpus} {lbm_loc} {name}.xml > {name}.txt'], fmt='%s') 
        print('write LBM end')    


    def create_new_folder(self):
        """
        Creates a new simulation directory in the 'sims' folder. It searches through numbered folders 
        'frac_i' (i ranging from 0 to 9999999999) and creates the first one that doesn't already exist.

        Parameters
        ----------
            None

        Returns
        -------
            i : int
                The number 'i' of the newly created folder.

        Notes
        ------
            None
        
        """
        
        for i in range(10000000000):
            try:
                mkdir(f'sims/frac_{i}')
                break
            except FileExistsError:
                continue
        return i
    

    def write_geom(self, frac):
        """
        Writes the 3D fracture geometry to a .dat file.

        The function first appends the shape of the input 'frac' array to the flattened 'frac' array, 
        then saves this information to a .dat file in the current simulation's directory.

        Parameters
        ----------
            frac : np.array
                The 3D fracture geometry to be written.

        Returns
        -------
            None

        Notes
        ------
            None
        
        """

        tofile = np.append([*frac.shape], frac.flatten('F'))
        np.savetxt(f'{self.folder_path}/frac.dat', tofile, fmt='%.0f')

        


    def save_pickle(self):
        """
        Saves the current instance of the write_MPLBM class as a pickled object.
        
        Parameters
        -----------
            None

        Returns
        -------
            None
            
        Raises
        ------
            PicklingError
                If the pickling process encounters an object that it cannot serialize.

        Notes
        ------
            None
        """
    

        with open(f'{self.folder_path}/lbm_config.pk', 'wb') as f:
            pickle.dump(self, f, protocol = pickle.HIGHEST_PROTOCOL) 
               


