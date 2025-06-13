import numpy as np 
import subprocess
import shutil 
import os
import glob 

def write_lagrit_script_convert_top_and_bottom_to_avs(self):
    """
    Generate a LaGriT input script to create top and bottom surfaces,
    and convert them to STL files.

    The script creates multiple computational meshes:
    - `mo_bottom`, `mo_top`, and `mo_quad` surface meshes
    - `mo_tri` for tetrahedral conversion
    - Dumps include `top.stl`, `bottom.stl`, and intermediate `.inp` files

    The generated script is saved as `convert_fracture_to_avs.lgi`.

    Assumes the instance has attributes:
        - `self.nx` (int): Number of grid cells in X
        - `self.ny` (int): Number of grid cells in Y
        - `self.lx` (float): Length in X direction
        - `self.ly` (float): Length in Y direction
    """

    lagrit_script = f"""

define / NX / {self.nx}
define / NY / {self.ny}

define / X1 / {self.lx:0.2f} 
define / Y1 / {self.ly:0.2f}  

# cmo/create/ mo_bottom / / / quad
# quadxy/ NX, NY / 0. 0. -1 / X1 0. -1 / X1 Y1 -1 / 0. Y1 -1
# createpts/brick/xyz/NX,NY,1/1 0 0 / connect

# cmo/create/ mo_top / / / quad
# quadxy/ NX, NY / 0. 0. 1 / X1 0. 1 / X1 Y1 1 / 0. Y1 1
# createpts/brick/xyz/NX,NY,1/1 0 0 / connect

# cmo / delete / mo_top
# cmo / delete / mo_bottom

cmo/create/ mo_quad / / / quad
quadxy/ NX, NY / 0. 0. 0. / X1 0. 0. / X1 Y1 0. / 0. Y1 0.
createpts/brick/xyz/NX,NY,1/1 0 0 / connect

resetpts / itp
cmo / status / brief

cmo / setatt / mo_quad / itetclr / 1 0 0 / 1
cmo / readatt / mo_quad / height / 1,0,0 / top.dat

resetpts / itp
cmo / status / brief
hextotet / 2 / mo_tri  / mo_quad
cmo / status / brief
cmo / copyatt / mo_tri / mo_tri / zic  / height 
dump / surf_top.inp / mo_tri
cmo / status / brief
cmo / delete / mo_quad
cmo / select / mo_tri 
dump / stl / top.stl


// cmo / delete / mo_tri
## load Bottom

cmo/create/ mo_quad / / / quad
quadxy/ NX, NY / 0. 0. 0. / X1 0. 0. / X1 Y1 0. / 0. Y1 0.
createpts/brick/xyz/NX,NY,1/1 0 0 / connect
cmo / setatt / mo_quad / itetclr / 1 0 0 / 1
cmo / readatt / mo_quad / height / 1,0,0 / bottom.dat

resetpts / itp
cmo / status / brief

hextotet / 2 / mo_tri  / mo_quad

cmo / status / brief

// dump / tri_with_top.inp / mo_tri

cmo / copyatt / mo_tri / mo_tri / zic  / height 

dump / surf_bottom.inp / mo_tri

cmo / delete / mo_quad

dump / stl / bottom.stl
cmo / delete / mo_tri

finish 

"""

    with open('convert_fracture_to_avs.lgi', 'w') as fp:
        fp.write(lagrit_script)
        fp.flush() 


def write_lagrit_script_interior_stl():
    """
    Generate a LaGriT input script for stacking and extracting STL surfaces
    from filled geometry between top and bottom surfaces.

    Operations include:
    - Stack creation and hex mesh filling
    - Extraction of surface meshes by `itetclr` labels
    - Export of STL files: `inlet.stl`, `outlet.stl`, `left.stl`, `right.stl`

    Saves the script content into a hardcoded output file, if added later.
    """

    lagrit_script = """
cmo create cmo_stack

stack / layers/avs/ surf_bottom.inp / surf_top.inp

// dump / stacked.inp / cmo_stack

stack / fill / cmohex / cmo_stack 

// dump / stacked.inp / cmohex

resetpts / itp
cmo / status / brief

cmo / setatt / cmohex / itetclr / 1 0 0 / 1

extract/surfmesh/ 1,0,0 / mo_surf / cmohex

// dump / external.inp / mo_surf

cmo / select / mo_surf
settets/normal

resetpts / itp
cmo / status / brief

// dump / external.inp / mo_surf

cmo / select / mo_surf 
eltset/ e_delete/ idface1 eq 1
rmpoint / element / eltset get e_delete
rmpoint / compress


cmo / select / mo_surf 
eltset/ e_delete/ idface1 eq 2
rmpoint / element / eltset get e_delete
rmpoint / compress


// dump / external.inp / mo_surf


## Grab surfaces and dump stl files
cmo / copy/ mo_inlet /mo_surf
cmo / select / mo_inlet 
eltset/ e_delete/ itetclr ne 4
rmpoint / element / eltset get e_delete
rmpoint / compress
cmo / status / brief 
hextotet / 2 / mo_tri  / mo_inlet
cmo / select / mo_tri 
// dump / inlet.inp / mo_tri
dump / stl / inlet.stl / mo_tri 
cmo / delete / mo_tri
cmo / delete / mo_inlet 

cmo / copy/ mo_inlet /mo_surf
cmo / select / mo_inlet 
eltset/ e_delete/ itetclr ne 6
rmpoint / element / eltset get e_delete
rmpoint / compress
cmo / status / brief 
hextotet / 2 / mo_tri  / mo_inlet
cmo / select / mo_tri 
// dump / outlet.inp / mo_tri
dump / stl / outlet.stl / mo_tri
cmo / delete / mo_tri
cmo / delete / mo_inlet

cmo / copy/ mo_inlet /mo_surf
cmo / select / mo_inlet 
eltset/ e_delete/ itetclr ne 3
rmpoint / element / eltset get e_delete
rmpoint / compress
cmo / status / brief 
hextotet / 2 / mo_tri  / mo_inlet
cmo / select / mo_tri 
// dump / left.inp / mo_tri
dump / stl / left.stl / mo_tri
cmo / delete / mo_tri
cmo / delete / mo_inlet


cmo / copy/ mo_inlet /mo_surf
cmo / select / mo_inlet 
eltset/ e_delete/ itetclr ne 5
rmpoint / element / eltset get e_delete
rmpoint / compress
cmo / status / brief 
hextotet / 2 / mo_tri  / mo_inlet
cmo / select / mo_tri 
// dump / right.inp / mo_tri
dump / stl / right.stl / mo_tri
cmo / delete / mo_tri
cmo / delete / mo_inlet

finish 

"""

    with open('convert_interior_to_stl.lgi', 'w') as fp:
        fp.write(lagrit_script)
        fp.flush() 


def write_lagrit_script_extract_exterior(self):
    """
    Generate a LaGriT script to extract and export STL surfaces for the
    top, bottom, and full external geometry of a 3D stacked mesh block.

    This function computes adjusted top and bottom surface elevations based on
    the `top` and `bottom` height arrays, then constructs the following:

    - STL files and AVS `.inp` mesh files for:
        - The bottom external surface
        - The top external surface
        - The full external mesh block
    - Volume stacking using `stack / fill` and `extract / surfmesh`
    - Mesh processing steps like `hextotet`, `geniee`, and clean-up of CMO objects

    The generated LaGriT script is saved as `make_exterior_blocks.lgi`.

    Assumes the following instance attributes are defined:
        self.nx (int): Number of grid points in the X direction.
        self.ny (int): Number of grid points in the Y direction.
        self.lx (float): Physical length of the domain in the X direction.
        self.ly (float): Physical length of the domain in the Y direction.
        self.top (np.ndarray): 2D or 1D array representing top surface elevation.
        self.bottom (np.ndarray): 2D or 1D array representing bottom surface elevation.

    Output Files:
        - bottom_external.stl, top_external.stl, whole_block.stl
        - bottom_external.inp, top_external.inp, whole_block_external.inp
        - blocks_empty.inp

    Raises:
        ValueError: If the shape of `top` or `bottom` does not match (nx * ny).

    Side Effects:
        Writes a LaGriT script file `make_exterior_blocks.lgi` to the current directory.
    """

    top = np.reshape(self.top, self.nx*self.ny)   
    bottom = np.reshape(self.bottom, self.nx*self.ny)   
    top_surface_height = max(max(top),max(bottom))*1.2
    bottom_surface_height = min(min(bottom),min(top))*1.2

    lagrit_script = f"""
define / NX / {self.nx}
define / NY / {self.ny}

define / X1 / {self.lx:0.2f} 
define / Y1 / {self.ly:0.2f}  

## make bottom of block
cmo/create/ mo_quad / / / quad
quadxy/ NX, NY / 0. 0. {bottom_surface_height:0.2f} / X1 0.  {bottom_surface_height:0.2f} & 
     / X1 Y1  {bottom_surface_height:0.2f} / 0. Y1  {bottom_surface_height:0.2f}
createpts/brick/xyz/NX,NY,1/1 0 0 / connect
cmo / setatt / mo_quad / itetclr / 1 0 0 / 1
hextotet / 2 / mo_tri  / mo_quad
dump / bottom.inp / mo_tri

cmo / delete / mo_tri
cmo / delete / mo_quad 

## Make top of block
cmo/create/ mo_quad / / / quad
quadxy/ NX, NY / 0. 0. {top_surface_height:0.2f} / X1 0. {top_surface_height:0.2f} &
    / X1 Y1 {top_surface_height:0.2f} / 0. Y1 {top_surface_height:0.2f}
createpts/brick/xyz/NX,NY,1/1 0 0 / connect
cmo / setatt / mo_quad / itetclr / 1 0 0 / 1
hextotet / 2 / mo_tri  / mo_quad

dump / top.inp / mo_tri

# Make stacks and dump exterior of the bottom
cmo create cmo_stack
stack / layers / avs / bottom.inp 1 / & 
     surf_bottom.inp 2 

# dump / stacked.inp / cmo_stack
stack / fill / cmohex / cmo_stack 
# dump / stacked_hex_bottom.inp / cmohex
resetpts / itp
cmo / status / brief
extract / surfmesh / 1,0,0 / mo_surf / cmohex
cmo / select / mo_surf
dump / bottom_external.inp / mo_surf

hextotet / 2 / mo_tri  / mo_surf
cmo / select / mo_tri 
geniee 
dump / stl / bottom_external.stl / mo_tri

# Clean up 
cmo  /delete / mo_tri
cmo  /delete / mo_surf
cmo / delete / cmohex 
cmo / delete / cmo_stack 

# Make stacks and dump exterior of the top block


cmo create cmo_stack
stack / layers / avs / surf_top.inp 1 / & 
     top.inp 2 

# dump / stacked.inp / cmo_stack
stack / fill / cmohex / cmo_stack 
#dump / stacked_hex_top.inp / cmohex
resetpts / itp
cmo / status / brief
extract / surfmesh / 1,0,0 / mo_surf / cmohex
cmo / select / mo_surf
dump / top_external.inp / mo_surf


hextotet / 2 / mo_tri  / mo_surf
cmo / select / mo_tri 
geniee 
dump / stl / top_external.stl / mo_tri


# Clean up 
cmo  /delete / mo_tri
cmo  /delete / mo_surf
cmo / delete / cmohex 
cmo / delete / cmo_stack 

## dump whole block to do 
cmo create cmo_stack
stack / layers / avs /  bottom.inp 1 / & 
     surf_bottom.inp 2 surf_top.inp 3 / & 
     top.inp 4
# dump / stacked.inp / cmo_stack
stack / fill / cmohex / cmo_stack 
# dump / blocks.inp / cmohex
resetpts / itp
cmo / status / brief

cmo / select / cmohex 
eltset/ e_delete/ itetclr eq 2
rmpoint / element / eltset get e_delete
rmpoint / compress

dump / blocks_empty.inp / cmohex
resetpts / itp
cmo / status / brief

resetpts / itp
cmo / status / brief
extract / surfmesh / 1,0,0 / mo_surf / cmohex
cmo / select / mo_surf
dump / whole_block_external.inp / mo_surf

hextotet / 2 / mo_tri  / mo_surf
cmo / select / mo_tri 
geniee 
dump / stl / whole_block.stl / mo_tri


finish 

"""
    with open('make_exterior_blocks.lgi', 'w') as fp:
        fp.write(lagrit_script)
        fp.flush() 

def cleanup():
    """
    Clean up the current working directory by organizing and removing specific files.

    Operations performed:
    1. Move all `.stl` files into a directory named `stl_files`.
    2. Move all `.lgi` files into a directory named `lagrit_scripts`.
    3. Delete all `.inp` and `.dat` files from the current directory.

    Directories (`stl_files` and `lagrit_scripts`) are created fresh. If they exist,
    they are removed and recreated to ensure they only contain current-session files.

    Raises:
        OSError: If a file or directory operation fails unexpectedly.
    """

    def move_files(pattern, target_dir):
        """
        Move files matching a pattern to a specified directory.

        Args:
            pattern (str): Glob pattern to match files (e.g., '*.stl').
            target_dir (str): Directory where matched files should be moved.
        """
        files_to_move = glob.glob(pattern)
        if os.path.exists(target_dir):
            try:
                shutil.rmtree(target_dir)
            except Exception as e:
                print(f"Error removing directory '{target_dir}': {e}")
                return
        try:
            os.mkdir(target_dir)
        except Exception as e:
            print(f"Error creating directory '{target_dir}': {e}")
            return

        for file in files_to_move:
            dest = os.path.join(target_dir, file)
            try:
                os.rename(file, dest)
            except Exception as e:
                print(f"Error moving file '{file}' to '{target_dir}': {e}")

    # Move STL and LGI files
    move_files("*.stl", "stl_files")
    move_files("*.lgi", "lagrit_scripts")

    # Remove .inp and .dat files
    for pattern in ["*.inp", "*.dat"]:
        for file in glob.glob(pattern):
            try:
                os.remove(file)
            except Exception as e:
                print(f"Error deleting file '{file}': {e}")



def _run_lagrit(self, script_file):
    """
    Run a LaGriT script with error checking.

    Args:
        script_file (str): Path to the LaGriT input script.

    Raises:
        FileNotFoundError: If the script file does not exist.
        RuntimeError: If LaGriT execution fails.
    """
    if not os.path.exists(script_file):
        raise FileNotFoundError(f"LaGriT script '{script_file}' not found.")
    try:
        subprocess.run(f"lagrit < {script_file}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"LaGriT execution failed on script: {script_file}") from e

def dump_stl(self, interior=True, exterior=False, cleanup_run=True):
    """
    Generate STL surface files representing a fracture mesh using LaGriT and optionally extract
    interior and/or exterior geometry.

    This method performs the following steps:
    1. Flattens and saves `top` and `bottom` surface data to `top.dat` and `bottom.dat`.
    2. Executes LaGriT scripts to convert those surfaces to STL format.
    3. Optionally extracts and exports:
       - Interior surfaces between the top and bottom using `convert_interior_to_stl.lgi`.
       - Exterior mesh geometry and block structure using `make_exterior_blocks.lgi`.
    4. Optionally cleans up intermediate files from the current directory.

    Args:
        interior (bool, optional): 
            If True, generates internal surface STL files 
            (e.g., `inlet.stl`, `outlet.stl`, etc.). Defaults to True.

        exterior (bool, optional): 
            If True, generates exterior STL representations of the block, 
            including top, bottom, and full volume. Defaults to False.

        cleanup_run (bool, optional): 
            If True, removes intermediate `.dat` and `.inp` files 
            and moves `.stl` files into `stl_files/`. Defaults to True.

    Raises:
        ValueError: If the `top` or `bottom` surface arrays cannot be reshaped to `nx * ny`.
        FileNotFoundError: If required LaGriT executable or inputs are missing.
        RuntimeError: If a LaGriT subprocess call fails.

    Side Effects:
        - Creates `top.dat` and `bottom.dat` files from reshaped surface arrays.
        - Generates and executes `.lgi` scripts using subprocess calls.
        - Produces multiple `.stl` and `.inp` files in the current working directory.
        - Optionally cleans up and organizes generated files via `cleanup()`.

    Note:
        Requires the LaGriT binary (`lagrit`) to be installed and available in the system's PATH.

    Example:
        >>> myfrac.dump_stl(interior=True, exterior=True)
    """

    print("--> Dumping fracture surfaces to STL")
    top = np.reshape(self.top, self.nx*self.ny)   
    np.savetxt("top.dat", top, delimiter=" ")
    bottom = np.reshape(self.bottom, self.nx*self.ny)   
    np.savetxt("bottom.dat", bottom, delimiter=" ")

    self.write_lagrit_script_convert_top_and_bottom_to_avs()
    self._run_lagrit('convert_fracture_to_avs.lgi')

    if interior:
        write_lagrit_script_interior_stl()
        self._run_lagrit('convert_interior_to_stl.lgi')

    if exterior:
        self.write_lagrit_script_extract_exterior()
        self._run_lagrit('make_exterior_blocks.lgi')

    if cleanup_run:
        cleanup()
    print("--> Dumping fracture surfaces to STL Complete")