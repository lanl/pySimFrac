import numpy as np 
import subprocess
import shutil 
import os
import glob 
# from stl import mesh


def write_lagrit_script_part1(self):

    lagrit_script = f"""

define / NX / {self.nx}
define / NY / {self.ny}

define / X1 / {self.lx:0.2f} 
define / Y1 / {self.ly:0.2f}  

cmo/create/ mo_bottom / / / quad
quadxy/ NX, NY / 0. 0. -1 / X1 0. -1 / X1 Y1 -1 / 0. Y1 -1
createpts/brick/xyz/NX,NY,1/1 0 0 / connect

// dump / bottom.inp / mo_bottom

cmo/create/ mo_top / / / quad
quadxy/ NX, NY / 0. 0. 1 / X1 0. 1 / X1 Y1 1 / 0. Y1 1
createpts/brick/xyz/NX,NY,1/1 0 0 / connect

// dump / top.inp / mo_top

cmo / delete / mo_top
cmo / delete / mo_bottom

cmo/create/ mo_quad / / / quad
quadxy/ NX, NY / 0. 0. 0. / X1 0. 0. / X1 Y1 0. / 0. Y1 0.
createpts/brick/xyz/NX,NY,1/1 0 0 / connect

// dump / quad.inp / mo_quad

resetpts / itp
cmo / status / brief

cmo / setatt / mo_quad / itetclr / 1 0 0 / 1
cmo / readatt / mo_quad / height / 1,0,0 / top.dat

resetpts / itp
cmo / status / brief

hextotet / 2 / mo_tri  / mo_quad

cmo / status / brief

// dump / tri_with_top.inp / mo_tri

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


# def write_lagrit_script_convert_exterior_to_stl(top_surface_height, bottom_surface_height):

#     lagrit_script = """

# read / top.inp / mo_top
# hextotet / 2 / mo_tri  / mo_top
# cmo / select / mo_tri 
# dump / stl / top_opposite.stl / mo_tri 
# cmo / delete / mo_tri
# cmo / delete / mo_top 

# read / surf_top.inp / mo_top
# hextotet / 2 / mo_tri  / mo_top
# cmo / select / mo_tri 
# dump / stl / top_surf.stl / mo_tri 
# cmo / delete / mo_tri
# cmo / delete / mo_top 

# read / top_external.inp / mo_surf 

# cmo / select / mo_surf
# settets/normal

# cmo / select / mo_surf 
# eltset/ e_delete/ idface1 eq 1
# rmpoint / element / eltset get e_delete
# rmpoint / compress

# cmo / select / mo_surf 
# eltset/ e_delete/ idface1 eq 2
# rmpoint / element / eltset get e_delete
# rmpoint / compress

# ## Grab surfaces and dump stl files
# # cmo / copy/ mo_inlet /mo_surf
# # cmo / select / mo_inlet 
# # eltset/ e_delete/ itetclr ne 1
# # rmpoint / element / eltset get e_delete
# # rmpoint / compress
# # cmo / status / brief 
# # hextotet / 2 / mo_tri  / mo_inlet
# # cmo / select / mo_tri 
# # // dump / bottom.inp / mo_tri
# # dump / stl / top_opposite.stl / mo_tri 
# # cmo / delete / mo_tri
# # cmo / delete / mo_inlet 

# cmo / copy/ mo_inlet /mo_surf
# cmo / select / mo_inlet 
# eltset/ e_delete/ itetclr ne 4
# rmpoint / element / eltset get e_delete
# rmpoint / compress
# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_inlet
# cmo / select / mo_tri 
# // dump / inlet.inp / mo_tri
# dump / stl / top_front.stl / mo_tri 
# cmo / delete / mo_tri
# cmo / delete / mo_inlet 

# cmo / copy/ mo_inlet /mo_surf
# cmo / select / mo_inlet 
# eltset/ e_delete/ itetclr ne 6
# rmpoint / element / eltset get e_delete
# rmpoint / compress
# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_inlet
# cmo / select / mo_tri 
# // dump / outlet.inp / mo_tri
# dump / stl / top_back.stl / mo_tri
# cmo / delete / mo_tri
# cmo / delete / mo_inlet

# cmo / copy/ mo_inlet /mo_surf
# cmo / select / mo_inlet 
# eltset/ e_delete/ itetclr ne 3
# rmpoint / element / eltset get e_delete
# rmpoint / compress
# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_inlet
# cmo / select / mo_tri 
# // dump / left.inp / mo_tri
# dump / stl / top_left.stl / mo_tri
# cmo / delete / mo_tri
# cmo / delete / mo_inlet


# cmo / copy/ mo_inlet /mo_surf
# cmo / select / mo_inlet 
# eltset/ e_delete/ itetclr ne 5
# rmpoint / element / eltset get e_delete
# rmpoint / compress
# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_inlet
# cmo / select / mo_tri 
# // dump / right.inp / mo_tri
# dump / stl / top_right.stl / mo_tri
# cmo / delete / mo_tri
# cmo / delete / mo_inlet

# finish 

# """
#     with open('top_exterior_to_stl.lgi', 'w') as fp:
#         fp.write(lagrit_script)
#         fp.flush() 


#     lagrit_script = """

# read / bottom.inp / mo_bottom
# hextotet / 2 / mo_tri  / mo_bottom
# cmo / select / mo_tri 
# dump / stl / bottom_opposite.stl / mo_tri 
# cmo / delete / mo_tri
# cmo / delete / mo_top 

# read / surf_bottom.inp / mo_top
# hextotet / 2 / mo_tri  / mo_top
# cmo / select / mo_tri 
# dump / stl / bottom_surf.stl / mo_tri 
# cmo / delete / mo_tri
# cmo / delete / mo_top 

# read / bottom_external.inp / mo_surf 

# cmo / select / mo_surf
# settets/normal

# cmo / select / mo_surf 
# eltset/ e_delete/ idface1 eq 1
# rmpoint / element / eltset get e_delete
# rmpoint / compress

# cmo / select / mo_surf 
# eltset/ e_delete/ idface1 eq 2
# rmpoint / element / eltset get e_delete
# rmpoint / compress

# ## Grab surfaces and dump stl files
# # cmo / copy/ mo_inlet /mo_surf
# # cmo / select / mo_inlet 
# # eltset/ e_delete/ itetclr ne 1
# # rmpoint / element / eltset get e_delete
# # rmpoint / compress
# # cmo / status / brief 
# # hextotet / 2 / mo_tri  / mo_inlet
# # cmo / select / mo_tri 
# # // dump / bottom.inp / mo_tri
# # dump / stl / bottom.stl / mo_tri 
# # cmo / delete / mo_tri
# # cmo / delete / mo_inlet 

# cmo / copy/ mo_inlet /mo_surf
# cmo / select / mo_inlet 
# eltset/ e_delete/ itetclr ne 4
# rmpoint / element / eltset get e_delete
# rmpoint / compress
# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_inlet
# cmo / select / mo_tri 
# // dump / inlet.inp / mo_tri
# dump / stl / bottom_front.stl / mo_tri 
# cmo / delete / mo_tri
# cmo / delete / mo_inlet 

# cmo / copy/ mo_inlet /mo_surf
# cmo / select / mo_inlet 
# eltset/ e_delete/ itetclr ne 6
# rmpoint / element / eltset get e_delete
# rmpoint / compress
# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_inlet
# cmo / select / mo_tri 
# // dump / outlet.inp / mo_tri
# dump / stl / bottom_back.stl / mo_tri
# cmo / delete / mo_tri
# cmo / delete / mo_inlet

# cmo / copy/ mo_inlet /mo_surf
# cmo / select / mo_inlet 
# eltset/ e_delete/ itetclr ne 3
# rmpoint / element / eltset get e_delete
# rmpoint / compress
# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_inlet
# cmo / select / mo_tri 
# // dump / left.inp / mo_tri
# dump / stl / bottom_left.stl / mo_tri
# cmo / delete / mo_tri
# cmo / delete / mo_inlet


# cmo / copy/ mo_inlet /mo_surf
# cmo / select / mo_inlet 
# eltset/ e_delete/ itetclr ne 5
# rmpoint / element / eltset get e_delete
# rmpoint / compress
# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_inlet
# cmo / select / mo_tri 
# // dump / right.inp / mo_tri
# dump / stl / bottom_right.stl / mo_tri
# cmo / delete / mo_tri
# cmo / delete / mo_inlet

# finish 

# """
#     with open('bottom_exterior_to_stl.lgi', 'w') as fp:
#         fp.write(lagrit_script)
#         fp.flush() 

#     #whole block    
#     lagrit_script = """
# read / bottom.inp / mo_bottom
# hextotet / 2 / mo_tri  / mo_bottom
# cmo / select / mo_tri 
# dump / stl / whole_block_bottom.stl / mo_tri 
# cmo / delete / mo_tri
# cmo / delete / mo_bottom 

# read / top.inp / mo_top
# hextotet / 2 / mo_tri  / mo_top
# cmo / select / mo_tri 
# dump / stl / whole_block_top.stl / mo_tri 
# cmo / delete / mo_tri
# cmo / delete / mo_top 

# read / whole_block_external.inp / mo_surf 

# cmo / select / mo_surf
# settets/normal
# ## remove surfaces for normal vector based dump
# cmo / select / mo_surf 
# eltset/ e_delete/ idface1 eq 1
# rmpoint / element / eltset get e_delete
# rmpoint / compress

# cmo / select / mo_surf 
# eltset/ e_delete/ idface1 eq 2
# rmpoint / element / eltset get e_delete
# rmpoint / compress

# ## Grab surfaces and dump stl files

# cmo / copy/ mo_inlet /mo_surf
# cmo / select / mo_inlet 
# eltset/ e_delete/ itetclr ne 4
# rmpoint / element / eltset get e_delete
# rmpoint / compress
# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_inlet
# cmo / select / mo_tri 
# // dump / inlet.inp / mo_tri
# dump / stl / whole_block_front.stl / mo_tri 
# cmo / delete / mo_tri
# cmo / delete / mo_inlet 

# cmo / copy/ mo_inlet /mo_surf
# cmo / select / mo_inlet 
# eltset/ e_delete/ itetclr ne 6
# rmpoint / element / eltset get e_delete
# rmpoint / compress
# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_inlet
# cmo / select / mo_tri 
# // dump / outlet.inp / mo_tri
# dump / stl / whole_block_back.stl / mo_tri
# cmo / delete / mo_tri
# cmo / delete / mo_inlet

# cmo / copy/ mo_inlet /mo_surf
# cmo / select / mo_inlet 
# eltset/ e_delete/ itetclr ne 3
# rmpoint / element / eltset get e_delete
# rmpoint / compress
# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_inlet
# cmo / select / mo_tri 
# // dump / left.inp / mo_tri
# dump / stl / whole_block_left.stl / mo_tri
# cmo / delete / mo_tri
# cmo / delete / mo_inlet

# cmo / copy/ mo_inlet /mo_surf
# cmo / select / mo_inlet 
# eltset/ e_delete/ itetclr ne 5
# rmpoint / element / eltset get e_delete
# rmpoint / compress
# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_inlet
# cmo / select / mo_tri 
# // dump / right.inp / mo_tri
# dump / stl / whole_block_right.stl / mo_tri
# cmo / delete / mo_tri
# cmo / delete / mo_inlet

# finish 

# """
#     with open('whole_block_exterior_to_stl.lgi', 'w') as fp:
#         fp.write(lagrit_script)
#         fp.flush() 

#     lagrit_script = f"""


# read / whole_block_external.inp / mo_surf 

# cmo / select / mo_surf
# settets/normal

# ## remove surfaces for normal vector based dump
# cmo / select / mo_surf 
# eltset/ e_delete/ idface1 ne 1
# rmpoint / element / eltset get e_delete
# rmpoint / compress

# pset / p_delete / attribute / zic / 1,0,0 / lt {bottom_surface_height*0.95:0.2f}
# rmpoint  / pset get p_delete
# rmpoint / compress

# dump / tmp.inp / mo_surf 

# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_surf
# cmo / select / mo_tri 
# dump / stl / whole_block_bottom_surface.stl / mo_tri
# cmo / delete / mo_tri
# cmo / delete / mo_surf

# cmo / status / brief 



# read / whole_block_external.inp / mo_surf 

# cmo / select / mo_surf
# settets/normal

# ## remove surfaces for normal vector based dump
# cmo / select / mo_surf 
# eltset/ e_delete/ idface1 ne 2
# rmpoint / element / eltset get e_delete
# rmpoint / compress

# pset / p_delete / attribute / zic / 1,0,0 / gt {top_surface_height*0.95:0.2f}
# rmpoint  / pset get p_delete
# rmpoint / compress

# dump / tmp.inp / mo_surf 

# cmo / status / brief 
# hextotet / 2 / mo_tri  / mo_surf
# cmo / select / mo_tri 
# dump / stl / whole_block_top_surface.stl / mo_tri
# cmo / delete / mo_tri
# cmo / delete / mo_surf

# cmo / status / brief 


# finish 

# """
#     with open('whole_block_surface_to_stl.lgi', 'w') as fp:
#         fp.write(lagrit_script)
#         fp.flush() 



# def combine_stl():

#     # # Load STL files
#     # blocks = ['top','bottom']
#     # for block in blocks:
#     #     mesh_right = mesh.Mesh.from_file(f'{block}_right.stl')
#     #     mesh_left = mesh.Mesh.from_file(f'{block}_left.stl')
#     #     mesh_front = mesh.Mesh.from_file(f'{block}_front.stl')
#     #     mesh_back = mesh.Mesh.from_file(f'{block}_back.stl')
#     #     mesh_opposite = mesh.Mesh.from_file(f'{block}_opposite.stl')
#     #     mesh_surf = mesh.Mesh.from_file(f'{block}_surf.stl')

#     #     # Combine data
#     #     combined_data = np.concatenate([mesh_right.data, mesh_left.data, \
#     #                                     mesh_front.data, mesh_back.data, \
#     #                                     mesh_opposite.data, mesh_surf.data ])

#     #     # Create new mesh
#     #     combined_mesh = mesh.Mesh(combined_data)
#     #     combined_mesh.save(f'{block}_block.stl')

#     blocks = ['whole_block']
#     for block in blocks:
#         mesh_right = mesh.Mesh.from_file(f'{block}_right.stl')
#         mesh_left = mesh.Mesh.from_file(f'{block}_left.stl')
#         mesh_front = mesh.Mesh.from_file(f'{block}_front.stl')
#         mesh_back = mesh.Mesh.from_file(f'{block}_back.stl')
#         mesh_top = mesh.Mesh.from_file(f'{block}_top.stl')
#         mesh_bottom = mesh.Mesh.from_file(f'{block}_bottom.stl')
#         mesh_top_surf = mesh.Mesh.from_file(f'{block}_top_surface.stl')
#         mesh_bottom_surf= mesh.Mesh.from_file(f'{block}_bottom_surface.stl')

#         # Combine data
#         combined_data = np.concatenate([mesh_right.data, mesh_left.data, \
#                                         mesh_front.data, mesh_back.data, \
#                                         mesh_top.data, mesh_bottom.data, \
#                                          mesh_top_surf.data, mesh_bottom_surf.data ])

#         # Create new mesh
#         combined_mesh = mesh.Mesh(combined_data)
#         combined_mesh.save(f'{block}.stl')


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

def dump_stl(self, interior = True, exterior = False, cleanup_run = True):

    print("--> Dumping fracture surfaces to STL")
    top = np.reshape(self.top, self.nx*self.ny)   
    np.savetxt("top.dat", top, delimiter=" ")
    bottom = np.reshape(self.bottom, self.nx*self.ny)   
    np.savetxt("bottom.dat", bottom, delimiter=" ")

    self.write_lagrit_script_part1()
    subprocess.call('lagrit <  convert_fracture_to_avs.lgi', shell = True)

    if interior:
        write_lagrit_script_interior_stl()
        subprocess.call('lagrit <  convert_interior_to_stl.lgi', shell = True)

    if exterior:
        self.write_lagrit_script_extract_exterior()
        subprocess.call('lagrit < make_exterior_blocks.lgi', shell = True)
    if cleanup_run:
        cleanup()
    print("--> Dumping fracture surfaces to STL Complete")