import numpy as np 
import subprocess
import shutil 
import os
import glob 

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
# ## combine meshes 
# read / surf_bottom.inp / mo_bottom 
# read / surf_top.inp / mo_top 

# addmesh / merge / mo_combined / mo_top / mo_bottom

# cmo / delete / mo_top
# cmo / delete / mo_bottom

# cmo / select / mo_combined 

# dump / fracture.inp / mo_combined
# dump / stl / fracture.stl

# finish 
"""

    with open('convert_fracture_to_avs.lgi', 'w') as fp:
        fp.write(lagrit_script)
        fp.flush() 


def write_lagrit_script_part2():
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

    with open('convert_to_stl.lgi', 'w') as fp:
        fp.write(lagrit_script)
        fp.flush() 
    
def cleanup():

    files_to_move = glob.glob("*stl")
    try:
        os.mkdir('stl_files')
    except Exception as e:
        print(e)
        print("--> Removing directory and re-making")
        shutil.rmtree('stl_files')
        os.mkdir('stl_files')
        pass 

    for f in files_to_move:
        os.rename(f, 'stl_files' + os.sep + f)

    files_to_delete = ["top.dat", "bottom.dat", "surf_top.inp" ,"surf_bottom.inp", "convert_fracture_to_avs.lgi", "convert_to_stl.lgi"]

    for f in files_to_delete:
        try:
            os.remove(f)
        except Exception as e:
            print(e)
            pass 

def dump_stl(self, cleanup_run = True):

    print("--> Dumping fracture surfaces to STL")
    top = np.reshape(self.top, self.nx*self.ny)   
    np.savetxt("top.dat", top, delimiter=" ")
    bottom = np.reshape(self.bottom, self.nx*self.ny)   
    np.savetxt("bottom.dat", bottom, delimiter=" ")

    self.write_lagrit_script_part1()
    subprocess.call('lagrit <  convert_fracture_to_avs.lgi', shell = True)

    write_lagrit_script_part2()
    subprocess.call('lagrit <  convert_to_stl.lgi', shell = True)

    if cleanup_run:
        cleanup()
    print("--> Dumping fracture surfaces to STL Complete")