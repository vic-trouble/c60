# TODO: set quality
# TODO: double bond in hexgon

import os
import sys

try:
    import bpy
    from mathutils import Vector
except ImportError:
    print('launching blender')
    os.system(r'"C:\Program Files\Blender Foundation\Blender\blender.exe" --python ' + __file__)
    sys.exit()

from math import pi, sqrt, sin, cos, tan

def clear_all():
    bpy.ops.object.select_by_type(type='MESH')
    bpy.ops.object.delete()

def create_atom(pos):
    bpy.ops.mesh.primitive_uv_sphere_add(location=pos)

def create_bound(pos1, pos2):
    dir = Vector(pos2[i] - pos1[i] for i in range(3))
    #bpy.context.object.rotation_mode = 'QUATERNION'
    length = dir.length
    origin = (Vector(pos1) + dir / 2).to_tuple()
    up = Vector((0, 0, 1))
    rotation = up.rotation_difference(dir).to_euler()
    bpy.ops.mesh.primitive_cylinder_add(radius=0.5, location=origin, depth=length, rotation=rotation)
    
    
def build_c60(bound_length):
    t = bound_length
    R5 = t * sqrt(10) * sqrt(5 + sqrt(5)) / 10
    
    # build base pentagon
    base_atoms = []
    for i in range(5):
        a = 2 * pi / 5 * i
        x, y = R5 * sin(a), R5 * cos(a)
        base_atoms.append((x, y, 0))
        create_atom(base_atoms[-1])
        
    for i in range(5):
        create_bound(base_atoms[i], base_atoms[(i + 1) % 5])
    
    
def render(filename):
    bpy.context.scene.render.filepath = filename
    bpy.context.scene.render.resolution_x = 1024
    bpy.context.scene.render.resolution_y = 768
    bpy.ops.render.render(write_still=True)

##########################################################################################

output = os.path.join(os.getcwd(), '1.png')
if os.path.exists(output):
    os.unlink(output)

try:        
    clear_all()
    build_c60(10)
    render(output)
except Exception as e:
    print(e, file=sys.stderr)
sys.exit()
