# TODO: set quality
# TODO: double bond in hexgon

import os
import sys

try:
    import bpy
    from mathutils import Vector, Matrix
except ImportError:
    print('launching blender')
    os.system(r'open ~/bin/blender-2.78a-OSX_10.6-x86_64/blender.app/Contents/MacOS/blender --args --python ' + __file__)
    sys.exit()

from math import pi, sqrt, sin, cos, tan, acos

def clear_all():
    bpy.ops.object.select_by_type(type='MESH')
    bpy.ops.object.delete()


def create_atom(pos: Vector, atom_size):
    bpy.ops.mesh.primitive_uv_sphere_add(location=pos.to_tuple(), size=atom_size)
    return bpy.context.object


def create_bound(pt1: Vector, pt2: Vector, atom_size):
    dir = pt2 - pt1
    #bpy.context.object.rotation_mode = 'QUATERNION'
    length = dir.length
    origin = (pt1 + dir / 2).to_tuple()
    up = Vector((0, 0, 1))
    rotation = up.rotation_difference(dir).to_euler()
    bpy.ops.mesh.primitive_cylinder_add(radius=atom_size / 2, location=origin, depth=length, rotation=rotation)
    return bpy.context.object


def join_meshes(meshes):
    bpy.ops.object.select_all(action='DESELECT')
    for obj in meshes:
        obj.select = True
    bpy.context.scene.objects.active = meshes[0]
    bpy.ops.object.join()
    return bpy.context.object


def build_pentagon(bound_length, atom_size):
    t = bound_length
    R5 = t * sqrt(10) * sqrt(5 + sqrt(5)) / 10
    base_atoms = []
    meshes = []
    for i in range(5):
        a = 2 * pi / 5 * i
        x, y = R5 * sin(a), R5 * cos(a)
        base_atoms.append(Vector((x, y, 0)))
        meshes.append(create_atom(base_atoms[-1], atom_size))
    for i in range(5):
        meshes.append(create_bound(base_atoms[i], base_atoms[(i + 1) % 5], atom_size))
    return (join_meshes(meshes), base_atoms)


def build_hexagon(bound_length, atom_size):
    t = bound_length
    R6 = t
    base_atoms = []
    meshes = []
    for i in range(6):
        a = 2 * pi / 6 * i
        x, y = R6 * sin(a), R6 * cos(a)
        base_atoms.append(Vector((x, y, 0)))
        meshes.append(create_atom(base_atoms[-1], atom_size))
    for i in range(6):
        meshes.append(create_bound(base_atoms[i], base_atoms[(i + 1) % 6], atom_size))
    return (join_meshes(meshes), base_atoms)


def build_c60(bound_length, atom_size):
    penta, penta_points = build_pentagon(bound_length, atom_size)
    for i in range(5):
        hexa, hexa_points = build_hexagon(bound_length, atom_size)
        #r = Matrix.Rotation(pi / 2, 4, Vector((0, 1, 0)))
        hexa.select = True
        bpy.ops.transform.translate(value=penta_points[i] - hexa_points[0])
        penta_edge = (penta_points[(i + 1) % 5] - penta_points[i]).normalized()
        hexa_edge = (hexa_points[1] - hexa_points[0]).normalized()
        diff = acos(penta_edge.normalized().dot(hexa_edge.normalized()))
        bpy.ops.transform.rotate(value=-diff, axis=penta_edge.cross(hexa_edge))
        bpy.ops.transform.rotate(value=pi + pi/5, axis=penta_edge)


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
    build_c60(2, 0.2)
    render(output)
except Exception as e:
    print(e, file=sys.stderr)
#sys.exit()
