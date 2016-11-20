# TODO: set quality
# TODO: double bond in hexgon
# TODO: rename bound -> bond

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


def build_pentagon(bound_length):
    t = bound_length
    R5 = t * sqrt(10) * sqrt(5 + sqrt(5)) / 10
    penta = []
    for i in range(5):
        a = 2 * pi / 5 * i
        x, y = R5 * sin(a), R5 * cos(a)
        penta.append(Vector((x, y, 0)))
    return penta


def build_hexagon(bound_length):
    t = bound_length
    R6 = t
    hexa = []
    for i in range(6):
        a = 2 * pi / 6 * i
        x, y = R6 * sin(a), R6 * cos(a)
        hexa.append(Vector((x, y, 0)))
    return hexa

def translate(shape, vec):
    m = Matrix.Translation(vec)
    for i in range(len(shape)):
        shape[i] = m * shape[i]


def rotate(shape, value, axis):
    m = Matrix.Rotation(value, 4, axis)
    for i in range(len(shape)):
        shape[i] = m * shape[i]


def flesh_out(shapes, bound_length, atom_size):
    for shape in shapes:
        for i in range(len(shape)):
            create_atom(shape[i], atom_size)
            create_bound(shape[i], shape[i - 1], atom_size)

def build_c60(bound_length, atom_size):
    # base penta
    penta = build_pentagon(bound_length)

    # first hexa belt
    hexa1 = []
    for i in range(5):
        hexa = build_hexagon(bound_length)
        penta_edge = (penta[(i + 1) % 5] - penta[i]).normalized()
        hexa_edge = (hexa[1] - hexa[0]).normalized()
        diff = acos(penta_edge.normalized().dot(hexa_edge.normalized()))
        rotate(hexa, -diff, penta_edge.cross(hexa_edge))
        rotate(hexa, pi + pi / 5, penta_edge) # TODO: pi/5 is heuristic
        translate(hexa, penta[i] - hexa[0])
        hexa1.append(hexa)

    flesh_out([penta] + hexa1, bound_length, atom_size)
    # interleaving penta
    #for i in range(1):
    #    inter_penta, ipenta_pts = build_pentagon(bound_length, atom_size)


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
