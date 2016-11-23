# TODO: set quality
# TODO: double bond in hexgon
# TODO: rename bound -> bond

import os
import sys
from _collections import defaultdict

try:
    import bpy
    from mathutils import Vector, Matrix
except ImportError:
    print('launching blender')
	# os.system(r'open ~/bin/blender-2.78a-OSX_10.6-x86_64/blender.app/Contents/MacOS/blender --args --python ' + __file__)
    os.system(r'"C:\Program Files\Blender Foundation\Blender\blender.exe" --python ' + __file__)
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


def align(shape, src_vec, dst_vec, rev=False):
    src_vec = src_vec.normalized()
    dst_vec = dst_vec.normalized()
    rotate(shape, acos(src_vec.dot(dst_vec)) * (-1 if rev else 1), src_vec.cross(dst_vec))


def flesh_out(shapes, bound_length, atom_size):
    eps = 0.05
    atoms = []
    bonds = []
    graph = defaultdict(list) # atom index -> [atom index]
    belong = defaultdict(list) # atom index -> [shape index]
    meshes = []
    for j, shape in enumerate(shapes):
        # fill-in helper collections
        indices = []
        for atom in shape:
            min_distance = None
            index = None
            for i, a in enumerate(atoms):
                distance = (a - atom).length
                if min_distance is None or distance < min_distance:
                    min_distance = distance
                    index = i
            if index is None or min_distance > eps:
                index = len(atoms)
                atoms.append(atom)
            indices.append(index)

        for i in indices:
            belong[i].append(j)

        # create primitive meshes
        for i in range(len(indices)):
            meshes.append(create_atom(atoms[indices[i]], atom_size))
            
            bond_i = indices[i]
            bond_j = indices[i - 1]
            graph[bond_i].append(bond_j)

            start = atoms[bond_i]
            end = atoms[bond_j]
            distance = 1000
            for b in bonds:
                distance = min(distance, (b[0] - start).length + (b[1] - end).length,
                                         (b[1] - start).length + (b[0] - end).length)
            if distance > eps:
                meshes.append(create_bound(start, end, atom_size))
                bonds.append((start, end))

    assert len(graph) == 60
    
    # annotate double-bonds
    hexa = [k for k in range(len(shapes)) if len(shapes[k]) == 6]
    for i in graph:
        for j in graph[i]:
            ij_shapes = list(set(belong[i]) & set(belong[j]))
            ij_hexa = [k in hexa for k in ij_shapes]
            print(ij_shapes, ij_hexa)
            if sum(ij_hexa) == 2:
                shape_index = ij_hexa.index(True)
                shape = shapes[ij_shapes[shape_index]]

                center = sum(shape, Vector((0, 0, 0))) / len(shape)

                start = atoms[i]
                end = atoms[j]
                a = 0.75
                b = 1 - a
                sub_start = start * a + end * b
                sub_end = start * b + end * a
                c = 0.85
                d = 1 - c
                sub_start = sub_start * c + center * d
                sub_end = sub_end * c + center * d
                meshes.append(create_bound(sub_start, sub_end, atom_size * 2 / 3))

    print('{} atoms, {} bonds'.format(len(atoms), len(bonds)))
    return join_meshes(meshes)


def build_c60(bound_length, atom_size):
    # base penta
    penta = build_pentagon(bound_length)

    # first hexa belt
    hexa1 = []
    unfold = 0.65
    for i in range(5):
        hexa = build_hexagon(bound_length)
        penta_edge = penta[(i + 1) % 5] - penta[i]
        hexa_edge = hexa[1] - hexa[0]
        align(hexa, penta_edge, hexa_edge, rev=True)
        rotate(hexa, pi + unfold, penta_edge)
        translate(hexa, penta[i] - hexa[0])
        hexa1.append(hexa)

    # interleaving penta
    ipenta1 = []
    for i in range(5):
        ipenta = build_pentagon(bound_length)
        ipenta_normal = (ipenta[1] - ipenta[0]).cross(ipenta[-1] - ipenta[0])

        hexa = hexa1[i]
        nhexa = hexa1[(i + 1) % 5]
        hexa_normal = (hexa[3] - hexa[2]).cross(nhexa[4] - nhexa[5])

        align(ipenta, ipenta_normal, hexa_normal)

        ipenta_edge = ipenta[1] - ipenta[0]
        hexa_edge = hexa[3] - hexa[2]
        align(ipenta, ipenta_edge, hexa_edge)

        translate(ipenta, hexa[2] - ipenta[0])
        ipenta1.append(ipenta)

    # interleaving hexa
    ihexa1 = []
    for i in range(5):
        hexa = build_hexagon(bound_length)
        hexa_normal = (hexa[1] - hexa[0]).cross(hexa[-1] - hexa[0])

        ipenta = ipenta1[i]
        nipenta = ipenta1[(i - 1 + 5) % 5]
        ipenta_normal = (ipenta[2] - ipenta[1]).cross(nipenta[4] - nipenta[3])
        align(hexa, hexa_normal, ipenta_normal)

        ipenta_edge = ipenta[1] - ipenta[2]
        hexa_edge = hexa[1] - hexa[0]
        align(hexa, hexa_edge, ipenta_edge)

        translate(hexa, ipenta[2] - hexa[0])
        ihexa1.append(hexa)

    # another interleaving hexa
    ihexa2 = []
    for i in range(5):
        ihexa = build_hexagon(bound_length)
        ihexa_normal = (ihexa[1] - ihexa[0]).cross(ihexa[-1] - ihexa[0])

        hexa = ihexa1[i]
        nhexa = ihexa1[(i + 1) % 5]
        hexa_normal = (hexa[5] - hexa[0]).cross(nhexa[4] - nhexa[3])
        align(ihexa, ihexa_normal, hexa_normal)

        hexa_edge = hexa[5] - hexa[0]
        ihexa_edge = ihexa[1] - ihexa[0]
        align(ihexa, ihexa_edge, hexa_edge)

        translate(ihexa, hexa[0] - ihexa[0])
        ihexa2.append(ihexa)

    # another interleaving penta
    ipenta2 = []
    for i in range(5):
        ipenta = build_pentagon(bound_length)
        ipenta_normal = (ipenta[1] - ipenta[0]).cross(ipenta[-1] - ipenta[0])

        hexa = ihexa2[i]
        phexa = ihexa2[i - 1]
        hexa_normal = (hexa[2] - hexa[1]).cross(phexa[4] - phexa[3])
        align(ipenta, ipenta_normal, hexa_normal)

        hexa_edge = hexa[2] - hexa[1]
        ipenta_edge = ipenta[0] - ipenta[1]
        align(ipenta, ipenta_edge, hexa_edge)

        translate(ipenta, hexa[2] - ipenta[0])
        ipenta2.append(ipenta)

    # second hexa belt
    hexa2 = []
    for i in range(5):
        hexa = build_hexagon(bound_length)
        hexa_normal = (hexa[1] - hexa[0]).cross(hexa[-1] - hexa[0])

        ipenta = ipenta2[i]
        nipenta = ipenta2[(i + 1) % 5]
        ipenta_normal = (ipenta[4] - ipenta[0]).cross(nipenta[4] - nipenta[3])
        align(hexa, hexa_normal, ipenta_normal)

        ipenta_edge = ipenta[4] - ipenta[0]
        hexa_edge = hexa[1] - hexa[0]
        align(hexa, hexa_edge, ipenta_edge)

        translate(hexa, ipenta[0] - hexa[0])
        hexa2.append(hexa)

    return flesh_out(hexa1 + ipenta1 + ihexa1 + ihexa2 + ipenta2 + hexa2, bound_length, atom_size)


def render(filename):
    bpy.context.scene.render.filepath = filename
    bpy.context.scene.render.resolution_x = 1600
    bpy.context.scene.render.resolution_y = 1200
    bpy.ops.render.render(write_still=True)

##########################################################################################

for i in (2,):
    output = os.path.join(os.getcwd(), '{}.png'.format(i))
    if os.path.exists(output):
        os.unlink(output)

    clear_all()
    build_c60(0.8, 0.08 * i)
    render(output)
#sys.exit()
