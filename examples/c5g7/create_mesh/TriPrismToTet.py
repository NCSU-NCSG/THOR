"""
  Script for reading a mesh made of triangular prisms, break it into
  tetrahedra, and write it in THOR format
"""
import sys
import os
import os.path
import math
import time

"""
  Determines if two floats are equal within tolerance
:x: first float
:y: second float
:tol: (optional) tolerance
:returns: True if x,y are within tol, False otherwise
"""
def fuzzy_equal(x, y, tol = 1.0e-8):
    if abs(x-y) <= tol:
        return True
    else:
        return False

"""
    Determines which boundary a side is on
::
"""
def which_boundary_am_i(tet_xyz, side, mesh_bbox):
    if side == 0:
        node1 = tet_xyz[1]
        node2 = tet_xyz[2]
        node3 = tet_xyz[3]
    elif side == 1:
        node1 = tet_xyz[0]
        node2 = tet_xyz[2]
        node3 = tet_xyz[3]
    elif side == 2:
        node1 = tet_xyz[0]
        node2 = tet_xyz[1]
        node3 = tet_xyz[3]
    elif side == 3:
        node1 = tet_xyz[0]
        node2 = tet_xyz[1]
        node3 = tet_xyz[2]
    else:
        sys.exit("Incorrect side id.")
    if fuzzy_equal(node1[0], mesh_bbox[0]) and fuzzy_equal(node2[0], mesh_bbox[0]) and\
       fuzzy_equal(node3[0], mesh_bbox[0]):
       return "-x"
    elif fuzzy_equal(node1[0], mesh_bbox[1]) and fuzzy_equal(node2[0], mesh_bbox[1]) and\
       fuzzy_equal(node3[0], mesh_bbox[1]):
       return "+x"
    elif fuzzy_equal(node1[1], mesh_bbox[2]) and fuzzy_equal(node2[1], mesh_bbox[2]) and\
       fuzzy_equal(node3[1], mesh_bbox[2]):
       return "-y"
    elif fuzzy_equal(node1[1], mesh_bbox[3]) and fuzzy_equal(node2[1], mesh_bbox[3]) and\
       fuzzy_equal(node3[1], mesh_bbox[3]):
       return "+y"
    elif fuzzy_equal(node1[2], mesh_bbox[4]) and fuzzy_equal(node2[2], mesh_bbox[4]) and\
       fuzzy_equal(node3[2], mesh_bbox[4]):
       return "-z"
    elif fuzzy_equal(node1[2], mesh_bbox[5]) and fuzzy_equal(node2[2], mesh_bbox[5]) and\
       fuzzy_equal(node3[2], mesh_bbox[5]):
       return "+z"
    else:
        print node1
        print node2
        print node3
        sys.exit("Could not classify this boundary side.")

"""
  Computes the bounding box around the mesh
:nodes: nodes of the mesh
"""
def mesh_bounding_box(nodes):
    xmin = sys.float_info.max
    xmax = -sys.float_info.max
    ymin = sys.float_info.max
    ymax = -sys.float_info.max
    zmin = sys.float_info.max
    zmax = -sys.float_info.max
    for node_id in nodes.keys():
        node = nodes[node_id]
        if node[0] < xmin: xmin = node[0]
        if node[0] > xmax: xmax = node[0]
        if node[1] < ymin: ymin = node[1]
        if node[1] > ymax: ymax = node[1]
        if node[2] < zmin: zmin = node[2]
        if node[2] > zmax: zmax = node[2]
    return [xmin, xmax, ymin, ymax, zmin, zmax]

"""
  Computes the adjacency list
:tet_list: array of all tets
:node_list: array of all nodes
:tet_2_node: map from tet to nodes
:node_2_tet: map listing the tets a node is in
"""
def get_adjaceny_list(tet_list, node_list ,tet_2_node,node_2_tet):
  adjacency_list = {}
  for tet in tet_list:
      adjacency_list[tet] = [-1,-1,-1,-1]
  common_nodes = {}
  for tet in tet_list:
      common_nodes[tet] = []
  for node in node_list:
      tets = node_2_tet[node]
      for tet1 in tets:
          for tet2 in tets:
              if tet1 != tet2:
                  common_nodes[tet1].append(tet2)
  for tet1 in common_nodes.keys():
      counter = {}
      for tet2 in common_nodes[tet1]:
          counter[tet2] = 0
      for tet2 in common_nodes[tet1]:
          counter[tet2] += 1
          if counter[tet2] == 3 and tet1 > tet2:
              nodes1 = tet_2_node[tet1]
              nodes2 = tet_2_node[tet2]
              for j in range(4):
                  node = nodes1[j]
                  if node not in nodes2: side1 = j
              for j in range(4):
                  node = nodes2[j]
                  if node not in nodes1: side2 = j
              adjacency_list[tet1][side1] = tet2
              adjacency_list[tet2][side2] = tet1
  """
  key = str(tet1) + str(tet2)
  #print key
  if key in common_node_counter.keys():
      common_node_counter[key] += 1
  else:
      common_node_counter[key] = 1
  if common_node_counter[key] == 3:
      nodes1 = tet_2_node[tet1]
      nodes2 = tet_2_node[tet2]
      for j in range(4):
          node = nodes1[j]
          if node not in nodes2: side1 = j
      for j in range(4):
          node = nodes2[j]
          if node not in nodes1: side2 = j
      adjacency_list[tet1][side1] = tet2
      adjacency_list[tet2][side2] = tet1
  """

  return adjacency_list

"""
  Reads a triangular prism mesh
:filename: filename of the mesh file
:format: selection of different formats
:returns:
"""
def readTriPrismMesh(filename, mesh_format = "simple"):
    if not os.path.isfile(filename):
        sys.exit("File "+filename+" does not exist.")
    lines = [line.rstrip('\n') for line in open(filename)]
    n_nodes = 0
    n_elem  = 0
    nodes = {}
    subdomains = {}
    src_regions = {}
    prisms = {}
    if mesh_format == "simple":
        try:
            # header
            n_nodes = int(lines[0])
            n_elem  = int(lines[1])

            # nodes
            for j in range(n_nodes):
                sline = lines[4+j].split()
                node_id = int(sline[0])
                x = float(sline[1])
                y = float(sline[2])
                z = float(sline[3])
                nodes[node_id] = [x,y,z]

            # regions
            for j in range(n_elem):
                sline = lines[4+n_nodes+j].split()
                elem_id = int(sline[0])
                sub_id = int(sline[1])
                src_id = int(sline[2])
                subdomains[elem_id] = sub_id
                src_regions[elem_id] = src_id

            # prisms
            for j in range(n_elem):
                sline = lines[4+n_elem+n_nodes+j].split()
                elem_id = int(sline[0])
                prisms[elem_id] = [int(sline[1]),int(sline[2]),int(sline[3]),\
                                   int(sline[4]),int(sline[5]),int(sline[6])]
        except:
            sys.exit("Error when reading mesh file "+filename+" with mesh format "+mesh_format)
    else:
        sys.exit("Mesh format "+mesh_format+" is not implemented.")

    return nodes, subdomains, src_regions, prisms

""" Converts a prismatic mesh into a tetmesh
:prism_elems: element to node dict of prismatic mesh
:prism_subdomains: subdomain dict for prismatic mesh
:prism_src_regions: source region dict for prismatic mesh
:returns: 3 dicts corresponding to the three inputs for tets
"""
def tetMeshFromTriPrism(prism_elems, prism_subdomains, prism_src_regions):

    # Map prism to tet ids
    tet_id_to_prism_id = {}
    prism_id_to_tet_id = {}
    tet = 1
    for prism in sorted(prism_elems.keys()):
        tet_id_to_prism_id[tet] = prism
        prism_id_to_tet_id[prism] = [tet]
        tet += 1
        tet_id_to_prism_id[tet] = prism
        prism_id_to_tet_id[prism].append(tet)
        tet += 1
        tet_id_to_prism_id[tet] = prism
        prism_id_to_tet_id[prism].append(tet)
        tet += 1

    # Assign the region and src ids
    tet_subdomains = {}
    tet_src_regions = {}
    for tet in sorted(tet_id_to_prism_id.keys()):
        prism = tet_id_to_prism_id[tet]
        tet_subdomains[tet] = prism_subdomains[prism]
        tet_src_regions[tet] = prism_src_regions[prism]

    # Assign the tet nodes
    tet_elems = {}
    for prism in sorted(prism_id_to_tet_id.keys()):
      prism_nodes = prism_elems[prism]
      tet1 = prism_id_to_tet_id[prism][0]
      tet2 = prism_id_to_tet_id[prism][1]
      tet3 = prism_id_to_tet_id[prism][2]
      smallest_pos = prism_nodes.index(min(prism_nodes))
      table = []
      if smallest_pos == 0:
        table = [1,2,3,4,5,6]
      elif smallest_pos == 1:
        table = [2,3,1,5,6,4]
      elif smallest_pos == 2:
        table = [3,1,2,6,4,5]
      elif smallest_pos == 3:
        table = [4,6,5,1,3,2]
      elif smallest_pos == 4:
        table = [5,4,6,2,1,3]
      elif smallest_pos == 5:
        table = [6,5,4,3,2,1]
      else:
        sys.exit("Invalid minimal index when finding Triangluar prism to tet map")

      i1 = prism_nodes[table[0]-1]
      i2 = prism_nodes[table[1]-1]
      i3 = prism_nodes[table[2]-1]
      i4 = prism_nodes[table[3]-1]
      i5 = prism_nodes[table[4]-1]
      i6 = prism_nodes[table[5]-1]

      if (min(i2,i6) < min(i3,i5)):
        tet_elems[tet1]  = ([i1,i2,i3,i6])
        tet_elems[tet2] = ([i1,i2,i6,i5])
        tet_elems[tet3] = ([i1,i5,i6,i4])
      else:
        tet_elems[tet1]   = ([i1,i2,i3,i5])
        tet_elems[tet2] = ([i1,i5,i3,i6])
        tet_elems[tet3] = ([i1,i5,i6,i4])

    return tet_elems, tet_subdomains, tet_src_regions

'''
  Writes a THOR Mesh
:nodes: dictionary of nodes
:elems: dictionary of elements (elem_id -> nodes)
:subdomains: dictionary of elem->subdomains region
:src_regions: dictionary of elem->src region
:returns: None
'''
def writeTHORMesh(nodes, elems, subdomains, src_regions, boundary_conditions, filename):
    fid = open(filename,"w")
    n_nodes = len(nodes.keys())
    n_elems = len(elems.keys())
    fid.write("%d \n" % (n_nodes))
    fid.write("%d \n 1\n 1\n" % (n_elems))

    # write the nodes
    for node in nodes.keys():
        string = '%10d %25.16e %25.16e %25.16e\n' % (node , nodes[node][0], nodes[node][1], nodes[node][2])
        fid.write(string)

    # write the tets
    for tet in elems.keys():
      fid.write("%d %d %d\n" % (tet, subdomains[tet], src_regions[tet]))

    for tet in elems.keys():
      fid.write("%d %d %d %d %d\n" % (tet,elems[tet][0],elems[tet][1],elems[tet][2],elems[tet][3]))

    # get the adjacency list
    node_2_tet = {}
    for node in nodes:
        node_2_tet[node] = []
    for tet in elems.keys():
        for j in range(4):
            n = elems[tet][j]
            node_2_tet[n].append(tet)

    start_adj_time = time.time()
    print "Entered adjacency list routine..."
    adjlist = get_adjaceny_list(elems.keys(), node_2_tet.keys(), elems, node_2_tet)
    print "Finished adjacency list ",time.time() - start_adj_time

    # write sides
    nb = 0
    for tet in adjlist.keys():
        for side in range(4):
            if adjlist[tet][side] == -1:
                nb += 1
    fid.write("%d\n" % (nb))

    mesh_bbox = mesh_bounding_box(nodes)
    print "Mesh BBOX"
    print mesh_bbox
    for tet in adjlist.keys():
        for side in range(4):
            if adjlist[tet][side] == -1:
                tet_xyz = [[] for j in range(4)]
                tet_xyz[0] = nodes[elems[tet][0]]
                tet_xyz[1] = nodes[elems[tet][1]]
                tet_xyz[2] = nodes[elems[tet][2]]
                tet_xyz[3] = nodes[elems[tet][3]]
                boundary_type = which_boundary_am_i(tet_xyz, side, mesh_bbox)
                fid.write("%d %d %d \n" % (tet,side,boundary_conditions[boundary_type]))

    # write adjacency list
    fid.write("%d\n" % (4 * n_elems))
    for tet in elems.keys():
      for side in range(4):
        if adjlist[tet][side] == -1:
          fid.write("%d %d %d %d\n" % (tet,side,0,0))
        else:
          neigh = adjlist[tet][side]
          for j in range(4):
            if adjlist[neigh][j] == tet: neigh_side = j
          fid.write("%d %d %d %d\n" % (tet,side,neigh,neigh_side))
    fid.close()

"""
 Some data
"""
# xmin, xmax, ymin, ymax, zmin, zmax
boundary_conditions = {"-x" : 0, "+x" : 1, "-y" : 0,"+y" : 1, "-z" : 0,"+z" : 1}
modify_z_coordinates = {57.12 : 64.26}

"""
 Start actual script
"""
start_time = time.time()
nodes, prism_subdomains, prism_src_regions, prism_elems = readTriPrismMesh("c5g7_3d_hex.dat", mesh_format = "simple")

for node_id in nodes.keys():
    node = nodes[node_id]
    for zold in modify_z_coordinates.keys():
        if fuzzy_equal(zold, node[2]):
            nodes[node_id][2] = modify_z_coordinates[zold]

print "Reading file took ", time.time() - start_time, " seconds"
tet_elems, tet_subdomains, tet_src_regions = tetMeshFromTriPrism(prism_elems, prism_subdomains, prism_src_regions)
all_subdomains = [0 for j in range(1000)]
for tet in tet_subdomains.keys():
    all_subdomains[tet_subdomains[tet]] = 1
print "-- Subdomains present in mesh --"
for j in range(1000):
    if all_subdomains[j] ==1:
        print j
print "Processing took ", time.time() - start_time, " seconds"
writeTHORMesh(nodes, tet_elems, tet_subdomains, tet_src_regions, boundary_conditions, "c5g7_3d_tet.dat")
print "Writing took ", time.time() - start_time, " seconds"
