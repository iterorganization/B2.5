#! /usr/bin/env python
# This utility converts edge CPO to IMAS edge_profiles IDS
#
# Quick start:
# cd /pfs/work/kosl/solps-iter-ids/modules/B2.5/src/ids
# Definition of the class structures is in file imas.py!
# To run this converter use:
# module use -a ~dkaljun/imas/etc/modulefiles
# module avail imas
# module load imas/develop/3/ual/develop
# module display imas/develop/3/ual/develop
# python cpo2ids.py
# python cpo2ids.py --shot=16151 --run=1000 --user=kosl --tokamak=aug --version=4.10a
# Resulting files are stored under $HOME/public/imasdb/solps-iter/3/0
# as ids_161511000.*
# Further documentation and notes:
# dd_doc
# imasdb -l # should show solps-iter with asterisk
# ls $UAL/pythonExamples
# cd ~/itmggd
# ls /pfs/home/kosl/itmggd/branches/4.10a/python/itmggd/examples
# less $UAL/python_pk/python2.6/ual/edge.py
# /pfs/work/kosl/bin/pycharm.sh &
#NODES: edge::grid::spaces::objects(1) -> ggd(:)/grid/space(:)/objects_per_dimension(1)/object(:)
  #x: edge/grid/spaces/objects(0)/nodes/geo(i,0) 	)
  #y: edge/grid/spaces/objects(0)/nodes/geo(i,1)	)-> ggd(:)/grid/space(:)/objects_per_dimension(0)/object(0)/geometry(:)
  #z: 0.0					 	)
import ual 
import imas
import numpy
import sys
import getopt

'''
This sample program will open an existing pulse file (shot 16151, run 1000) and will
read the stored edge_profiles IDSs.

It will then output the content of some fields of the edge_profiles IDSs.
'''

# This routine reads an edge_profiles IDS from the database, filling
# some fields of the IDS


ual.ualdef.DEVEL_DEBUG = 0
ual.ualdef.PRINT_DEBUG = 0
ual.ualdef.VERBOSE_DEBUG = 0


def read_ids():
    """Class IMAS is the main class for the UAL.

    It contains a set of field classes, each corresponding to an IDS
    defined in the UAL The parameters passed to this creator define the
    shot and run number. The second pair of arguments defines the
    reference shot and run and is used when the a new database is
    created, as in this example.

    """
    print 'Reading IDS: '
    my_ids_obj = imas.ids(shot, run, shot, run)


    my_ids_obj.open()  # Open the database
    #my_ids_obj.open_env(user, tokamak, version)

    if my_ids_obj.isConnected():
        print 'open OK!'
    else:
        print 'open FAILED!'
        sys.exit()

    my_ids_obj.edge_profiles.get()

    # Example for reading data from IDS and CPO
    print 'IDS'
    # print my_ids_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0]
    # print 'CPO'
    # print edge.fluid

    # Write output to text:
    # f = open("cpo2ids_output.txt", "w")
    # f.write(str(edge.species))
    # f.close()

    my_ids_obj.close()

def write_ids():
    print 'Writing IDS: '
    time = 1
    interp = 1
    #treename = "ids"

    imas_obj = imas.ids(shot, run, shot, run)

    imas_obj.create() #Create the data entry
    #imas_obj.create_env(user, tokamak, version)

    if imas_obj.isConnected():
        print 'Creation of data entry OK!'
    else:
        print 'Creation of data entry FAILED!'
        sys.exit()

    imas_obj.edge_profiles.profiles_1d.resize(1)
    imas_obj.edge_profiles.ggd.resize(1)
    imas_obj.edge_profiles.putNonTimed()

    # ----Allocating space for Nodes, faces, faces---- #
    imas_obj.edge_profiles.ggd[0].grid.space.resize(1)
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension.resize(3)

    imas_obj.edge_profiles.ids_properties.homogeneous_time = 1 # Mandatory to define this property (?)

    num_nodes_all = len(edge.grid.spaces[0].objects[0].geo)
    imas_obj.edge_profiles.time.resize(num_nodes_all * 2)   #(?) Not sure how big should it be
    ##print "num_nodes_all: ", num_nodes_all


    imas_obj.edge_profiles.ggd[0].grid.space[0].coordinates_type.resize(1)

    num_subgrids = len(edge.grid.subgrids)
    subgrid_name = numpy.array([])

    ##print "num_subgrids: ", num_subgrids

    # First we must know how many subgrids of each class we have. We have to know that in order to properly resize our .object structure.
    # Structure can be resized when needed, but then all data already inside gets deleted. Also putSlice() command can be used only once, to use it again then we
    # would probably have to open the database again (test runs have shown so)

    #Allocating .object space for each type of class

    num_obj_class_0 = 0
    num_obj_class_1 = 0
    num_obj_class_2 = 0

    for i in range(num_subgrids):
        subgrid_class = edge.grid.subgrids[i].list[0].cls[0]
        if subgrid_class == 0:
            num_obj_class_0 += 1
        elif subgrid_class == 1:
            num_obj_class_1 += 1
        else:
            num_obj_class_2 += 1

    ##print "Number of objects class 0: ", num_obj_class_0
    ##print "Number of objects class 1: ", num_obj_class_1
    ##print "Number of objects class 2: ", num_obj_class_2

    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object.resize(num_obj_class_0) #nodes
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[1].object.resize(num_obj_class_1) #edges
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[2].object.resize(num_obj_class_2) #faces/cells
    imas_obj.edge_profiles.ggd[0].grid.grid_subset.resize(num_subgrids)

    obj_class_0_id = -1
    obj_class_1_id = -1
    obj_class_2_id = -1
    num_cells_all_nodes = 0

    for i in range(num_subgrids):
        subgrid_id = i
        subgrid_name = numpy.append(subgrid_name, edge.grid.subgrids[i].id)
        subgrid_class = edge.grid.subgrids[i].list[0].cls[0]
        ##print "Subgrid base id: ", i+1

        imas_obj.edge_profiles.ggd[0].grid.grid_subset[i].identifier.name = subgrid_name[i]
        imas_obj.edge_profiles.ggd[0].grid.grid_subset[i].identifier.index = subgrid_id + 1
        imas_obj.edge_profiles.ggd[0].grid.grid_subset[i].element.resize(1)
        imas_obj.edge_profiles.ggd[0].grid.grid_subset[i].element[0].object.resize(1)
        imas_obj.edge_profiles.ggd[0].grid.grid_subset[i].element[0].object[0].space = 0 + 1
        imas_obj.edge_profiles.ggd[0].grid.grid_subset[i].element[0].object[0].dimension = subgrid_class + 1

        indset_found = len(edge.grid.subgrids[i].list[0].indset)

        range0 = 0;
        range1 = 0;
        jIndex = 0;
        ind_num = len(edge.grid.subgrids[i].list[0].ind)
        if indset_found > 0:
            range_size = len(edge.grid.subgrids[i].list[0].indset[0].range)
            if range_size > 1:
                range_found = 1
                range0 = edge.grid.subgrids[i].list[0].indset[0].range[0]-1
                range1 = edge.grid.subgrids[i].list[0].indset[0].range[1]
                start_index = range0
                end_index = range1
                cpo_array_size = end_index-start_index
            else:
                print "Range is either EMPTY or it doesn't exist"
        else :
            range_found = 0
            start_index  = edge.grid.subgrids[i].list[0].ind[0,0]-1
            end_index = edge.grid.subgrids[i].list[0].ind[ind_num-1,0]-1
            # print "Range not found"
            cpo_array_size = ind_num

        ##print 'subgrid', '{0:3}'.format(i+1), 'name: ', '{0:18}'.format(subgrid_name[i]), " class: ", '{0:1}'.format(subgrid_class), 'cpo_array_size: ', cpo_array_size

        index_array = numpy.array([], dtype='int')

        if range_found == 1:
            for j in range(start_index, end_index):
                index_array = numpy.insert(index_array, j-start_index, j)
        else:
            for j in range(0, ind_num):
                jIndex = edge.grid.subgrids[i].list[0].ind[j,0] - 1 #because in Fortran indiced go 1,2,3... while in C go 0, 1, 2
                index_array = numpy.insert(index_array, j, jIndex)

        size = len(index_array)

        # ----Allocating space and writing to IDS: Geometry and Nodes for class 0 -> Nodes---- #
        if subgrid_class == 0:
            obj_class_0_id += 1
            imas_obj.edge_profiles.ggd[0].grid.grid_subset[i].element[0].object[0].index = obj_class_0_id+1

            if subgrid_name[i] == "Nodes": #because we have geometry only in nodes. Geometry of other subgrids is taken using subgrid indices in connection with Nodes geometry (for now)
                imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_0_id].geometry.resize(size * 2)
                # it should create (x,y) -> x*y elements = size
                imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_0_id].nodes.resize(size)

                for j in range(size):
                    # There are size*2 = 7176 geometry entries. [x1, x2, ... xn, y1, y2, ...yn] -> Fortran notation (?)
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[obj_class_0_id].geometry[j] = \
                        edge.grid.spaces[0].objects[0].geo[j, 0]
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[obj_class_0_id].geometry[size + j] = \
                        edge.grid.spaces[0].objects[0].geo[j, 1]
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[obj_class_0_id].nodes[j] = index_array[j]+1 #to put in Fortran style: indices start with 1
                    imas_obj.edge_profiles.time[j] = time;
            else:
                for j in range(size):
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_0_id].geometry.resize(size * 2)
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_0_id].nodes.resize(size)

                    k = index_array[j]
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[obj_class_0_id].geometry[j] = \
                        edge.grid.spaces[0].objects[0].geo[k, 0]
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[obj_class_0_id].geometry[size + j] = \
                        edge.grid.spaces[0].objects[0].geo[k, 1]
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[obj_class_0_id].nodes[j] = index_array[j] + 1  # to put in Fortran style: indices start with 1
                    imas_obj.edge_profiles.time[j] = time;
                    # print imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[obj_class_0_id].nodes[j]

        # ----Allocating space and writing to IDS: Geometry and Nodes for class 1 -> Edges---- #
        elif subgrid_class == 1:
            obj_class_1_id += 1
            imas_obj.edge_profiles.ggd[0].grid.grid_subset[i].element[0].object[0].index = obj_class_1_id +1

            imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_1_id].nodes.resize(size * 2)  # each edge consists of 2 nodes->7044*2
            for j in range(size):
                # imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_1_id].geometry.resize(size * 4)

                k = index_array[j]
                # imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[obj_class_1_id].geometry[j] = \
                #     edge.grid.spaces[0].objects[0].geo[k, 0]
                # imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[obj_class_1_id].geometry[size + j] = \
                #     edge.grid.spaces[0].objects[0].geo[k, 1]

                imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_1_id].nodes[j] = \
                    edge.grid.spaces[0].objects[1].boundary[k,0]

                imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_1_id].nodes[size+j] = \
                    edge.grid.spaces[0].objects[1].boundary[k,1]

                #imas_obj.edge_profiles.time[j] = time;

        # ----Allocating space and writing to IDS: Geometry and Nodes for class 2 -> Cells---- #
        elif subgrid_class == 2:
            obj_class_2_id += 1
            imas_obj.edge_profiles.ggd[0].grid.grid_subset[i].element[0].object[0].index = obj_class_2_id +1

            imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_2_id].nodes.resize(size*4) #each cell consisty of 4 nodes

            if subgrid_name[i] == "Cells":
                num_cells_all_nodes = size
                all_cells_array = numpy.array([], dtype='int')

                #Reading Cells geometry data from CPO and writing them to IDS database in more orderly form
                for j in range(size):
                    node_idx = numpy.array([0,0,0,0], dtype='int')
                    free_edge = numpy.array([0,0,0], dtype='int')
                    edge_idx = edge.grid.spaces[0].objects[2].boundary[j,0] - 1
                    free_edge[0] = edge.grid.spaces[0].objects[2].boundary[j, 1] - 1
                    free_edge[1] = edge.grid.spaces[0].objects[2].boundary[j, 2] - 1
                    free_edge[2] = edge.grid.spaces[0].objects[2].boundary[j, 3] - 1

                    node_idx[0] = edge.grid.spaces[0].objects[1].boundary[edge_idx,0]-1
                    last_idx = 1
                    node_idx[last_idx] = edge.grid.spaces[0].objects[1].boundary[edge_idx, 1]-1

                    for loop_count in range (0,3):
                        if last_idx < 3:
                            for m in range(0,3):
                                edge_idx = free_edge[m]
                                if edge_idx < 0:
                                    continue
                                node1 = edge.grid.spaces[0].objects[1].boundary[edge_idx, 0]-1
                                node2 = edge.grid.spaces[0].objects[1].boundary[edge_idx, 1]-1
                                # print "{0:4}".format(j),'{0:5}'.format(node1), '{0:5}'.format(node2), '{0:5}'.format(edge_idx), '{0:5}'.format(loop_count)

                                if node_idx[last_idx] == node1:
                                    free_edge[m] = -1
                                    last_idx += 1
                                    node_idx[last_idx]= node2
                                    break
                                if node_idx[last_idx] == node2:
                                    free_edge[m] = -1
                                    last_idx += 1
                                    node_idx[last_idx] = node1
                                    break
                        assert loop_count < 3
                        loop_count += 1
                    # print k, "| ", node_idx[0], " | ", node_idx[1], " | ",node_idx[2], " | ",node_idx[3]
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_2_id].nodes[j]= node_idx[0] + 1 #+1 because of Fortran indices notation (1,2,3...)
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_2_id].nodes[size + j] = node_idx[1] + 1
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_2_id].nodes[2*size + j] = node_idx[2] + 1
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_2_id].nodes[3*size + j] = node_idx[3] + 1
            else:
                for j in range(size):

                    k = index_array[j]
                    # we get our indices from previously found cell nodes in "Cells" subgrid
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_2_id].nodes[j] = \
                        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[0].nodes[k]
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_2_id].nodes[size+j] = \
                        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[0].nodes[num_cells_all_nodes+k]
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_2_id].nodes[2*size+j] = \
                        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[0].nodes[2*num_cells_all_nodes+k]
                    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[obj_class_2_id].nodes[3*size+j] = \
                        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[subgrid_class].object[0].nodes[3*num_cells_all_nodes+k]

    # ----Allocating space and writing to IDS: Electron density--- #
    num_ne_subgrid = len(edge.fluid.ne.value)
    num_ne_subgrid_new = num_ne_subgrid + 4 # +4 because we have to add also values for new subgrids CORE, SOL, inner divertor and outer divertor
    # print "num_ne_subgrid: ", num_ne_subgrid
    # print "num_ne_subgrid_new: ", num_ne_subgrid_new
    imas_obj.edge_profiles.ggd[0].electrons.density.resize(num_ne_subgrid_new)


    for j in range(num_ne_subgrid_new):
        if j < num_ne_subgrid:
            imas_obj.edge_profiles.ggd[0].electrons.density[j].grid_subset_index = edge.fluid.ne.value[j].subgrid
            subgrid_base_id = imas_obj.edge_profiles.ggd[0].electrons.density[j].grid_subset_index
            ne_subgrid_scalar_size = len(edge.fluid.ne.value[j].scalar)
            # print "subgrid: ", subgrid_base_id, "ne_subgrid_scalar_size", ne_subgrid_scalar_size
            imas_obj.edge_profiles.ggd[0].electrons.density[j].values.resize(ne_subgrid_scalar_size)
            for k in range (ne_subgrid_scalar_size):
                imas_obj.edge_profiles.ggd[0].electrons.density[j].values[k] = edge.fluid.ne.value[j].scalar[k]
        else: # In CPO we got scalars for subgrids CORE, SOL, inner/outer divertor via indices, stored in edge.grid.subgrids[:].list[0].ind[0...N].
              # In IDS we don't have this option (or it was not found) and one way is to directly store scalars under new electron density subgrids values
            imas_obj.edge_profiles.ggd[0].electrons.density[j].grid_subset_index = j+2 #TRY TO FIX THIS! currently this goes hand in hand, but for the long run it should be done in a better way

            subgrid_base_id = imas_obj.edge_profiles.ggd[0].electrons.density[j].grid_subset_index
            ne_subgrid_scalar_size = len(edge.grid.subgrids[j+2-1].list[0].ind)
            # print "subgrid: ", subgrid_base_id, "ne_subgrid_scalar_size", ne_subgrid_scalar_size
            imas_obj.edge_profiles.ggd[0].electrons.density[j].values.resize(ne_subgrid_scalar_size)
            for k in range(ne_subgrid_scalar_size):
                scalar_index = edge.grid.subgrids[j + 2 - 1].list[0].ind[k]-1 #becuase we read indices in Fortran notation, we need indicesin C++ notation
                imas_obj.edge_profiles.ggd[0].electrons.density[j].values[k] = edge.fluid.ne.value[0].scalar[scalar_index]
                # print scalar_index, " ", imas_obj.edge_profiles.ggd[0].electrons.density[j].values[k]

    # ----Allocating space and writing to IDS: Electron temperature--- #
    num_te_subgrid = len(edge.fluid.te.value)
    num_te_subgrid_new = num_te_subgrid + 4 # +4 because we have to add also values for new subgrids CORE, SOL, inner divertor and outer divertor
    # print "num_te_subgrid: ", num_te_subgrid
    # print "num_te_subgrid_new: ", num_te_subgrid_new
    imas_obj.edge_profiles.ggd[0].electrons.temperature.resize(num_te_subgrid_new)
    for j in range(num_te_subgrid_new):
        if j < num_te_subgrid:
            imas_obj.edge_profiles.ggd[0].electrons.temperature[j].grid_subset_index = edge.fluid.te.value[j].subgrid
            te_subgrid_scalar_size = len(edge.fluid.te.value[j].scalar)
            subgrid_base_id = imas_obj.edge_profiles.ggd[0].electrons.temperature[j].grid_subset_index
            # print "subgrid: ", subgrid_base_id, "te_subgrid_scalar_size", te_subgrid_scalar_size
            imas_obj.edge_profiles.ggd[0].electrons.temperature[j].values.resize(te_subgrid_scalar_size)
            for k in range(te_subgrid_scalar_size):
                imas_obj.edge_profiles.ggd[0].electrons.temperature[j].values[k] = edge.fluid.te.value[j].scalar[k]
        else:  # In CPO we got scalars for subgrids CORE, SOL, inner/outer divertor via indices, stored in edge.grid.subgrids[:].list[0].ind[0...N].
            # In IDS we don't have this option (or it was not found) and one way is to directly store scalars under new electron temperature subgrids values
            imas_obj.edge_profiles.ggd[0].electrons.temperature[j].grid_subset_index = j + 2  # TRY TO FIX THIS! currently this goes hand in hand, but for the long run it should be done in a better way
            subgrid_base_id = imas_obj.edge_profiles.ggd[0].electrons.temperature[j].grid_subset_index
            te_subgrid_scalar_size = len(edge.grid.subgrids[j + 2 - 1].list[0].ind)
            # print "subgrid: ", subgrid_base_id, "te_subgrid_scalar_size", te_subgrid_scalar_size
            imas_obj.edge_profiles.ggd[0].electrons.temperature[j].values.resize(te_subgrid_scalar_size)
            for k in range(te_subgrid_scalar_size):
                scalar_index = edge.grid.subgrids[j + 2 - 1].list[0].ind[k] - 1  # becuase we read indices in Fortran notation, we need indicesin C++ notation
                imas_obj.edge_profiles.ggd[0].electrons.temperature[j].values[k] = edge.fluid.te.value[0].scalar[scalar_index]
                # print scalar_index, " ", imas_obj.edge_profiles.ggd[0].electrons.temperature[j].values[k]

    # ----Allocating space and writing to IDS: Ion density--- #
    ni_species_num = len(edge.fluid.ni)
    ##print "ni_species_num", ni_species_num
    imas_obj.edge_profiles.ggd[0].ion.resize(ni_species_num)
    for n in range (ni_species_num):
        num_ni_subgrid = len(edge.fluid.ni[n].value)
        num_ni_subgrid_new = num_ni_subgrid + 4 # +4 because we have to add also values for new subgrids CORE, SOL, inner divertor and outer divertor
        # print "num_ni_subgrid: ", num_ni_subgrid
        # print "num_ni_subgrid_new", num_ni_subgrid_new

        # Writing ion charge names
        imas_obj.edge_profiles.ggd[0].ion[n].label = edge.species[n].label

        imas_obj.edge_profiles.ggd[0].ion[n].density.resize(num_ni_subgrid_new)
        for j in range(num_ni_subgrid_new):
            if j < num_ni_subgrid:
                imas_obj.edge_profiles.ggd[0].ion[n].density[j].grid_subset_index = edge.fluid.ni[n].value[j].subgrid
                ni_subgrid_scalar_size = len(edge.fluid.ni[n].value[j].scalar)
                subgrid_base_id = imas_obj.edge_profiles.ggd[0].ion[n].density[j].grid_subset_index
                # print "subgrid: ", subgrid_base_id, "ni_subgrid_scalar_size", ni_subgrid_scalar_size
                imas_obj.edge_profiles.ggd[0].ion[n].density[j].values.resize(ni_subgrid_scalar_size)
                for k in range(ni_subgrid_scalar_size):
                    imas_obj.edge_profiles.ggd[0].ion[n].density[j].values[k] = edge.fluid.ni[n].value[j].scalar[k]
            else: # In CPO we got scalars for subgrids CORE, SOL, inner/outer divertor via indices, stored in edge.grid.subgrids[:].list[0].ind[0...N]
                # In IDS we don't have this option (or it was not found) and one way is to directly store scalars under new ion density subgrids values
                imas_obj.edge_profiles.ggd[0].ion[n].density[j].grid_subset_index = j + 2  # TRY TO FIX THIS! currently this goes hand in hand, but for the long run it should be done in a better way
                subgrid_base_id = imas_obj.edge_profiles.ggd[0].ion[n].density[j].grid_subset_index
                ni_subgrid_scalar_size = len(edge.grid.subgrids[j + 2 - 1].list[0].ind)
                # print "subgrid: ", subgrid_base_id, "ni_subgrid_scalar_size", ni_subgrid_scalar_size
                imas_obj.edge_profiles.ggd[0].ion[n].density[j].values.resize(ni_subgrid_scalar_size)
                for k in range(ni_subgrid_scalar_size):
                    scalar_index = edge.grid.subgrids[j + 2 - 1].list[0].ind[k] - 1  # becuase we read indices in Fortran notation, we need indices in C++ notation
                    imas_obj.edge_profiles.ggd[0].ion[n].density[j].values[k] = edge.fluid.ni[n].value[0].scalar[scalar_index]
                    # print scalar_index, " ", imas_obj.edge_profiles.ggd[0].ion[n].density[j].values[k]

    # ----Allocating space and writing to IDS: Ion temperature--- #
    ti_species_num = len(edge.fluid.ti)
    ##print "ti_species_num", ti_species_num
    #imas_obj.edge_profiles.ggd[0].ion.resize(ti_species_num)  #it has the same number of species as ion density

    for n in range(ti_species_num):
        num_ti_subgrid = len(edge.fluid.ti[n].value)
        num_ti_subgrid_new = num_ti_subgrid + 4 # +4 because we have to add also values for new subgrids CORE, SOL, inner divertor and outer divertor
        # print "num_ti_subgrid: ", num_ti_subgrid
        # print "num_ti_subgrid_new: ", num_ti_subgrid_new
        imas_obj.edge_profiles.ggd[0].ion[n].temperature.resize(num_te_subgrid_new)
        for j in range(num_ti_subgrid_new):
            if j < num_ti_subgrid:
                imas_obj.edge_profiles.ggd[0].ion[n].temperature[j].grid_subset_index = edge.fluid.ti[n].value[j].subgrid
                ti_subgrid_scalar_size = len(edge.fluid.ti[n].value[j].scalar)
                subgrid_base_id = imas_obj.edge_profiles.ggd[0].ion[n].temperature[j].grid_subset_index
                # print "subgrid: ", subgrid_base_id, "ti_subgrid_scalar_size", ti_subgrid_scalar_size
                imas_obj.edge_profiles.ggd[0].ion[n].temperature[j].values.resize(ti_subgrid_scalar_size)
                for k in range(ti_subgrid_scalar_size):
                    imas_obj.edge_profiles.ggd[0].ion[n].temperature[j].values[k] = edge.fluid.ti[n].value[j].scalar[k]
            else: # In CPO we got scalars for subgrids CORE, SOL, inner/outer divertor via indices, stored in edge.grid.subgrids[:].list[0].ind[0...N]
                # In IDS we don't have this option (or it was not found) and one way is to directly store scalars under new ion temperature subgrids values
                imas_obj.edge_profiles.ggd[0].ion[n].temperature[j].grid_subset_index = j + 2  # TRY TO FIX THIS! currently this goes hand in hand, but for the long run it should be done in a better way
                subgrid_base_id = imas_obj.edge_profiles.ggd[0].ion[n].temperature[j].grid_subset_index
                ti_subgrid_scalar_size = len(edge.grid.subgrids[j + 2 - 1].list[0].ind)
                # print "subgrid: ", subgrid_base_id, "ti_subgrid_scalar_size", ti_subgrid_scalar_size
                imas_obj.edge_profiles.ggd[0].ion[n].temperature[j].values.resize(ti_subgrid_scalar_size)
                for k in range(ti_subgrid_scalar_size):
                    scalar_index = edge.grid.subgrids[j + 2 - 1].list[0].ind[k] - 1  # becuase we read indices in Fortran notation, we need indices in C++ notation
                    imas_obj.edge_profiles.ggd[0].ion[n].temperature[j].values[k] = edge.fluid.ti[n].value[0].scalar[scalar_index]
                    # print scalar_index, " ", imas_obj.edge_profiles.ggd[0].ion[n].density[j].values[k]

    imas_obj.edge_profiles.putSlice()
    imas_obj.close()

def read_edge_cpo():

        cpo = ual.itm(shot, run, 'edge')
        cpo.open_env(user, tokamak, version)

        if cpo.isConnected():
                print 'open OK!'
        else:
                print 'open FAILED!'
                sys.exit()

        cpo.edgeArray.get()
        return cpo.edgeArray.array[0]

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], "srutvh", ["shot=", "run=", "user=", "tokamak=", "version=", "help"])

        for opt, arg in opts:
            #print opt, arg
            if opt in ("-s", "--shot"):
                shot = int(arg)
            elif opt in ("-r", "--run"):
                run = int(arg)
            elif opt in ("-u", "--user"):
                user = arg
            elif opt in ("-t", "--tokamak"):
                tokamak = arg
            elif opt in ("-v", "--version"):
                version = arg

            if opt in ("-h", "--help"):
                print "In order to run cpo2ids shot, run, user, tokamak and version variables must be defined. Example (terminal): python cpo2ids.py --shot=16151 --run=1000 --user=kosl --tokamak=aug --version=4.10a"
                try:
                    shot, run, user, tokamak, version
                except:
                    sys.exit()

    except getopt.GetoptError:
        print ('Supplied option not recognized!')
        print ('For help: cpo2ids.py -h / --help')
        sys.exit(2)

    edge = read_edge_cpo()

    write_ids()

    read_ids()







    
      

