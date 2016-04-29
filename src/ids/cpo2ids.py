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
# imasdb solps-iter
# python cpo2ids.py
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

#EDGES: edge::grid::spaces::objects(1) -> ggd(:)/grid/space(:)/objects_per_dimension(0)/object(1)
#FACES: edge::grid::spaces::objects(2) -> ggd(:)/grid/space(:)/objects_per_dimension(0)/object(2)

import ual 
import imas
import numpy
import sys

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
    # my_ids_obj = ual.ids(8148, 12, 8148, 12)
    my_ids_obj = imas.ids(16151, 1000, 16151, 1000)

    my_ids_obj.open()  # Open the database

    if my_ids_obj.isConnected():
        print 'open OK!'
    else:
        print 'open FAILED!'
        sys.exit()

    my_ids_obj.edge_profiles.get()

    # Example for reading data from IDS and CPO
    # print 'IDS'
    # print my_ids_obj.edge_profiles.ggd[0].electrons.temperature[4]
    # print 'CPO'
    # print edge.fluid.te.value[4]

    my_ids_obj.close()



def write_ids():
    time = 1
    interp = 1

    imas_obj = imas.ids(16151, 1000, 16151, 1000)

    imas_obj.create() #Create the data entry

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
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension.resize(1)
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object.resize(3)
    imas_obj.edge_profiles.ids_properties.homogeneous_time = 1 # Mandatory to define this property

    imas_obj.edge_profiles.ggd[0].grid.space[0].coordinates_type.resize(1)

    # ----Allocating space and writing to IDS: Geometry and Nodes---- #
    num_nodes = len(edge.grid.spaces[0].objects[0].geo)
    print 'CPO num_nodes:', num_nodes
    imas_obj.edge_profiles.time.resize(num_nodes*2);
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry.resize(num_nodes*2) #it should create (x,y) -> x*y elements = size
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].nodes.resize(num_nodes)

    for i in range(num_nodes):  # There are num_nodes*2 = 7176 geometry entries. [x1, y1, x2, y2, ..., xn, yn]
        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry[i * 2] = \
            edge.grid.spaces[0].objects[0].geo[i, 0]
        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry[(i * 2) + 1] = \
            edge.grid.spaces[0].objects[0].geo[i, 1]

        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].nodes[i] = i

        imas_obj.edge_profiles.time[i] = time;

    # ----Allocating space and writing to IDS: Edges---- #
    num_edges = len(edge.grid.spaces[0].objects[1].boundary)
    print 'CPO num_edges:', num_edges
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[1].boundary.resize(1)
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[1].boundary[0].neighbours.resize(num_edges*2)

    for i in range(num_edges):
        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[1].boundary[0].neighbours[i * 2] = \
            edge.grid.spaces[0].objects[1].boundary[i, 0]
        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[1].boundary[0].neighbours[
            (i * 2) + 1] = \
            edge.grid.spaces[0].objects[1].boundary[i, 1]

    # ----Allocating space and writing to IDS: Faces---- #
    num_faces = len(edge.grid.spaces[0].objects[2].boundary)
    print 'CPO num_faces:', num_faces
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[2].boundary.resize(1)
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[2].boundary[0].neighbours.resize(num_faces*4)

    for i in range(num_faces):
        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[2].boundary[0].neighbours[i * 4] = \
            edge.grid.spaces[0].objects[2].boundary[i, 0]
        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[2].boundary[0].neighbours[
            (i * 4) + 1] = \
            edge.grid.spaces[0].objects[2].boundary[i, 1]
        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[2].boundary[0].neighbours[
            (i * 4) + 2] = \
            edge.grid.spaces[0].objects[2].boundary[i, 2]
        imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[2].boundary[0].neighbours[
            (i * 4) + 3] = \
            edge.grid.spaces[0].objects[2].boundary[i, 3]

    # ----Allocating space and writing to IDS: Electron density---- #
    num_ne_subset = len(edge.fluid.ne.value)
    imas_obj.edge_profiles.ggd[0].electrons.density.resize(num_ne_subset)

    for i in range(num_ne_subset):
        num_ne_scalars = edge.fluid.ne.value[i].scalar.size
        imas_obj.edge_profiles.ggd[0].electrons.density[i].values.resize(num_ne_scalars)
        imas_obj.edge_profiles.ggd[0].electrons.density[i].grid_subset_index = edge.fluid.ne.value[i].subgrid
        for j in range(num_ne_scalars):
            imas_obj.edge_profiles.ggd[0].electrons.density[i].values[j] = edge.fluid.ne.value[i].scalar[j]

    # ----Allocating space and writing to IDS: Electron temperature---- #
    num_te_subset = len(edge.fluid.te.value)
    imas_obj.edge_profiles.ggd[0].electrons.temperature.resize(num_te_subset)
    for i in range(num_te_subset):
        num_te_scalars = edge.fluid.te.value[i].scalar.size
        imas_obj.edge_profiles.ggd[0].electrons.temperature[i].values.resize(num_te_scalars)
        imas_obj.edge_profiles.ggd[0].electrons.temperature[i].grid_subset_index = edge.fluid.te.value[i].subgrid

        for j in range(num_te_scalars):
            imas_obj.edge_profiles.ggd[0].electrons.temperature[i].values[j] = edge.fluid.te.value[i].scalar[j]

    # ----Allocating space and writing to IDS: Ion density ---- #
    num_ni_species = len(edge.fluid.ni)
    imas_obj.edge_profiles.ggd[0].ion.resize(num_ni_species)
    print 'num_ni_species:', num_ni_species
    for k in range(num_ni_species):
        num_ni_subset = len(edge.fluid.ni[k].value)
        imas_obj.edge_profiles.ggd[0].ion[k].density.resize(num_ni_subset)
        for i in range(num_ni_subset):
            num_ni_scalars = edge.fluid.ni[k].value[i].scalar.size
            imas_obj.edge_profiles.ggd[0].ion[k].density[i].values.resize(num_ni_scalars)
            imas_obj.edge_profiles.ggd[0].ion[k].density[i].grid_subset_index = edge.fluid.ni[k].value[i].subgrid

            for j in range(num_ni_scalars):
                imas_obj.edge_profiles.ggd[0].ion[k].density[i].values[j] = edge.fluid.ni[k].value[i].scalar[j]

    # ----Allocating space and writing to IDS: Ion temperature ---- #
    num_ti_species = len(edge.fluid.ti)
    # imas_obj.edge_profiles.ggd[0].ion.resize(num_ti_species) #It was already allocated in Ion density num_ni_species == num_ti_species
    for k in range(num_ti_species):
        num_ti_subset = len(edge.fluid.ti[k].value)
        imas_obj.edge_profiles.ggd[0].ion[k].temperature.resize(num_ti_subset)
        for i in range(num_ti_subset):
            num_ti_scalars = edge.fluid.ti[k].value[i].scalar.size
            imas_obj.edge_profiles.ggd[0].ion[k].temperature[i].values.resize(num_ti_scalars)
            imas_obj.edge_profiles.ggd[0].ion[k].temperature[i].grid_subset_index = edge.fluid.ti[k].value[i].subgrid

            for j in range(num_ti_scalars):
                imas_obj.edge_profiles.ggd[0].ion[k].temperature[i].values[j] = edge.fluid.ti[k].value[i].scalar[j]


    imas_obj.edge_profiles.putSlice();

    imas_obj.close()


def read_edge_cpo():
        cpo = ual.itm(16151, 1000, 'edge')
        #ual.itm(16151, 1000, 16151, 1000)
        cpo.open_env("kosl", "aug", "4.10a")
        

        if cpo.isConnected():
                print 'open OK!'
        else:
                print 'open FAILED!'
                sys.exit()

        cpo.edgeArray.get()
        return cpo.edgeArray.array[0]


if __name__ == '__main__':        
    edge = read_edge_cpo()

    write_ids()

    read_ids()







    
      

