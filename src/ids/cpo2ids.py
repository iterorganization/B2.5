#! /usr/bin/env python
"""
-----------------------------------------------------------------------------
DESCRIPTION
This utility converts edge CPO to IMAS edge_profiles IDS. It opens an
existing CPO, reads it and writes data to IDS.
Quick start:
Definition of the class structures is in file imas.py!
To run this converter use:
module load imas/3.7.4/ual/3.4.0
imasdb solps-iter
Run example (terminal):
     python cpo2ids.py --ishot=16151 --irun=1000 --iuser=g2penkod
     --idevice=aug --iversion=4.10a --oshot=16151 --orun=1000
     --ouser=g2penkod --odevice=solps-iter --oversion=3
Resulting files are stored under $HOME/public/imasdb/solps-iter/3/0
as ids_161511000.*
Further documentation and notes:
dd_doc
imasdb -l # should show solps-iter with asterisk
----------------------------------------------------------------------------
"""

import ual
import imas
import numpy
import sys
import getopt

ual.ualdef.DEVEL_DEBUG = 0
ual.ualdef.PRINT_DEBUG = 0
ual.ualdef.VERBOSE_DEBUG = 0


def read_ids():
    """Class IMAS is the main class for the UAL.

    It contains a set of field classes, each corresponding to an IDS
    defined in the UAL The parameters passed to this creator define the
    shot and run number. The second pair of arguments defines the
    reference shot and run and is used when the a new database is
    created, as in this example
    """

    print('Reading IDS: ')
    my_ids_obj = imas.ids(oshot, orun, oshot, orun)

    # my_ids_obj.open()  """ Open the database """
    my_ids_obj.open_env(ouser, odevice, oversion)

    if my_ids_obj.isConnected():
        print('IDS open OK!')
    else:
        print('IDS open FAILED!')
        sys.exit()

    my_ids_obj.edge_profiles.get()

    """ Example for reading data from IDS and CPO """
    print('IDS')
    # print(my_ids_obj.edge_profiles.ggd[0].grid.identifier.description)
    # print(my_ids_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0])
    # print my_ids_obj.edge_profiles.ggd[0].grid
    # print(my_ids_obj.edge_profiles.ggd[0])
    # print('CPO')
    my_ids_obj.close()

def write_ids():
    print('Writing IDS: ')
    time = 1
    interp = 1
    #treename = "ids"

    """--- CREATE NEW IDS ----"""
    imas_obj = imas.ids(oshot, orun, oshot, orun)

    # imas_obj.create() #Create the data entry
    imas_obj.create_env(ouser, odevice, oversion)

    if imas_obj.isConnected():
        print('Creation of data entry OK!')
    else:
        print('Creation of data entry FAILED!')
        sys.exit()

    """ "Shortcut" variable to imas_obj.edge_profiles data tree node"""
    edge_profiles=imas_obj.edge_profiles

    """Set and fill IDS with mandatory data"""
    edge_profiles.profiles_1d.resize(1)
    edge_profiles.ggd.resize(1)
    edge_profiles.putNonTimed()
    edge_profiles.ggd[0].grid.space.resize(1)
    edge_profiles.ggd[0].grid.space[0].objects_per_dimension.resize(3)
    edge_profiles.ids_properties.homogeneous_time = 1
    edge_profiles.time.resize(1)

    """Get number of coordinate types in CPO edge"""
    num_coordtype = len(edge.grid.spaces[0].coordtype)
    """Set coordinates_type"""
    edge_profiles.ggd[0].grid.space[0].coordinates_type.resize(num_coordtype)
    """Fill coordinates_type"""
    for j in range(num_coordtype):
        edge_profiles.ggd[0].grid.space[0].coordinates_type[j]= \
            edge.grid.spaces[0].coordtype[j, 0]

    """Put IDS grid description"""
    grid_description = "This is CPO" + \
        " shot=" + str(ishot) + " run=" + str(irun) + " user=" + str(iuser) + \
        " device=" + str(idevice) + " version=" + str(iversion) + \
        " converted to IDS" + \
        " shot=" + str(oshot) + " run=" + str(orun) + " user=" + str(ouser) + \
        " device=" + str(odevice) + " version=" + str(oversion) + "."

    edge_profiles.ggd[0].grid.identifier.description = grid_description

    """--- OBJECTS OF ALL DIMENSIONS (0D, 1D, 2D) ---"""
    num_dim = 3
    """Get number of 0D objects - grid nodes  out from grid subset
    holding list of all 0D objects"""
    #num_obj_0D_all = edge.grid.subgrids[1].list[0].indset[0].range[1]
    num_obj_0D_all = len(edge.grid.spaces[0].objects[0].geo)
    """Get number of 1D objects - edges  out from grid subset holding
    list of all 1D objects"""
    # num_obj_1D_all = edge.grid.subgrids[2].list[0].indset[0].range[1]
    num_obj_1D_all = len(edge.grid.spaces[0].objects[1].boundary)
    """Get number of 2D objects - cells/faces  out from grid subset
    holding list of all 2D objects"""
    # num_obj_2D_all = edge.grid.subgrids[0].list[0].indset[0].range[1]
    num_obj_2D_all = len(edge.grid.spaces[0].objects[2].boundary)

    """ "Shortcut" variable to ...space[0] node"""
    ids_space = edge_profiles.ggd[0].grid.space[0]
    """Set space for 0D objects - grid nodes"""
    ids_space.objects_per_dimension[0].object.resize(num_obj_0D_all)
    """Set space for 1D objects - edges"""
    ids_space.objects_per_dimension[1].object.resize(num_obj_1D_all)
    """Set space for 2D objects - faces\2D cells"""
    ids_space.objects_per_dimension[2].object.resize(num_obj_2D_all)

    """Set and fill 0D objects"""
    """ "Shortcut" variable to ...objects_per_dimension[0] node"""
    ids_dim_0D = \
        edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0]
    for i in range(num_obj_0D_all):
        """Set geometry of the 0D objects"""
        ids_dim_0D.object[i].geometry.resize(num_coordtype)
        """Set nodes of the 0Dobjects"""
        ids_dim_0D.object[i].nodes.resize(1)
        """Fill geometry of the 0D objects"""
        for j in range(num_coordtype):
            ids_dim_0D.object[i].geometry[j] = \
                edge.grid.spaces[0].objects[0].geo[i, j]
        """Fill nodes of the 0D objects """
        """(in Fortran notation. i_Fortran = i_Python + 1) """
        ids_dim_0D.object[i].nodes[0] = i + 1

    """Set and fill 1D objects"""
    """ "Shortcut" variable to ...objects_per_dimension[1] node"""
    ids_dim_1D = \
        edge_profiles.ggd[0].grid.space[0].objects_per_dimension[1]
    num_gridNodes_1D    = 2
    num_boundary_1D = 2
    for i in range(num_obj_1D_all):
        # Set nodes of the 1D objects
        ids_dim_1D.object[i].nodes.resize(num_gridNodes_1D)
        """Set boundary of the 1D objects"""
        ids_dim_1D.object[i].boundary.resize(num_boundary_1D)
        for j in range(num_boundary_1D):
            ids_dim_1D.object[i].boundary[j].index = \
                edge.grid.spaces[0].objects[1].boundary[i,j]
            # Fill nodes of the 1D objects
            ids_dim_1D.object[i].nodes[j] = \
                edge.grid.spaces[0].objects[1].boundary[i,j]

    """Set and fill 2D objects"""
    """ "Shortcut" variable to ...objects_per_dimension[2] node"""
    ids_dim_2D = \
        edge_profiles.ggd[0].grid.space[0].objects_per_dimension[2]
    obj_2d_all_array = numpy.array([], dtype='int')
    ids_dim_2D.object.resize(num_obj_2D_all)

    """Get 2D objects geometry (nodes) data from CPO, sort them into more
    orderly form  and put them into the IDS"""
    num_gridNodes_2D    = 4
    num_boundary_2D = 4
    # print("i, edge_idx, free_edge[0], free_edge[1], free_edge[2], node_idx[0], node_idx[last_idx]")
    for i in range(num_obj_2D_all):
        ids_dim_2D.object[i].nodes.resize(num_gridNodes_2D)

        node_idx = numpy.array([0,0,0,0], dtype='int')
        free_edge = numpy.array([0,0,0], dtype='int')
        edge_idx = edge.grid.spaces[0].objects[2].boundary[i,0] - 1
        free_edge[0] = edge.grid.spaces[0].objects[2].boundary[i, 1] - 1
        free_edge[1] = edge.grid.spaces[0].objects[2].boundary[i, 2] - 1
        free_edge[2] = edge.grid.spaces[0].objects[2].boundary[i, 3] - 1

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
        """Fill grid nodes of the 2D objects"""
        ids_dim_2D.object[i].nodes[0]= node_idx[0] + 1
        """+1 because of Fortran indices notation (1,2,3...)"""
        ids_dim_2D.object[i].nodes[1] = node_idx[3] + 1
        ids_dim_2D.object[i].nodes[2] = node_idx[2] + 1
        ids_dim_2D.object[i].nodes[3] = node_idx[1] + 1
        """Set boundary of the 2D objects"""
        ids_dim_2D.object[i].boundary.resize(num_boundary_2D)
        for j in range(num_boundary_2D):
            ids_dim_2D.object[i].boundary[j].index = \
                edge.grid.spaces[0].objects[2].boundary[i,j]

    """--- GRID SUBSETS ---"""

    """Get number of all existing grid subsets on CPO"""
    num_gridSubset = len(edge.grid.subgrids)

    """Set space for grid subsets"""
    edge_profiles.ggd[0].grid.grid_subset.resize(num_gridSubset)

    """Get number of grid_subset for each dimension from the CPO"""
    """Set placeholders"""
    num_gridSubset_dim_1 = 0
    num_gridSubset_dim_2 = 0
    num_gridSubset_dim_3 = 0

    for i in range(num_gridSubset):
        gridSubset_dim = edge.grid.subgrids[i].list[0].cls[0]
        if gridSubset_dim == 0:
            num_gridSubset_dim_1 += 1
        elif gridSubset_dim == 1:
            num_gridSubset_dim_2 += 1
        else:
            num_gridSubset_dim_3 += 1

    # print('Number of 0D grid subsets: ' + str(num_gridSubset_dim_1))
    # print('Number of 1D grid subsets: ' + str(num_gridSubset_dim_2))
    # print('Number of 2D grid subsets: ' + str(num_gridSubset_dim_3))

    for i in range(num_gridSubset):
        """"Shortcut" variable to grid_subset substructure"""
        ids_grid_subset = edge_profiles.ggd[0].grid.grid_subset[i]
        """Get index of the grid subset"""
        gridSubset_ind   = i + 1
        """Get name of the grid subset"""
        gridSubset_name = edge.grid.subgrids[i].id
        """Get class/dimension of objects which form the grid subset.
	In CPO edge its marked with 0,1,2..., in IDS we want it in
	Fortran notation """
        gridSubset_dim  = edge.grid.subgrids[i].list[0].cls[0] + 1

        """Put the grid subset name"""
        ids_grid_subset.identifier.name = gridSubset_name
        """Put the grid subset index"""
        ids_grid_subset.identifier.index = gridSubset_ind
        """Put the grid subset dimension"""
        ids_grid_subset.dimension = gridSubset_dim
        """Set elements of the grid subset"""
        """In this case grid subsets contain only 1 element"""
        ids_grid_subset.element.resize(1)

        """Get the object indices from either range of indices or explicit
        list of indices"""
        gridSubset_indset_found = len(edge.grid.subgrids[i].list[0].indset)

        ind_range_start = 0;
        ind_range_end   = 0;
        object_index = 0;
        num_obj_index = len(edge.grid.subgrids[i].list[0].ind)
        if gridSubset_indset_found > 0:
            object_index_range_size = len(edge.grid.subgrids[i].list[0]. \
                indset[0].range)
            if object_index_range_size > 1:
                range_found = 1
                ind_range_start = edge.grid.subgrids[i].list[0]. \
                    indset[0].range[0]-1
                ind_range_end = edge.grid.subgrids[i].list[0]. \
                    indset[0].range[1]
            else:
                print("Range is either EMPTY or it doesn't exist")
        else :
            range_found = 0
            ind_range_start  = edge.grid.subgrids[i].list[0].ind[0,0] - 1
            ind_range_end = edge.grid.subgrids[i].list[0]. \
                ind[num_obj_index-1,0] - 1

        gridSubset_object_index_array = numpy.array([], dtype='int')

        if range_found == 1:
            for j in range(ind_range_start, ind_range_end):
                gridSubset_object_index_array = numpy.insert(
                    gridSubset_object_index_array, j-ind_range_start, j)
        else:
            for j in range(0, num_obj_index):
                """python_index = Fortran_index - 1"""
                object_index = edge.grid.subgrids[i].list[0].ind[j,0] - 1
                gridSubset_object_index_array = numpy.insert(
                    gridSubset_object_index_array, j, object_index)

        """Set grid subset elements"""
        """IDS documentation: "One scalar value is provided per ELEMENT in """
        """the grid subset."""
        """In CPO some grid subsets have set values for them """
        """(1 node -> 1 value; 1 2D cell -> 1 value). """
        """Taking into account the IDS documentation, in our case our grid """
        """subsets consists of ELEMENTS containing 1 OBJECT (since we have """
        """n values available we need n elements)"""
        num_gridSubset_el = len(gridSubset_object_index_array)
        ids_grid_subset.element.resize(num_gridSubset_el)

        for j in range(num_gridSubset_el):
            """Set elements object for the grid subset element"""
            ids_grid_subset.element[j].object.resize(1)
            ids_grid_subset.element[j].object[0].space      = 1
            ids_grid_subset.element[j].object[0].dimension  = gridSubset_dim + 1
            ids_grid_subset.element[j].object[0].index = \
                gridSubset_object_index_array[j] + 1

    """---- Set and put to IDS: Electron density (ne) --- """
    num_ne_gridSubset = len(edge.fluid.ne.value)
    num_ne_gridSubset_new = num_ne_gridSubset + 4
    """+4 because we add also
    values for new gridSubsets
    CORE, SOL, inner divertor
    and outer divertor
    """
    edge_profiles.ggd[0].electrons.density.resize(num_ne_gridSubset_new)

    for j in range(num_ne_gridSubset_new):
        ids_ne = edge_profiles.ggd[0].electrons.density[j]
        if j < num_ne_gridSubset:

            ids_ne.grid_subset_index = edge.fluid.ne.value[j].subgrid
            gridSubset_index = ids_ne.grid_subset_index
            ne_gridSubset_scalar_size = len(edge.fluid.ne.value[j].scalar)
            ids_ne.values.resize(ne_gridSubset_scalar_size)
            for k in range (ne_gridSubset_scalar_size):
                ids_ne.values[k] = edge.fluid.ne.value[j].scalar[k]
        else:
            """In CPO we got scalars for gridSubsets CORE, SOL, inner/outer
            divertor via indices, stored in
            edge.grid.subgrids[:].list[0].ind[0...N].
            In IDS we don't have this option (or it was not found) and one
            way is to directly store scalars under new electron density
            gridSubsets values
            """
            ids_ne.grid_subset_index = j + 2
            """TODO: currently hardcoded,
            (to get grid subset
            indices 7,8,9,10)
            since this grid subset
            electron density values
            are not set in given CPO.
            Idea: combine search of grid
            subset name and getting the
            index
            """
            gridSubset_index = ids_ne.grid_subset_index
            ne_gridSubset_scalar_size = \
                len(edge.grid.subgrids[j + 2 - 1].list[0].ind)
            ids_ne.values.resize(ne_gridSubset_scalar_size)
            for k in range(ne_gridSubset_scalar_size):
                """Read scalar index in Fortran notation (1,2,3,...)
                but use it in python notation (0,1,2,...)
                """
                scalar_index = edge.grid.subgrids[j + 2 - 1].list[0].ind[k]
                ids_ne.values[k] = \
                    edge.fluid.ne.value[0].scalar[scalar_index - 1]

    """---- Set and put to IDS: Electron Temperature (te) ----"""
    num_te_gridSubset = len(edge.fluid.te.value)
    num_te_gridSubset_new = num_te_gridSubset + 4
    """ +4 because we add also
    values for new gridSubsets
    CORE, SOL, inner divertor
    and outer divertor
    """
    edge_profiles.ggd[0].electrons. \
        temperature.resize(num_te_gridSubset_new)
    for j in range(num_te_gridSubset_new):
        ids_te = edge_profiles.ggd[0].electrons.temperature[j]
        if j < num_te_gridSubset:
            ids_te.grid_subset_index = edge.fluid.te.value[j].subgrid
            te_gridSubset_scalar_size = len(edge.fluid.te.value[j].scalar)
            gridSubset_index = ids_te.grid_subset_index
            ids_te.values.resize(te_gridSubset_scalar_size)
            for k in range(te_gridSubset_scalar_size):
                ids_te.values[k] = edge.fluid.te.value[j].scalar[k]
        else:
            """In CPO we got scalars for gridSubsets CORE, SOL, inner/outer
            divertor via indices, stored in
            edge.grid.subgrids[:].list[0].ind[0...N].
            In IDS we don't have this option (or it was not found) and one
            way is to directly store scalars under new electron temperature
            gridSubsets values
            """
            ids_te.grid_subset_index = j + 2
            """TODO: currently hardcoded,
            (to get grid subset
            indices 7,8,9,10)
            since this grid subset
            electron density values
            are not set in given CPO.
            Idea: combine search of grid
            subset name and getting the
            index
            """
            gridSubset_index = ids_te.grid_subset_index
            te_gridSubset_scalar_size = \
                len(edge.grid.subgrids[j + 2 - 1].list[0].ind)
            ids_te.values.resize(te_gridSubset_scalar_size)
            for k in range(te_gridSubset_scalar_size):
                """Read scalar index in Fortran notation (1,2,3,...)
                but use it in python notation (0,1,2,...)
                """
                scalar_index = edge.grid.subgrids[j + 2 - 1].list[0].ind[k]
                ids_te.values[k] = \
                    edge.fluid.te.value[0].scalar[scalar_index - 1]

    """---- Set and put to IDS: Ion Density (ni) ----"""
    ni_species_num = len(edge.fluid.ni)
    """print "ni_species_num", ni_species_num"""
    edge_profiles.ggd[0].ion.resize(ni_species_num)
    for n in range (ni_species_num):
        num_ni_gridSubset = len(edge.fluid.ni[n].value)
        num_ni_gridSubset_new = num_ni_gridSubset + 4
        """+4 because we add also
        values for new gridSubsets
        CORE, SOL, inner divertor
        and outer divertor
        """

        ids_ion = edge_profiles.ggd[0].ion[n]

        """Put ion charge names. This label goes together with ion density and
        also ion temperature
        """
        ids_ion.label = edge.species[n].label

        ids_ion.density.resize(num_ni_gridSubset_new)
        for j in range(num_ni_gridSubset_new):
            ids_ni = edge_profiles.ggd[0].ion[n].density[j]
            if j < num_ni_gridSubset:
                ids_ni.grid_subset_index = edge.fluid.ni[n].value[j].subgrid
                ni_gridSubset_scalar_size = \
                    len(edge.fluid.ni[n].value[j].scalar)
                gridSubset_base_id = ids_ni.grid_subset_index
                ids_ni.values.resize(ni_gridSubset_scalar_size)
                for k in range(ni_gridSubset_scalar_size):
                    ids_ni.values[k] = edge.fluid.ni[n].value[j].scalar[k]
            else:
                """In CPO we got scalars for gridSubsets CORE, SOL,
                inner/outer divertor via indices, stored in
                edge.grid.gridSubsets[:].list[0].ind[0...N]
                In IDS we don't have this option (or it was not found)
                and one way is to directly store scalars under new ion
                density gridSubsets values
                """
                ids_ni.grid_subset_index = j + 2
                """TODO: currently hardcoded,
                (to get grid subset
                indices 7,8,9,10)
                since this grid subset
                electron density values
                are not set in given CPO.
                Idea: combine search of grid
                subset name and getting the
                index
                """
                gridSubset_base_id = ids_ni.grid_subset_index
                ni_gridSubset_scalar_size = \
                    len(edge.grid.subgrids[j + 2 - 1].list[0].ind)
                ids_ni.values.resize(ni_gridSubset_scalar_size)
                for k in range(ni_gridSubset_scalar_size):
                    """Read scalar index in Fortran notation (1,2,3,...)
                    but use it in python notation (0,1,2,...)
                    """
                    scalar_index = edge.grid.subgrids[j + 2 - 1].list[0].ind[k]
                    ids_ni.values[k] = \
                        edge.fluid.ni[n].value[0].scalar[scalar_index - 1]

    """---- Set and put to IDS: Ion Temperature (ni) ----"""
    ti_species_num = len(edge.fluid.ti)
    # edge_profiles.ggd[0].ion.resize(ti_species_num)
    """it has the same number of species as ion density"""

    for n in range(ti_species_num):
        num_ti_gridSubset = len(edge.fluid.ti[n].value)
        num_ti_gridSubset_new = num_ti_gridSubset + 4
        """+4 because we add also
        values for new gridSubsets
        CORE, SOL, inner divertor
        and outer divertor
        """
        edge_profiles.ggd[0].ion[n]. \
            temperature.resize(num_te_gridSubset_new)
        for j in range(num_ti_gridSubset_new):
            ids_ti = edge_profiles.ggd[0].ion[n].temperature[j]
            if j < num_ti_gridSubset:
                ids_ti.grid_subset_index = edge.fluid.ti[n].value[j].subgrid
                ti_gridSubset_scalar_size = len(edge.fluid.ti[n].value[j].scalar)
                gridSubset_base_id = ids_ti.grid_subset_index
                ids_ti.values.resize(ti_gridSubset_scalar_size)
                for k in range(ti_gridSubset_scalar_size):
                    ids_ti.values[k] = edge.fluid.ti[n].value[j].scalar[k]
            else:
                """In CPO we got scalars for gridSubsets CORE, SOL, inner/outer
                divertor via indices, stored in
                edge.grid.gridSubsets[:].list[0].ind[0...N]
                In IDS we don't have this option (or it was not found)
                and one way is to directly store scalars under new ion
                temperature gridSubsets values
                """
                ids_ti.grid_subset_index = j + 2
                """TODO: currently hardcoded,
                (to get grid subset
                indices 7,8,9,10)
                since this grid subset
                electron density values
                are not set in given CPO.
                Idea: combine search of grid
                subset name and getting the
                index
                """
                gridSubset_base_id = ids_ti.grid_subset_index
                ti_gridSubset_scalar_size = \
                    len(edge.grid.subgrids[j + 2 - 1].list[0].ind)
                ids_ti.values.resize(ti_gridSubset_scalar_size)
                for k in range(ti_gridSubset_scalar_size):
                    """Read scalar index in Fortran notation (1,2,3,...)
                    but use it in python notation (0,1,2,...)
                    """
                    scalar_index = \
                        edge.grid.subgrids[j + 2 - 1].list[0].ind[k]
                    ids_ti.values[k] =  \
                        edge.fluid.ti[n].value[0].scalar[scalar_index - 1]

    """Write put data do IDS"""
    edge_profiles.putSlice()
    """Close IDS"""
    imas_obj.close()

def read_edge_cpo():
        cpo = ual.itm(ishot, irun)
        cpo.open_env(iuser, idevice, iversion)

        if cpo.isConnected():
                print('CPO open OK!')
        else:
                print('CPO open FAILED!')
                sys.exit()

        cpo.edgeArray.get()
        return cpo.edgeArray.array[0]

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "srutvh", \
            ["ishot=", "irun=", "iuser=", "idevice=", "iversion=", \
             "oshot=", "orun=", "ouser=", "odevice=", "oversion=", \
            "help"])
        """ Input CPO parameters (shot, run, user, device, itm version) and
        output IDS parameters (shot, run, user, device, imas version)
        """
        for opt, arg in opts:
            if opt in ('-is', '--ishot'):
                ishot = int(arg)
            elif opt in ('-ir', '--irun'):
                irun = int(arg)
            elif opt in ('-iu', '--iuser'):
                iuser = arg
            elif opt in ('-id', '--idevice'):
                idevice = arg
            elif opt in ('-iv', '--iversion'):
                iversion = arg
            elif opt in ('-os', '--oshot'):
                oshot = int(arg)
            elif opt in ('-or', '--orun'):
                orun = int(arg)
            elif opt in ('-ou', '--ouser'):
                ouser = arg
            elif opt in ('-od', '--odevice'):
                odevice = arg
            elif opt in ('-ov', '--oversion'):
                oversion = arg

            if opt in ('-h', '--help'):
                print("In order to run cpo2ids input CPO and output IDS " + \
                    "shot, run, user, device and " + \
                    "version variables must be defined. Example (terminal): " + \
                    "python cpo2ids.py --ishot=16151 --irun=1000 " + \
                    "--iuser=g2penkod --idevice=aug --iversion=4.10a " + \
                    "--oshot=16151 --orun=1000 --ouser=g2penkod " + \
                    "--odevice=solps-iter --oversion=3" )
                try:
                    ishot, irun, iuser, idevice, iversion, \
                    oshot, orun, ouser, odevice, oversion
                except:
                    sys.exit()

    except getopt.GetoptError:
        print('Supplied option not recognized!')
        print('For help: cpo2ids.py -h / --help')
        sys.exit(2)

    edge = read_edge_cpo()

    write_ids()

    read_ids()
