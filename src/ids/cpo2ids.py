#! /usr/bin/env python
# cd /pfs/work/kosl/solps-iter-ids/modules/B2.5/src/ids
# Definition of the class structures in file imas.py
# module use -a ~dkaljun/imas/etc/modulefiles
# module avail imas
# module load imas/develop/3/ual/develop
# module display imas/develop/3/ual/develop
# imasdb solps-iter
# dd_doc
# imasdb -l # should show solps-iter with asterisk
# ls $UAL/pythonExamples
# cd ~/itmggd
# ls /pfs/home/kosl/itmggd/branches/4.10a/python/itmggd/examples
# less $UAL/python_pk/python2.6/ual/edge.py
# /pfs/work/kosl/bin/pycharm.sh &

#CPO -> IDS/IMAS
#edge
    #datainfo
	#dataprovider	-> /
	#putdate	-> /
	#source		-> /
	#comment	-> ids_properties/comment
	#cocos		-> /
	#id		-> /
	#isref		-> /
	#whatref	-> /
	    #user		-> /
	    #machine		-> /
	    #shot		-> /
	    #run		-> /
	    #occurence		-> /
	#putinfo	-> /
	    #putmethod		-> /
	    #putaccess		-> /
	    #putlocation	-> /
	    #rights		-> /
    #grid	-> ggd(:)/grid
	#uid		-> /
	#id		-> ggd(:)/grid/identifier/index
	#spaces	-> ggd(:)/grid/space(:)
	    #geotype	-> ggd(:)/grid/space(:)/geometry_type
	    #geotypeid	-> ggd(:)/grid/space(:)/geometry_type/index
	    #coordtype	-> ggd(:)/grid/space(:)/coordinates_type
	    #objects	-> ggd(:)/grid/space(:)/objects_per_dimension(:)/object(:)
		#boundary	-> ggd(:)/grid/space(:)/objects_per_dimension(:)/object(:)/boundary(:)
		#neighbour	-> ggd(:)/grid/space(:)/objects_per_dimension(:)/object(:)/boundary(:)/neighbours
		#geo		-> ggd(:)/grid/space(:)/objects_per_dimension(:)/object(:)/geometry
		#measure 	-> ggd(:)/grid/space(:)/objects_per_dimension(:)/object(:)/measure
	    #xpoints 	-> ggd(:)/grid/space(:)/objects_per_dimension(:)/object(:)/nodes  (??????)
	#subgrids-> ggd(:)/grid/grid_subset(:)
	    #id 	-> ggd(:)/grid/grid_subset(:)/identifier/index
	    #list 	-> ggd(:)/grid/grid_subset(:)/element(:)/object(:)	(??????)
		#cls		-> /
		#indset		-> /
		#ind		-> /
	#metric	-> ...


#NODES: edge::grid::spaces::objects(1) -> ggd(:)/grid/space(:)/objects_per_dimension(1)/object(:)
  #x: edge/grid/spaces/objects(1)/nodes/geo(i,0) 	) 
  #y: edge/grid/spaces/objects(1)/nodes/geo(i,1)	)-> ggd(:)/grid/space(:)/objects_per_dimension(1)/object(:)/geometry(x,y,z)
  #z: 0.0					 	)

#EDGES: edge::grid::spaces::objects(2) -> ggd(:)/grid/space(:)/objects_per_dimension(2)/object(:)
#CELLS: edge::grid::spaces::objects(3) -> ggd(:)/grid/space(:)/objects_per_dimension(3)/object(:)




import ual 
import imas
import numpy
import sys

'''
This sample program will open an existing pulse file (shot 8148, run 12) and will
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
    #my_ids_obj = ual.ids(8148, 12, 8148, 12)
    my_ids_obj = imas.ids(13, 1, 13, 1)

    my_ids_obj.open()  # Open the database

    if my_ids_obj.isConnected():
        print 'open OK!'
    else:
        print 'open FAILED!'
        sys.exit()

    my_ids_obj.edge_profiles.get()

    print '   my_ids_obj = '
    print my_ids_obj.edge_profiles




    my_ids_obj.close()

#read_ids()

def write_ids():
    shot = 13
    time = 1
    interp = 1

    imas_obj = imas.ids(1314,1,1314,1)

    imas_obj.create() #Create the data entry

    if imas_obj.isConnected():
        print 'Creation of data entry OK!'
    else:
        print 'Creation of data entry FAILED!'
        sys.exit()

    imas_obj.edge_profiles.ggd.resize(1)
    imas_obj.edge_profiles.ggd[0].grid.space.resize(1)
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension.resize(1)
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object.resize(1)

    #----test

    imas_obj.edge_profiles.ids_properties.homogeneous_time = 1 # Mandatory to define this property
    imas_obj.edge_profiles.ids_properties.comment = 'This is a test ids put using put_slice'
    imas_obj.edge_profiles.putNonTimed()
    imas_obj.edge_profiles.time.resize(1)


    imas_obj.edge_profiles.profiles_1d.resize(1)

    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry.resize(1, 2)
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry[0,0]= edge.grid.spaces[0].objects[0].geo[0,0]
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry[0,1] = edge.grid.spaces[0].objects[0].geo[0,1]

    print imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry.size
    print imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry[0, 0], \
    imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry[0, 1]
    imas_obj.edge_profiles.time[0] = 1.2

    print imas_obj.edge_profiles.time[0]
    imas_obj.edge_profiles.putSlice()

    #---


    # imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry.resize(num_nodes/2,2) #it creates (x,y) -> x*y elements = size
    #
    # print num_nodes
    # print imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry.size
    #
    #
    # imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry[0,0]= edge.grid.spaces[0].objects[0].geo[0,0]
    # imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry[0,1] = edge.grid.spaces[0].objects[0].geo[0,1]
    #
    # #ggd(:)/grid/space(:)/objects_per_dimension(1)/object(:)/geometry(x,y,z)
    #
    #
    # imas_obj.edge_profiles.time.resize(num_nodes/2);
    # print "time_size: ", imas_obj.edge_profiles.time.size
    #
    # imas_obj.edge_profiles.putNonTimed();
    #
    # for i in range (num_nodes/2): #3588 ok, 3589 not ok, 14352/4=3588
    #     imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry[i, 0] = edge.grid.spaces[0].objects[0].geo[i,0]
    #     imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry[i, 1] = edge.grid.spaces[0].objects[0].geo[i,1]
    #
    #     imas_obj.edge_profiles.time[i] = 0;
    #
    # imas_obj.edge_profiles.putSlice();
    #
    #
    # print imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry
    # print imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry[0,0], imas_obj.edge_profiles.ggd[0].grid.space[0].objects_per_dimension[0].object[0].geometry[0,1]



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

    #NODES: edge::grid::spaces::objects(1) -> ggd(:)/grid/space(:)/objects_per_dimension(1)/object(:)
      #x: edge/grid/spaces/objects(1)/nodes/geo(i,0) 	) 
      #y: edge/grid/spaces/objects(1)/nodes/geo(i,1)	)-> ggd(:)/grid/space(:)/objects_per_dimension(1)/object(:)/geometry(x,y,z)
      #z: 0.0					 	)
    
    num_nodes = edge.grid.spaces[0].objects[0].geo.size #it reads number of all elements, not nodes. NUmber of elements = number of nodes * 2
    # for i in range (num_nodes/2):
    #    #print x, y
    #    print edge.grid.spaces[0].objects[0].geo[i,0], edge.grid.spaces[0].objects[0].geo[i,1]

    # print edge.grid.spaces[0].objects[0].geo
    # print edge.grid.spaces[0].objects[0].geo[0, 0], edge.grid.spaces[0].objects[0].geo[0, 1]
    # print edge.grid.spaces[0].objects[0].geo[1, 0], edge.grid.spaces[0].objects[0].geo[1, 1]
    # print edge.grid.spaces[0].objects[0].geo[3587, 0], edge.grid.spaces[0].objects[0].geo[3587, 1]
    # print edge.grid.spaces[0].objects[0].geo[3588, 0], edge.grid.spaces[0].objects[0].geo[3588, 1]
    write_ids()
    #read_ids()
    
    
      

