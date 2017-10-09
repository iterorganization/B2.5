Interface Data Structure for B2.5
=================================

B2_ual_write and B2_ual_write_gsl code
-----------------

IMPORTANT: 
In order for b2_ual_write_gsl code to compile it requires beforehand locally  
compiled GSL library for Fortran90 in ${HOME}/ggd directory (ITER project 
repository [GGD](https://git.iter.org/projects/IMEX/repos/ggd), branch feature/IDS)!

Compile and set the environment:
   $ cd $HOME/solps-iter
   $ tcsh
   $ source setup.csh
   $ cd modules/B2.5
   $ make ids

Note: At the time of writing this manual imas module 
imas/3.8.0/ual/3.5.0 was used.

Use the script:

In terminal navigate to directory containing the b2fgmtry and b2fstate 
files and run the command

    $ $HOME/solps-iter/modules/B2.5/builds/standalone.ITER.ifort64/b2_ual_write.exe

or

	$ $HOME/solps-iter/modules/B2.5/builds/standalone.ITER.ifort64/b2_ual_write_gsl.exe


CPO2IDS code
------------

Using this code is possible on  EUROFusion Gateway where both CPO and IDS support is
available.

Set the environment: 
    $ module load imas/3.7.4/ual/3.4.0
    $ imasdb solps-iter

See Python code for usage or run:

    $ python cpo2ids.py --help

