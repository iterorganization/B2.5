Interface Data Structure for B2.5
=================================

B2.5 is a submodule of SOLPS-ITER project ([link](https://git.iter.org/projects/BND/repos/solps-iter/browse)).

b2_ual_write_b2mod
------------------

b2_ual_write_b2mod code is used to generate b2_ual_write_b2mod.exe
(main program), which is a post-processor for B2.
The code reads the plasma grid
geometry (full geometry descriptions of all available grid subsets)
and plasma state (electron density/temperature, ion temperature, velocity etc.).
The code then writes the obtained data to IDS database
with the use of b2mod scripts that utilize IMAS GGD Grid Service
Library routines.

#### Compiling and setting the environment:

In order compile the B2.5 writer code use the commands below while in
the SOLPS-ITER project main directory:

    $ tcsh
    $ source setup.csh
    $ cd modules/B2.5
    $ make ids

Note: At the time of writing this manual IMAS module
imas/3.13.0/ual/3.6.3 was used. This IMAS module provides also GGD support
as it includes IMAS GGD library routines (Fortran90).

Note: b2_ual_write and b2_ual_write_gsl are OUTDATED codes, but were left in
      the repository for documentation purposes and as and extra examples.

#### Documentation

To compile and see the B2.5 IDS documentation please run the following
commands in the terminal while in the SOLPS-ITER project main directory:

    $ tcsh
    $ source setup.csh
    $ cd $SOLPSTOP/modules/B2.5/documentation/
    $ make doc
    $ firefox source/Doxygen/html/index.html

#### Running the code:

The examples are available on ITER portal:
[B2.5 examples link](https://portal.iter.org/departments/POP/CM/IMAS/Forms/AllItems.aspx?RootFolder=%2Fdepartments%2FPOP%2FCM%2FIMAS%2FSOLPS-ITER%2FExamples)

In terminal navigate to directory containing the case required data files
(b2fgmtry, b2fstate etc.) and run the following command:

    $ $SOLPSTOP/modules/B2.5/builds/standalone.$HOST_NAME.$COMPILER/b2_ual_write_b2mod.exe --shot <shot> --run <run> --username <username> --device <device> --version <version> --step <number of steps>

The arguments marked with < ... > are the parameters of the IDS database
where the data is to be stored:
 - shot: The shot number of the database being created
 - run:  The run number of the database being created
 - username: Creator/owner of the IMAS IDS database
 - device: Device name of the IMAS IDS database (i. e. solps-iter, iter, aug)
 - version: Major version of the IMAS IDS database
 - step: Number of steps to be processed with b2mn_step() routine

Example of the command:

    $ $SOLPSTOP/modules/B2.5/builds/standalone.$HOST_NAME.$COMPILER/b2_ual_write_b2mod.exe --shot 1512 --run 6 --username penkod --device solps-iter --version 3 --step 250

CPO2IDS code
------------

Using this code is possible on  EUROFusion Gateway where both CPO and IDS support is
available.

Note: At the time of writing this manual IMAS module
imas/3.7.4/ual/3.4.0 was used.

#### Setting  the environment:

    $ module load imas/3.7.4/ual/3.4.0
    $ imasdb solps-iter

#### Running the code

    $ python cpo2ids.py --ishot=16151 --irun=1000 --iuser=g2penkod --idevice=aug --iversion=4.10a --oshot=16151 --orun=1000 --ouser=g2penkod --odevice=solps-iter --oversion=3"

where the arguments are:
- ishot = input IDS shot number
- irun = input IDS run number
- iuser = input IDS user name
- idevice = input IDS device name
- iversion = IMAS major version of input IDS
- oshot = output IDS shot number
- orun = output IDS run number
- ouser = output IDS user name
- odevice = output IDS device name
- oversion = IMAS major version of output IDS

Running the command

    $ python cpo2ids.py --help

will display the command example above.
