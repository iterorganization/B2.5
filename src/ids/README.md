Interface Data Structure for B2.5
=================================

b2_ual_write_b2mod
------------------

#### Compiling and setting the environment:

   $ cd $HOME/solps-iter
   $ tcsh
   $ source setup.csh
   $ cd modules/B2.5
   $ make ids

Note: At the time of writing this manual IMAS module
imas/3.13.0/ual/3.6.3 was used. This IMAS module provides also GGD support
as it includes IMAS GGD library routines (Fortran90).

Note: b2_ual_write and b2_ual_write_gsl are OUTDATED codes, but were left in
      the repository for documentation purposes and as and extra examples.

#### Running the code:

The examples are available on ITER portal:
[B2.5 examples](https://portal.iter.org/departments/POP/CM/IMAS/Forms/AllItems.aspx?RootFolder=%2Fdepartments%2FPOP%2FCM%2FIMAS%2FSOLPS-ITER%2FExamples)

In terminal navigate to directory containing the case required data files
(b2fgmtry, b2fstate etc.) and run the following command:

    $ $HOME/solps-iter/modules/B2.5/builds/standalone.ITER.ifort64/b2_ual_write_b2mod.exe

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
