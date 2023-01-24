B2.5 Algorithmic Differentiation documentation
======================

The source files contained in this folder are not compiled with any standard B2.5(-EIRENE) executable.
They are only needed as additional subroutines when a new version of the differentiated code is created.
The source files in tangent and adjoint folders are instead used to compile the tangent and adjoint AD codes respectively. They have been put in these subdirectories otherwise the Make file would try to use them for compilaton of standalone B2.5 or coupled B2.5-EIRENE executables.

The files with name "files_to_*.txt" are needed when differentiating the code, to tell Tapenade which source files
to keep/exclude from the differentiation. The splitting between files_to_exclude and files_to_ignore is needed by the automatic postprocessing of the differentiated files.

solpsGeneralLib: file where some basic functions/subroutine not shown to Tapenade are defined, only used for differentiation

Tangent AD
-------------------------------------------

Source files differentiated in multidirectional tangent mode for standalone B2.5. To compile and run:

    $ gmake listobj depend b25_tgt
    $ b2run -tgt b2mn

By default, running this version of the code will calculate the derivative of ANY cost function wrt three basic parameters
 - parm_dna(1): particle diffusion coefficient for species 1
 - parm_hce:  electron heat diffusion coefficient
 - parm_hci(1): ion heat diffusion coefficient for species 1

If any other coefficient is needed, then b2mn_d.F90 must be modified accordingly, as suggested by the commented example in the file itself, recalling that "nbdirs" must be set to the desired number of independent variables.

The code will automatically print in run.log the value of the cost function(s) and their derivatives.

b2mod_diffsizes.F acts similarly to DIMENSIONS.F and stores the maximum number of directions (nbdirsmax) allowed, i.e. the maximum number of indpendent variables with respect to which one can calculate the derivative of the cost function. Increase nbdirsmax and recompile if needed.

Current parameters for which derivative calculation is possible (default parameters):
 - Transport coefficients: parm_dna parm_hce parm_hci parm_vsa parm_sig parm_alf 
 - Radially-varying and ballooned transport coefficients: tdata b2tqna_ballooning b2tqna_ballooning_rescale
 - Boundary conditions: enepar enipar conpar mompar potpar enkpar 
 - Recycling coefficients: b2recyc 
 - k-model switches: keps_cd b2sikt_fac_sheath b2sikt_fac_sheath_core b2sikt_fac_diss b2sikt_fac_diss_core b2tfhi_fsigkt keps_heat keps_heat_i keps_sig keps_alf keps_visc keps_dkt keps_dzt keps_shear b2sikt_fac_vis_RS b2tfhi_fflokt b2tfhi_fconkt b2tfhi_fflozt b2tfhi_fconzt b2tfhi_fkt_hie b2tfhe_vis_kt 
 - Reaction rates: rtlsa rtlra rtlcx rtlqa
 - Optimization parameters: sigma par_opt_phys

To compile and run the optimization version, for example with PETCs/TAO optimization library:

    $ game listobj depend
    $ set_tao
    $ gmake b25_tgt
    $ b2run -opt_tgt b2mn

Note that b2mn_d.F90 need not be modified to run the optimization with a certain set of parameters! One need only to specify them with b2.optimization.parameters.

Adjoint AD
-------------------------------------------

Source files differentiated in adjoint mode for standalone B2.5. To compile and run:

    $ gmake listobj depend b25_adj
    $ b2run -adj b2mn

By default, running this version of the code will calculate the derivative of the FIRST cost function wrt all the default parameters indicated above
At each iteration however it will print out in run.log only the derivative of:
 - parm_dna(1): particle diffusion coefficient for species 1
 - parm_hce:  electron heat diffusion coefficient
 - parm_hci(1): ion heat diffusion coefficient for species 1
 - conpar(1,1,1): continuty BC for species 1 at boundary 1
 - enepar(1,1): electron heat BC at boundary 1
 - enipar(1,1): ion heat BC at boundary 1
 - b2recyc(1,2): recycling coefficient for species 1 at stratum 2

At the end of the run it will eventually print the gradient of all the default parameters for species [0-1], boundaries [1-6], and strata [1-5]. If more printout is needed during or at the end of the run then b2mod_driver_diff.F90 and b2mod_main_diff.F90 need to be modified accordingly.

To compile and run the optimization version, for example with PETCs/TAO optimization library:

    $ game listobj depend
    $ set_tao
    $ gmake b25_adj
    $ b2run -opt_adj b2mn

