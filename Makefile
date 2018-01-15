# Test whether necessary environment variables are defined; if not, exit
ifndef HOST_NAME
  $(error HOST_NAME not defined)
endif
ifndef COMPILER
  $(error COMPILER not defined)
endif

SRCB2   = ${PWD}
SRCDIR  = ${SRCB2}/src
DOCDIR  = ${SRCDIR}/documentation
PYTHON  = python

MAKES = ${SRCB2}/Makefile
# Include global SOLPS compiler settings
ifndef SOLPS_CPP
ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER} ] && echo yes || echo no ),yes)
  include ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}
  MAKES += ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}
  ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local ] && echo yes || echo no ),yes)
    include ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local
    MAKES += ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local
  endif
else
  $(warning ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER} not found.)
endif
else
  MAKES += ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}
  ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local ] && echo yes || echo no ),yes)
    include ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local
    MAKES += ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local
  endif
endif

ifdef USE_EIRENE
  ifndef SOLPSTOP
    $(error SOLPSTOP not defined, but trying to compile with Eirene)
  endif
endif
ifdef SOLPS_MPI
  ifndef USE_MPI
    USE_MPI = -DUSE_MPI
  endif
endif
ifdef SOLPS_OPENMP
  ifndef USE_OPENMP
    USE_OPENMP = -D_OPENMP
  endif
endif

# Default prefix for OBJDIR: standalone
PREF_OBJDIR = standalone
ifdef USE_EIRENE
PREF_OBJDIR = couple_SOLPS-ITER
endif

# Extension for OBJDIR if mpi and/or debug options are used
ifdef USE_MPI
EXT_MPI = .mpi
endif
ifdef USE_OPENMP
EXT_OPENMP = .openmp
endif
ifdef USE_IMPGYRO
EXT_IMPGYRO = .ig
endif
ifdef SOLPS_DEBUG
EXT_DEBUG = .debug
endif

# Directory where objectcode/binaries will be created
OBJDIR = ${SRCB2}/builds/${PREF_OBJDIR}.${HOST_NAME}.${COMPILER}${EXT_OPENMP}${EXT_MPI}${EXT_IMPGYRO}${EXT_DEBUG}

# If compiling with Eirene, look in default place for Eirene sources/lib
ifdef USE_EIRENE
  SRCEIR = ${SOLPSTOP}/modules/Eirene/src
  EIRDIR = ${SOLPSTOP}/modules/Eirene/builds/couple_SOLPS-ITER.${HOST_NAME}.${COMPILER}${EXT_MPI}${EXT_IMPGYRO}${EXT_DEBUG}
endif

ifeq ($(shell [ -e ${OBJDIR}/LISTOBJ ] && echo yes || echo no ),yes)
  include ${OBJDIR}/LISTOBJ
endif
  include ${SRCB2}/config/compile
  MAKES += ${SRCB2}/config/compile ${SRCB2}/config/config.${HOST_NAME}.${COMPILER}
ifeq ($(shell [ -e ${SRCB2}/config/config.${HOST_NAME}.${COMPILER}.local ] && echo yes || echo no ),yes)
  include ${SRCB2}/config/config.${HOST_NAME}.${COMPILER}.local
  MAKES+ = ${SRCB2}/config/config.${HOST_NAME}.${COMPILER}.local
endif

# Add external includes first
INCLUDE =
TAGSLIST =
ifdef NCDIR
INCLUDE += -I${NCDIR}/include
endif

ifdef USE_IMPGYRO
INCLUDE += ${MPI_CPP}
else
ifdef USE_MPI
INCLUDE += ${MPI_CPP}
else
INCLUDE += -I${SRCDIR}/mpi_dummy
endif
endif

ifdef MDSPLUS_DIR
ifndef LD_MDSPLUS
LD_MDSPLUS=-L${MDSPLUS_DIR}/lib -lMdsLib_client
endif
INCLUDE += -I${MDSPLUS_DIR}/include
endif

# If compiling with Paraview Catalyst
ifdef LD_CATALYST
SRCCAT = ${SRCDIR}/catalyst
INCLUDE += $(shell paraview-config --include)
endif

# Local includes
BASEDIR = ${OBJDIR}
INCLOCAL = ${SRCDIR}/include.local
MODLOCAL = ${SRCDIR}/modules.local
SRCLOCAL = ${SRCB2}/src.local
SOLPS4 = ${SOLPSTOP}/modules/solps4-5/src/solps4_B2
EIR4 = ${SOLPSTOP}/modules/solps4-5/src/Eirene_modules
ifeq ($(shell [ -d ${MODLOCAL} ] && echo yes || echo no ),yes)
INCLUDE += -I${MODLOCAL}
TAGSLIST += ${MODLOCAL}/*.F
endif
ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
INCLUDE += -I${SRCLOCAL}
TAGSLIST += ${SRCLOCAL}/*.F
endif
ifeq ($(shell [ -d ${INCLOCAL} ] && echo yes || echo no ),yes)
INCLUDE += -I${INCLOCAL}
TAGSLIST += ${INCLOCAL}/*.*
endif
INCLUDE += -I${SRCDIR}/common -I${SRCDIR}/include
SOLPS4INCLUDE = -I${SOLPSTOP}/modules/solps4-5/src/B2_include
TAGSLIST += ${SRCDIR}/include/*.* ${SRCDIR}/common/*.* ${SRCDIR}/common/COUPLE/*.F ${SRCDIR}/*/*.F ${SRCDIR}/*/*.F90 ${DOCDIR}/*.xml

DEFINES = ${B25_DEFINES} ${SOLPS_CPP}
ifdef USE_MPI
DEFINES += ${USE_MPI}
else
SOLPS4INCLUDE += -I${SOLPSTOP}/modules/solps4-5/src/Eirene_commons
endif
ifdef USE_OPENMP
DEFINES += ${USE_OPENMP}
endif
ifdef USE_IMPGYRO
DEFINES += ${USE_IMPGYRO}
endif
ifdef PERFMON
DEFINES += ${PERFMON}
endif
ifdef USE_EIRENE
DEFINES += ${USE_EIRENE}
endif
ifdef SOLPS_DEBUG
DEFINES += -DDBG
endif

# Switches to disable individual OpenMP regions, for debugging
# uncomment to activate
#
# It is not necessary to define all these flags to completely disable the
# OpenMP parallelization, in order to that, just compile without OpenMP
# compiler options (ifort -qopenmp or similar)

# DEFINES += -DNO_OPENMP_B2XPFE
# DEFINES += -DNO_OPENMP_B2SIFRTF
# DEFINES += -DNO_OPENMP_B2SIHS
# DEFINES += -DNO_OPENMP_B2SPCX
# DEFINES += -DNO_OPENMP_B2SPEL
# DEFINES += -DNO_OPENMP_B2SQCX
# DEFINES += -DNO_OPENMP_B2SQEL
# DEFINES += -DNO_OPENMP_B2SRDT
# DEFINES += -DNO_OPENMP_B2SRST
# DEFINES += -DNO_OPENMP_B2STCX
# DEFINES += -DNO_OPENMP_B2STEL
# DEFINES += -DNO_OPENMP_MYBLAS
# DEFINES += -DNO_OPENMP_SFILL
DEFINES += -DNO_OPENMP_B2NEWS_UNDERSCORE_LOOP1
DEFINES += -DNO_OPENMP_B2NEWS_UNDERSCORE_LOOP2
DEFINES += -DNO_OPENMP_B2NEWS_UNDERSCORE_LOOP3
DEFINES += -DNO_OPENMP_B2NEWS_UNDERSCORE_LOOP4
DEFINES += -DNO_OPENMP_B2NPMO

VHEAD =
ifeq ($(shell [ -d ${MODLOCAL} ] && echo yes || echo no ),yes)
VHEAD+=${MODLOCAL}:
endif
ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
VHEAD+=${SRCLOCAL}:
endif
# Required spacing variables for substitution
space :=
space +=
empty :=
#VPATH=$(subst $(space),$(empty),${VHEAD}${SRCDIR}/modules:${SRCDIR}/b2aux:${SRCDIR}/convert:${SRCDIR}/documentation:${SRCDIR}/driver:${SRCDIR}/equations:${SRCDIR}/input:${SRCDIR}/output:${SRCDIR}/postprocessing:${SRCDIR}/preprocessing:${SRCDIR}/solvers:${SRCDIR}/sources:${SRCDIR}/transport:${SRCDIR}/utility:${SRCDIR}/b2plot:${SRCDIR}/user)
VPATH=${VHEAD}${SRCDIR}/modules:${SRCDIR}/b2aux:${SRCDIR}/convert:${SRCDIR}/documentation:${SRCDIR}/driver:${SRCDIR}/equations:${SRCDIR}/input:${SRCDIR}/output:${SRCDIR}/postprocessing:${SRCDIR}/preprocessing:${SRCDIR}/solvers:${SRCDIR}/sources:${SRCDIR}/transport:${SRCDIR}/utility:${SRCDIR}/b2plot:${SRCDIR}/user
FPATH:=${VPATH}
VPATH += :${SRCDIR}/ids:${SRCDIR}/test
FFPATH += :${SRCDIR}/ids
FFPATH += :${SRCDIR}/modules

MODLIST =
MODLISTF =
MODLISTF90 =
ifeq ($(shell [ -d ${MODLOCAL} ] && echo yes || echo no ),yes)
MODLIST += ${MODLOCAL}/*.F
MODLISTF += ${MODLOCAL}/*.F
endif
ifdef LD_CATALYST
MODLIST += ${SRCDIR}/catalyst/*.F90
MODLISTF90 += ${SRCDIR}/catalyst/*.F90
endif
MODLIST += ${SRCDIR}/*/b2mod_*.F ${SRCDIR}/*/b2mod_*.F90 ${SRCDIR}/ids/*.F90
MODLISTF += ${SRCDIR}/*/b2mod_*.F
MODLISTF90 += ${SRCDIR}/*/b2mod_*.F90 ${SRCDIR}/ids/*.F90

ifeq ($(shell [ -d ${SOLPS4} ] && echo yes || echo no ),yes)
S4LIST = ${SOLPS4}/*.F
endif
ifeq ($(shell [ -d ${EIR4} ] && echo yes || echo no ),yes)
ifdef USE_EIRENE
E4LIST = ${EIR4}/precision.F ${EIR4}/parmmod.F ${EIR4}/braeir.F ${EIR4}/ccoupl.F ${EIR4}/clgin.F ${EIR4}/eirdiag.F ${EIR4}/ceirsrt.F
else
E4LIST = ${EIR4}/*.F
endif
endif

ifdef LD_CATALYST
ALLOBJS = ${OBJS:%.o=${OBJDIR}/%.o} ${OBJDIR}/cxxAdaptor.o
else
ALLOBJS = ${OBJS:%.o=${OBJDIR}/%.o}
endif

PROG_GE = b2pl.exe
PROG_GR = b2yg.exe b2yi.exe b2ym.exe b2yn.exe b2yp.exe b2yq.exe b2yr.exe b2ymb.exe b2yrp.exe b2ydm.exe
PROG_MN = b2mn.exe b2mnastra.exe
PROG_XD = b2xd.exe
PROG_OE = b2ag.exe b2co.exe b2fu.exe b2ts.exe b2uf.exe b2ye.exe b2yt.exe calc_atomic_data.exe
PROG_OT = b2ah.exe b2ai.exe b2ar.exe b2yi_gnuplot.exe b2yh.exe b2yv.exe b2fgmtry_mod.exe
#PROG_90 = check_b2_output.exe
PROG_OP = b2op.exe
PROG_OQ = b2mn_opt.exe
PROG_MD = b2md.exe b2rd.exe
PROG_ID = b2_ual_write.exe b2_ual_write_gsl.exe b2_ual_write_b2mod.exe

EXCLUDELIST = ${patsubst %.exe, %.o, ${PROG_GE} ${PROG_GR} ${PROG_MN} ${PROG_XD} ${PROG_OE} ${PROG_OT} ${PROG_90} ${PROG_MD} ${PROG_OP} ${PROG_OQ} ${PROG_ID}}
EXELIST = ${patsubst %.exe, %.o, ${PROG_GE} ${PROG_GR} ${PROG_MN} ${PROG_XD} ${PROG_OE} ${PROG_OT} ${PROG_MD} ${PROG_OP} ${PROG_OQ}}
EX90LIST = ${patsubst %.exe, %.o, ${PROG_90} ${PROG_ID}}

GEEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_GE}}
GREXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_GR}}
XDEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_XD}}
MNEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_MN}}
OEEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_OE}}
OTEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_OT}}
O9EXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_90}}
OPEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_OP}}
OQEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_OQ}}
MDEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_MD}}
IDEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_ID}}

.PHONY: DEFAULT NOPLOT ALL VERSION clean depend listobj tags echo local force

DEFAULT: VERSION ${MNEXE} ${OEEXE} ${OTEXE} ${O9EXE} ${GEEXE} ${GREXE}
ALL: VERSION ${MNEXE} ${OEEXE} ${OTEXE} ${O9EXE} ${GEEXE} ${GREXE} ${XDEXE}
NOPLOT: VERSION ${MNEXE} ${OEEXE} ${OTEXE} ${O9EXE}
ifdef MDSPLUS_DIR
DEFAULT: ${MDEXE}
ALL: ${MDEXE}
NOPLOT: ${MDEXE}
endif
ifdef PAR_OPT
DEFAULT: ${OPEXE} ${OQEXE}
ALL: ${OPEXE} ${OQEXE}
NOPLOT: ${OPEXE} ${OQEXE}
endif
ifdef IMAS_VERSION
DEFAULT: ${IDEXE}
ids: ${IDEXE}
NOPLOT: ${IDEXE}
endif
MAIN: VERSION ${MNEXE}

ifdef USE_EIRENE
VPATH+=${SRCEIR}/modules:${SRCEIR}/interfaces/couple_SOLPS-ITER
MODLIST+=${SRCEIR}/modules/*.f ${SRCEIR}/modules/*.[fF]90 ${SRCEIR}/interfaces/couple_SOLPS-ITER/eirmod_*.f ${SRCEIR}/interfaces/couple_SOLPS-ITER/eirmod_*.F90
MODLISTF+=${SRCEIR}/modules/*.f ${SRCEIR}/interfaces/couple_SOLPS-ITER/eirmod_*.f
MODLISTF90+=${SRCEIR}/modules/*.[fF]90 ${SRCEIR}/interfaces/couple_SOLPS-ITER/eirmod_*.F90
MNEXTRA=${EIRDIR}/libeirene.a ${EIRDIR}/libgr_dummy.a ${EIRDIR}/ioflush.o
else
# MNEXTRA=${EIRDIR}/eirmod_balanced_strategy.o ${EIRDIR}/eirmod_braeir.o ${EIRDIR}/eirmod_brascl.o ${EIRDIR}/eirmod_braspoi.o ${EIRDIR}/eirmod_cadgeo.o ${EIRDIR}/eirmod_cai.o ${EIRDIR}/eirmod_calstr_buffered.o ${EIRDIR}/eirmod_caprmc.o ${EIRDIR}/eirmod_ccflux.o ${EIRDIR}/eirmod_ccona.o ${EIRDIR}/eirmod_ccoupl.o ${EIRDIR}/eirmod_ccrm.o ${EIRDIR}/eirmod_cestim.o ${EIRDIR}/eirmod_cfplk.o ${EIRDIR}/eirmod_cgeom.o ${EIRDIR}/eirmod_cgrid.o ${EIRDIR}/eirmod_cgrptl.o ${EIRDIR}/eirmod_cinit.o ${EIRDIR}/eirmod_clast.o ${EIRDIR}/eirmod_clgin.o ${EIRDIR}/eirmod_clogau.o ${EIRDIR}/eirmod_comnnl.o ${EIRDIR}/eirmod_comprt.o ${EIRDIR}/eirmod_comsig.o ${EIRDIR}/eirmod_comsou.o ${EIRDIR}/eirmod_comspl.o ${EIRDIR}/eirmod_comusr.o ${EIRDIR}/eirmod_comxs.o ${EIRDIR}/eirmod_coutau.o ${EIRDIR}/eirmod_cpes.o ${EIRDIR}/eirmod_cpl3d.o ${EIRDIR}/eirmod_cplmsk.o ${EIRDIR}/eirmod_cplot.o ${EIRDIR}/eirmod_cpolyg.o ${EIRDIR}/eirmod_crand.o ${EIRDIR}/eirmod_crech.o ${EIRDIR}/eirmod_cref.o ${EIRDIR}/eirmod_crefmod.o ${EIRDIR}/eirmod_csdvi.o ${EIRDIR}/eirmod_csdvi_bgk.o ${EIRDIR}/eirmod_csdvi_cop.o ${EIRDIR}/eirmod_cspei.o ${EIRDIR}/eirmod_cspez.o ${EIRDIR}/eirmod_cstep.o ${EIRDIR}/eirmod_ctetra.o ${EIRDIR}/eirmod_ctext.o ${EIRDIR}/eirmod_ctrcei.o ${EIRDIR}/eirmod_ctrig.o ${EIRDIR}/eirmod_ctsurf.o ${EIRDIR}/eirmod_cupd.o ${EIRDIR}/eirmod_cvarusr.o ${EIRDIR}/eirmod_czt1.o ${EIRDIR}/eirmod_eirbra.o ${EIRDIR}/eirmod_eirdiag.o ${EIRDIR}/eirmod_infcop.o ${EIRDIR}/eirmod_module_avltree.o ${EIRDIR}/eirmod_mpi.o ${EIRDIR}/eirmod_octree.o ${EIRDIR}/eirmod_parmmod.o ${EIRDIR}/eirmod_precision.o
# EXCLUDELIST += ${patsubst ${OBJDIR}/%.o, %.o, ${MNEXTRA} }
endif
ifdef LD_CATALYST
VPATH += :${SRCDIR}/catalyst
FFPATH += :${SRCDIR}/catalyst
endif

MODULES = ${patsubst %.F %.f %.F90 %.f90,%.o,${shell echo ${MODLIST} } }
MODMODS = ${MODULES:%.o=${OBJDIR}/%.${MOD}}
MODOBJS = ${MODULES:%.o=${OBJDIR}/%.o}
SOLPS4OBJS = ${patsubst ${SOLPS4}/%.F,${OBJDIR}/%.o,${shell echo ${S4LIST} } }
EIR4MODS = ${patsubst ${EIR4}/%.F,${OBJDIR}/%.${MOD},${shell echo ${E4LIST} } }
EIR4OBJS = ${patsubst ${EIR4}/%.F,${OBJDIR}/%.o,${shell echo ${E4LIST} } }

ifeq (${MOD},o)
LIBOBJS = $(filter-out ${MODOBJS},${ALLOBJS})
else
LIBOBJS = ${ALLOBJS}
endif

${DOCDIR}/b2cdci.F: ${DOCDIR}/b2input.xml ${DOCDIR}/b2cdci.py
	-cd ${DOCDIR}; ${PYTHON} b2cdci.py || echo "! Error building b2cdci.F from b2input.xml" > ${DOCDIR}/b2cdci.F

${DOCDIR}/b2cdcn.F: ${DOCDIR}/b2input.xml ${DOCDIR}/b2cdcn.py
	-cd ${DOCDIR}; ${PYTHON} b2cdcn.py || echo "! Error building b2cdcn.F from b2input.xml" > ${DOCDIR}/b2cdcn.F

ifdef USE_EIRENE
${OBJDIR}/libgr_dummy.a:
	ln -sf ${EIRDIR}/libgr_dummy.a ${OBJDIR}

${OBJDIR}/libeirene.a:
	ln -sf ${EIRDIR}/libeirene.a ${OBJDIR}

${OBJDIR}/ioflush.o:
	ln -sf ${EIRDIR}/ioflush.o ${OBJDIR}

ifneq (${MOD},o)
${OBJDIR}/eirmod_extrab25.${MOD}:
	ln -sf ${EIRDIR}/eirmod_extrab25.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_balanced_strategy.${MOD}:
	ln -sf ${EIRDIR}/eirmod_balanced_strategy.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_braeir.${MOD}:
	ln -sf ${EIRDIR}/eirmod_braeir.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_brascl.${MOD}:
	ln -sf ${EIRDIR}/eirmod_brascl.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_braspoi.${MOD}:
	ln -sf ${EIRDIR}/eirmod_braspoi.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cadgeo.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cadgeo.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cai.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cai.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_calstr_buffered.${MOD}:
	ln -sf ${EIRDIR}/eirmod_calstr_buffered.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_caprmc.${MOD}:
	ln -sf ${EIRDIR}/eirmod_caprmc.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_ccflux.${MOD}:
	ln -sf ${EIRDIR}/eirmod_ccflux.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_ccona.${MOD}:
	ln -sf ${EIRDIR}/eirmod_ccona.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_ccoupl.${MOD}:
	ln -sf ${EIRDIR}/eirmod_ccoupl.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_ccrm.${MOD}:
	ln -sf ${EIRDIR}/eirmod_ccrm.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cestim.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cestim.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cfplk.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cfplk.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cgeom.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cgeom.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cgrid.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cgrid.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cgrptl.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cgrptl.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cinit.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cinit.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_clast.${MOD}:
	ln -sf ${EIRDIR}/eirmod_clast.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_clgin.${MOD}:
	ln -sf ${EIRDIR}/eirmod_clgin.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_clogau.${MOD}:
	ln -sf ${EIRDIR}/eirmod_clogau.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_clmsur.${MOD}:
	ln -sf ${EIRDIR}/eirmod_clmsur.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_comnnl.${MOD}:
	ln -sf ${EIRDIR}/eirmod_comnnl.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_comprt.${MOD}:
	ln -sf ${EIRDIR}/eirmod_comprt.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_comsig.${MOD}:
	ln -sf ${EIRDIR}/eirmod_comsig.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_comsou.${MOD}:
	ln -sf ${EIRDIR}/eirmod_comsou.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_comspl.${MOD}:
	ln -sf ${EIRDIR}/eirmod_comspl.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_comusr.${MOD}:
	ln -sf ${EIRDIR}/eirmod_comusr.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_comxs.${MOD}:
	ln -sf ${EIRDIR}/eirmod_comxs.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_coutau.${MOD}:
	ln -sf ${EIRDIR}/eirmod_coutau.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cpes.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cpes.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cpl3d.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cpl3d.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cplmsk.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cplmsk.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cplot.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cplot.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cpolyg.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cpolyg.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_crand.${MOD}:
	ln -sf ${EIRDIR}/eirmod_crand.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_crech.${MOD}:
	ln -sf ${EIRDIR}/eirmod_crech.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cref.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cref.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_crefmod.${MOD}:
	ln -sf ${EIRDIR}/eirmod_crefmod.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_csdvi.${MOD}:
	ln -sf ${EIRDIR}/eirmod_csdvi.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_csdvi_bgk.${MOD}:
	ln -sf ${EIRDIR}/eirmod_csdvi_bgk.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_csdvi_cop.${MOD}:
	ln -sf ${EIRDIR}/eirmod_csdvi_cop.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cspei.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cspei.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cspez.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cspez.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cstep.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cstep.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_ctetra.${MOD}:
	ln -sf ${EIRDIR}/eirmod_ctetra.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_ctext.${MOD}:
	ln -sf ${EIRDIR}/eirmod_ctext.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_ctrcei.${MOD}:
	ln -sf ${EIRDIR}/eirmod_ctrcei.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_ctrig.${MOD}:
	ln -sf ${EIRDIR}/eirmod_ctrig.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_ctsurf.${MOD}:
	ln -sf ${EIRDIR}/eirmod_ctsurf.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cupd.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cupd.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_cvarusr.${MOD}:
	ln -sf ${EIRDIR}/eirmod_cvarusr.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_czt1.${MOD}:
	ln -sf ${EIRDIR}/eirmod_czt1.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_eirbra.${MOD}:
	ln -sf ${EIRDIR}/eirmod_eirbra.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_eirdiag.${MOD}:
	ln -sf ${EIRDIR}/eirmod_eirdiag.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_infcop.${MOD}: ${OBJDIR}/eirmod_cplot.${MOD}
	ln -sf ${EIRDIR}/eirmod_infcop.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_module_avltree.${MOD}:
	ln -sf ${EIRDIR}/eirmod_module_avltree.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_mpi.${MOD}:
	ln -sf ${EIRDIR}/eirmod_mpi.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_octree.${MOD}:
	ln -sf ${EIRDIR}/eirmod_octree.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_parmmod.${MOD}:
	ln -sf ${EIRDIR}/eirmod_parmmod.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_precision.${MOD}:
	ln -sf ${EIRDIR}/eirmod_precision.${MOD} ${OBJDIR}
endif

${OBJDIR}/eirmod_extrab25.o:
	ln -sf ${EIRDIR}/eirmod_extrab25.o ${OBJDIR}

${OBJDIR}/eirmod_balanced_strategy.o:
	ln -sf ${EIRDIR}/eirmod_balanced_strategy.o ${OBJDIR}

${OBJDIR}/eirmod_braeir.o:
	ln -sf ${EIRDIR}/eirmod_braeir.o ${OBJDIR}

${OBJDIR}/eirmod_brascl.o:
	ln -sf ${EIRDIR}/eirmod_brascl.o ${OBJDIR}

${OBJDIR}/eirmod_braspoi.o:
	ln -sf ${EIRDIR}/eirmod_braspoi.o ${OBJDIR}

${OBJDIR}/eirmod_cadgeo.o:
	ln -sf ${EIRDIR}/eirmod_cadgeo.o ${OBJDIR}

${OBJDIR}/eirmod_cai.o:
	ln -sf ${EIRDIR}/eirmod_cai.o ${OBJDIR}

${OBJDIR}/eirmod_calstr_buffered.o:
	ln -sf ${EIRDIR}/eirmod_calstr_buffered.o ${OBJDIR}

${OBJDIR}/eirmod_caprmc.o:
	ln -sf ${EIRDIR}/eirmod_caprmc.o ${OBJDIR}

${OBJDIR}/eirmod_ccflux.o:
	ln -sf ${EIRDIR}/eirmod_ccflux.o ${OBJDIR}

${OBJDIR}/eirmod_ccona.o:
	ln -sf ${EIRDIR}/eirmod_ccona.o ${OBJDIR}

${OBJDIR}/eirmod_cestim.o:
	ln -sf ${EIRDIR}/eirmod_cestim.o ${OBJDIR}

${OBJDIR}/eirmod_ccoupl.o:
	ln -sf ${EIRDIR}/eirmod_ccoupl.o ${OBJDIR}

${OBJDIR}/eirmod_ccrm.o:
	ln -sf ${EIRDIR}/eirmod_ccrm.o ${OBJDIR}

${OBJDIR}/eirmod_cfplk.o:
	ln -sf ${EIRDIR}/eirmod_cfplk.o ${OBJDIR}

${OBJDIR}/eirmod_cgeom.o:
	ln -sf ${EIRDIR}/eirmod_cgeom.o ${OBJDIR}

${OBJDIR}/eirmod_cgrid.o:
	ln -sf ${EIRDIR}/eirmod_cgrid.o ${OBJDIR}

${OBJDIR}/eirmod_cgrptl.o:
	ln -sf ${EIRDIR}/eirmod_cgrptl.o ${OBJDIR}

${OBJDIR}/eirmod_cinit.o:
	ln -sf ${EIRDIR}/eirmod_cinit.o ${OBJDIR}

${OBJDIR}/eirmod_clast.o:
	ln -sf ${EIRDIR}/eirmod_clast.o ${OBJDIR}

${OBJDIR}/eirmod_clgin.o:
	ln -sf ${EIRDIR}/eirmod_clgin.o ${OBJDIR}

${OBJDIR}/eirmod_clogau.o:
	ln -sf ${EIRDIR}/eirmod_clogau.o ${OBJDIR}

${OBJDIR}/eirmod_clmsur.o:
	ln -sf ${EIRDIR}/eirmod_clmsur.o ${OBJDIR}

${OBJDIR}/eirmod_comnnl.o:
	ln -sf ${EIRDIR}/eirmod_comnnl.o ${OBJDIR}

${OBJDIR}/eirmod_comprt.o:
	ln -sf ${EIRDIR}/eirmod_comprt.o ${OBJDIR}

${OBJDIR}/eirmod_comsig.o:
	ln -sf ${EIRDIR}/eirmod_comsig.o ${OBJDIR}

${OBJDIR}/eirmod_comsou.o:
	ln -sf ${EIRDIR}/eirmod_comsou.o ${OBJDIR}

${OBJDIR}/eirmod_comspl.o:
	ln -sf ${EIRDIR}/eirmod_comspl.o ${OBJDIR}

${OBJDIR}/eirmod_comusr.o:
	ln -sf ${EIRDIR}/eirmod_comusr.o ${OBJDIR}

${OBJDIR}/eirmod_comxs.o:
	ln -sf ${EIRDIR}/eirmod_comxs.o ${OBJDIR}

${OBJDIR}/eirmod_coutau.o:
	ln -sf ${EIRDIR}/eirmod_coutau.o ${OBJDIR}

${OBJDIR}/eirmod_cpes.o:
	ln -sf ${EIRDIR}/eirmod_cpes.o ${OBJDIR}

${OBJDIR}/eirmod_cpl3d.o:
	ln -sf ${EIRDIR}/eirmod_cpl3d.o ${OBJDIR}

${OBJDIR}/eirmod_cplmsk.o:
	ln -sf ${EIRDIR}/eirmod_cplmsk.o ${OBJDIR}

${OBJDIR}/eirmod_cplot.o:
	ln -sf ${EIRDIR}/eirmod_cplot.o ${OBJDIR}

${OBJDIR}/eirmod_cpolyg.o:
	ln -sf ${EIRDIR}/eirmod_cpolyg.o ${OBJDIR}

${OBJDIR}/eirmod_crand.o:
	ln -sf ${EIRDIR}/eirmod_crand.o ${OBJDIR}

${OBJDIR}/eirmod_crech.o:
	ln -sf ${EIRDIR}/eirmod_crech.o ${OBJDIR}

${OBJDIR}/eirmod_cref.o:
	ln -sf ${EIRDIR}/eirmod_cref.o ${OBJDIR}

${OBJDIR}/eirmod_crefmod.o:
	ln -sf ${EIRDIR}/eirmod_crefmod.o ${OBJDIR}

${OBJDIR}/eirmod_csdvi.o:
	ln -sf ${EIRDIR}/eirmod_csdvi.o ${OBJDIR}

${OBJDIR}/eirmod_csdvi_bgk.o:
	ln -sf ${EIRDIR}/eirmod_csdvi_bgk.o ${OBJDIR}

${OBJDIR}/eirmod_csdvi_cop.o:
	ln -sf ${EIRDIR}/eirmod_csdvi_cop.o ${OBJDIR}

${OBJDIR}/eirmod_cspei.o:
	ln -sf ${EIRDIR}/eirmod_cspei.o ${OBJDIR}

${OBJDIR}/eirmod_cspez.o:
	ln -sf ${EIRDIR}/eirmod_cspez.o ${OBJDIR}

${OBJDIR}/eirmod_cstep.o:
	ln -sf ${EIRDIR}/eirmod_cstep.o ${OBJDIR}

${OBJDIR}/eirmod_ctetra.o:
	ln -sf ${EIRDIR}/eirmod_ctetra.o ${OBJDIR}

${OBJDIR}/eirmod_ctext.o:
	ln -sf ${EIRDIR}/eirmod_ctext.o ${OBJDIR}

${OBJDIR}/eirmod_ctrcei.o:
	ln -sf ${EIRDIR}/eirmod_ctrcei.o ${OBJDIR}

${OBJDIR}/eirmod_ctrig.o:
	ln -sf ${EIRDIR}/eirmod_ctrig.o ${OBJDIR}

${OBJDIR}/eirmod_ctsurf.o:
	ln -sf ${EIRDIR}/eirmod_ctsurf.o ${OBJDIR}

${OBJDIR}/eirmod_cupd.o:
	ln -sf ${EIRDIR}/eirmod_cupd.o ${OBJDIR}

${OBJDIR}/eirmod_cvarusr.o:
	ln -sf ${EIRDIR}/eirmod_cvarusr.o ${OBJDIR}

${OBJDIR}/eirmod_czt1.o:
	ln -sf ${EIRDIR}/eirmod_czt1.o ${OBJDIR}

${OBJDIR}/eirmod_eirbra.o:
	ln -sf ${EIRDIR}/eirmod_eirbra.o ${OBJDIR}

${OBJDIR}/eirmod_infcop.o: ${OBJDIR}/eirmod_cplot.o
	ln -sf ${EIRDIR}/eirmod_infcop.o ${OBJDIR}

${OBJDIR}/eirmod_module_avltree.o:
	ln -sf ${EIRDIR}/eirmod_module_avltree.o ${OBJDIR}

${OBJDIR}/eirmod_mpi.o:
	ln -sf ${EIRDIR}/eirmod_mpi.o ${OBJDIR}

${OBJDIR}/eirmod_octree.o:
	ln -sf ${EIRDIR}/eirmod_octree.o ${OBJDIR}

${OBJDIR}/eirmod_parmmod.o:
	ln -sf ${EIRDIR}/eirmod_parmmod.o ${OBJDIR}

${OBJDIR}/eirmod_precision.o:
	ln -sf ${EIRDIR}/eirmod_precision.o ${OBJDIR}
else
${OBJDIR}/eirmod_balanced_strategy.${MOD}:
	touch ${OBJDIR}/eirmod_balanced_strategy.${MOD}

${OBJDIR}/eirmod_braeir.${MOD}:
	touch ${OBJDIR}/eirmod_braeir.${MOD}

${OBJDIR}/eirmod_ccona.${MOD}:
	touch ${OBJDIR}/eirmod_ccona.${MOD}

${OBJDIR}/eirmod_ccoupl.${MOD}:
	touch ${OBJDIR}/eirmod_ccoupl.${MOD}

${OBJDIR}/eirmod_cestim.${MOD}:
	touch ${OBJDIR}/eirmod_cestim.${MOD}

${OBJDIR}/eirmod_cinit.${MOD}:
	touch ${OBJDIR}/eirmod_cinit.${MOD}

${OBJDIR}/eirmod_clogau.${MOD}:
	touch ${OBJDIR}/eirmod_clogau.${MOD}

${OBJDIR}/eirmod_comprt.${MOD}:
	touch ${OBJDIR}/eirmod_comprt.${MOD}

${OBJDIR}/eirmod_comusr.${MOD}:
	touch ${OBJDIR}/eirmod_comusr.${MOD}

${OBJDIR}/eirmod_coutau.${MOD}:
	touch ${OBJDIR}/eirmod_coutau.${MOD}

${OBJDIR}/eirmod_cpes.${MOD}:
	touch ${OBJDIR}/eirmod_cpes.${MOD}

${OBJDIR}/eirmod_ctrcei.${MOD}:
	touch ${OBJDIR}/eirmod_ctrcei.${MOD}

${OBJDIR}/eirmod_ctrig.${MOD}:
	touch ${OBJDIR}/eirmod_ctrig.${MOD}

${OBJDIR}/eirmod_eirbra.${MOD}:
	touch ${OBJDIR}/eirmod_eirbra.${MOD}

${OBJDIR}/eirmod_eirdiag.${MOD}:
	touch ${OBJDIR}/eirmod_eirdiag.${MOD}

${OBJDIR}/eirmod_extrab25.${MOD}:
	touch ${OBJDIR}/eirmod_extrab25.${MOD}

${OBJDIR}/eirmod_infcop.${MOD}:
	touch ${OBJDIR}/eirmod_infcop.${MOD}

${OBJDIR}/eirmod_mpi.${MOD}:
	touch ${OBJDIR}/eirmod_mpi.${MOD}

${OBJDIR}/eirmod_parmmod.${MOD}:
	touch ${OBJDIR}/eirmod_parmmod.${MOD}

${OBJDIR}/eirmod_precision.${MOD}:
	ln -s ${OBJDIR}/precision.${MOD} ${OBJDIR}/eirmod_precision.${MOD}
endif

${MNEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${LDLIBES} ${LD_CATALYST} ${LDOPTSend}

${OEEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${LDLIBES} ${LDOPTSend}

${OPEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MAKES}
	${LD} ${LDOPTS} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${LDLIBES} ${LDOPTSend}

${OQEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${LDLIBES} ${LD_CATALYST} ${LDOPTSend}

${OTEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MAKES}
	${LD} ${LDOPTS} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${LDLIBES} ${LDOPTSend}

${GEEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${GRLIBES} ${LDLIBES} ${LDOPTSend}

${GREXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MAKES}
	${LD} ${LDOPTS} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${GRLIBES} ${LDLIBES} ${LDOPTSend}

${XDEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${OBJDIR}/libsolps4.a ${MAKES}
	${LD} ${LDOPTS} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${OBJDIR}/libsolps4.a ${LCPP} ${GRLIBES} ${LDLIBES} ${LDEXTRA} ${LDOPTSend}

${MDEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${LDLIBES} ${LD_MDSPLUS} ${LDOPTSend}

${IDEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${IMASLIBS} ${LDLIBES} ${LD_CATALYST} ${LDOPTSend}

${OBJDIR}/libb2.a: ${LIBOBJS} ${SRCDIR}/include/git_version_B25.h ${DOCDIR}/b2cdci.F ${DOCDIR}/b2cdcn.F
	@${BLD} $@ ${LIBOBJS}

${SOLPS4OBJS}: ${OBJDIR}/%.o: ${SOLPS4}/%.F
	@- /bin/rm -f ${OBJDIR}/$*.f ${OBJDIR}/$*.o
ifeq ($(strip $(CPP)),)
	${FC} ${FCOPTS} ${FFLAGSEXTRA} ${DEFINES} ${EQUIVS} ${INCLUDE} ${SOLPS4INCLUDE} -c $<
else
ifeq ($(strip $(SED)),)
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} ${SOLPS4INCLUDE} $< ${OBJDIR}/$*.f
else
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} ${SOLPS4INCLUDE} $< | ${SED} > ${OBJDIR}/$*.f
endif
	${FC} ${FCOPTS} ${FFLAGSEXTRA} -c ${MODINCLUDE} ${INCMODS} ${OBJDEST} ${OBJDIR}/$*.f
endif
	@if [ -f $*.o ] ; then /bin/mv $*.o ${OBJDIR}/ ; fi

${EIR4OBJS}: ${OBJDIR}/%.o: ${EIR4}/%.F
	${FC} ${FCOPTS} ${DEFINES} ${EQUIVS} ${INCLUDE} ${OBJDEST} -c $?

${EIR4MODS}: ${OBJDIR}/%.${MOD}: ${EIR4}/%.F
	@- /bin/rm -f ${OBJDIR}/$*.f ${OBJDIR}/$*.o
ifeq ($(strip $(CPP)),)
	${FC} ${FCOPTS} ${FFLAGSEXTRA} ${DEFINES} ${EQUIVS} ${INCLUDE} ${SOLPS4INCLUDE} -c $<
else
ifeq ($(strip $(SED)),)
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} ${SOLPS4INCLUDE} $< ${OBJDIR}/$*.f
else
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} ${SOLPS4INCLUDE} $< | ${SED} > ${OBJDIR}/$*.f
endif
	${FC} ${FCOPTS} ${FFLAGSEXTRA} -c ${MODINCLUDE} ${INCMODS} ${SOLPS4INCLUDE} ${OBJDEST} ${OBJDIR}/$*.f
endif

${OBJDIR}/libsolps4.a: ${SOLPS4OBJS} ${EIR4MODS}
	${BLD} $@ ${SOLPS4OBJS}

${OBJDIR}/b2rw.o: ${OBJDIR}/eirdiag.${MOD}
${OBJDIR}/default.o: ${OBJDIR}/ceirsrt.${MOD}

ifneq (${MOD},o)
${OBJDIR}/adsp.${MOD}: ${OBJDIR}/cadgeo.${MOD} ${OBJDIR}/clogau.${MOD} ${OBJDIR}/comusr.${MOD} ${OBJDIR}/cpes.${MOD} ${OBJDIR}/ctrcei.${MOD} ${OBJDIR}/comprt.${MOD}
${OBJDIR}/avltree.${MOD}: ${OBJDIR}/ccona.${MOD}
${OBJDIR}/caprmc.${MOD}: ${OBJDIR}/cgrid.${MOD} ${OBJDIR}/comxs.${MOD} ${OBJDIR}/comsou.${MOD}
${OBJDIR}/ccflux.${MOD}: ${OBJDIR}/ctrig.${MOD} ${OBJDIR}/cgeom.${MOD}
${OBJDIR}/ccrm.${MOD}: ${OBJDIR}/cestim.${MOD} ${OBJDIR}/csdvi.${MOD} ${OBJDIR}/czt1.${MOD} ${OBJDIR}/photon.${MOD}
${OBJDIR}/comxs.${MOD}: ${OBJDIR}/cupd.${MOD}
${OBJDIR}/cpes.${MOD}: ${OBJDIR}/comprt.${MOD} ${OBJDIR}/ctrcei.${MOD} ${OBJDIR}/eirmod_precision.${MOD}
${OBJDIR}/cupd.${MOD}: ${OBJDIR}/comsig.${MOD}
${OBJDIR}/eirdiag.${MOD}: ${OBJDIR}/precision.${MOD} ${OBJDIR}/parmmod.${MOD} ${OBJDIR}/braeir.${MOD} ${OBJDIR}/ccoupl.${MOD} ${OBJDIR}/clgin.${MOD}
${OBJDIR}/eirgrid_lib.${MOD}: ${OBJDIR}/eirmap.${MOD}
endif
${OBJDIR}/adsp.o: ${OBJDIR}/cadgeo.o ${OBJDIR}/clogau.o ${OBJDIR}/comusr.o ${OBJDIR}/cpes.o ${OBJDIR}/ctrcei.o ${OBJDIR}/comprt.o
${OBJDIR}/avltree.o: ${OBJDIR}/ccona.o
${OBJDIR}/caprmc.o: ${OBJDIR}/cgrid.o ${OBJDIR}/comxs.o ${OBJDIR}/comsou.o
${OBJDIR}/ccflux.o: ${OBJDIR}/ctrig.o ${OBJDIR}/cgeom.o
${OBJDIR}/ccrm.o: ${OBJDIR}/cestim.o ${OBJDIR}/csdvi.o ${OBJDIR}/czt1.o ${OBJDIR}/photon.o
${OBJDIR}/comxs.o: ${OBJDIR}/cupd.o
${OBJDIR}/cpes.o: ${OBJDIR}/comprt.o ${OBJDIR}/ctrcei.o ${OBJDIR}/eirmod_precision.o
${OBJDIR}/cupd.o: ${OBJDIR}/comsig.o
${OBJDIR}/eirdiag.o: ${OBJDIR}/precision.o ${OBJDIR}/parmmod.o ${OBJDIR}/braeir.o ${OBJDIR}/ccoupl.o ${OBJDIR}/clgin.o
${OBJDIR}/eirgrid_lib.o: ${OBJDIR}/eirmap.o

# target 'clean' cleans up the directory.
clean :
	-mkdir ${OBJDIR}/.delete
	-mv -i ${OBJDIR}/*.o ${OBJDIR}/*.f ${OBJDIR}/*.f90 ${OBJDIR}/*.a ${OBJDIR}/*.exe ${SRCDIR}/include/git_version_B25.h ${OBJDIR}/LISTOBJ ${OBJDIR}/dependencies ${DOCDIR}/b2cdci.F ${DOCDIR}/b2cdcn.F ${OBJDIR}/.delete >& /dev/null
ifneq (${MOD},o)
	-mv -i ${OBJDIR}/*.${MOD} ${OBJDIR}/.delete >& /dev/null
endif
	-rm -rf ${OBJDIR}/.delete &

depend: ${OBJDIR}/LISTOBJ ${B2OBJS:.o=.F} ${B2F90OBJS:.o=.F90} ${EXELIST:.o=.F} ${EX90LIST:.o=.F90}
	@`which makedepend` -p'$${OBJDIR}/' ${DEFINES} -f- ${INCLUDE} $^ | \
	sed 's,^$${OBJDIR}/[^ ][^ ]*/,\$${OBJDIR}/,' | \
        sed 's,: ${SOLPSTOP},: $${SOLPSTOP},' > ${OBJDIR}/dependencies
	@echo '# 1' >> ${OBJDIR}/dependencies
ifneq (${MOD},o)
	@`which makedepend` -p'$${OBJDIR}/' ${DEFINES} -f- ${INCLUDE} ${MODLIST} -o.${MOD} | \
	sed 's,^$${OBJDIR}/[^ ][^ ]*/,\$${OBJDIR}/,' | \
        sed 's,: ${SOLPSTOP},: $${SOLPSTOP},' >> ${OBJDIR}/dependencies
	@echo '# 2' >> ${OBJDIR}/dependencies
endif
ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
	@egrep -aiH '^ {6,}use ' ${SRCLOCAL}/*.F | grep -v 'IGNORE' | awk '{sub("\\.F:",".o:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 3' >> ${OBJDIR}/dependencies
endif
	@egrep -aiH '^ {6,}use ' ${SRCDIR}/*/*.F | grep -v 'IGNORE' | awk '{sub("\\.F:",".o:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 4a' >> ${OBJDIR}/dependencies
	@egrep -aiH '^ {0,}use ' ${SRCDIR}/*/*.F90 | grep -v 'IGNORE' | awk '{sub("\\.F90:",".o:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 4b' >> ${OBJDIR}/dependencies
ifneq (${MOD},o)
	@egrep -aiH '^ {6,}use ' ${MODLISTF} | grep -v 'IGNORE' | awk '{sub("\\.F:",".${MOD}:",$$1);sub("\\.f:",".${MOD}:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 5a' >> ${OBJDIR}/dependencies
ifdef MODLISTF90
	@egrep -aiH '^ {0,}use ' ${MODLISTF90} | grep -v 'IGNORE' | awk '{sub("\\.F90:",".${MOD}:",$$1);sub("\\.f90:",".${MOD}:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 5b' >> ${OBJDIR}/dependencies
endif
endif

tags:
	rm -f ${SRCB2}/TAGS ; etags -o ${SRCB2}/TAGS ${TAGSLIST}

listobj: ${OBJDIR}/dependencies ${DOCDIR}/b2cdci.F ${DOCDIR}/b2cdcn.F
ifdef USE_EIRENE
	@rm -f ${OBJDIR}/LISTOBJ; touch ${OBJDIR}/LISTOBJ; l="OBJS ="; \
	for d in `echo "${FPATH}" | tr : \ `; do \
		l="$$l `(cd $$d > /dev/null; echo *.F)`"; \
	done; \
	for d in `echo "${FFPATH}" | tr : \ `; do \
		l="$$l `(cd $$d > /dev/null; echo *.F90)`"; \
	done; \
	l="$$l `(cd ${SRCEIR}/modules > /dev/null; echo *.f)`"; \
	l="$$l `(cd ${SRCEIR}/interfaces/couple_SOLPS-ITER > /dev/null; echo eirmod_*.F90 eirmod_*.f)`"; \
	E="-e 's/\.F90/\.o/g' -e 's/\.F/\.o/g' -e 's/\.f/\.o/g'" ; for f in ${EXCLUDELIST}; do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$l" | eval sed "$$E" > ${OBJDIR}/LISTOBJ
else
	@rm -f ${OBJDIR}/LISTOBJ; touch ${OBJDIR}/LISTOBJ; l="OBJS ="; \
	for d in `echo "${FPATH}" | tr : \ `; do \
		l="$$l `(cd $$d > /dev/null; echo *.F)`"; \
	done; \
	for d in `echo "${FFPATH}" | tr : \ `; do \
		l="$$l `(cd $$d > /dev/null; echo *.F90)`"; \
	done; \
	E="-e 's/\.F90/\.o/g' -e 's/\.F/\.o/g'" ; for f in ${EXCLUDELIST}; do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$l" | eval sed "$$E" > ${OBJDIR}/LISTOBJ
endif
	@l="B2OBJS ="; \
	for d in `echo "${FPATH}" | tr : \ `; do \
		l="$$l `(cd $$d > /dev/null; echo *.F)`"; \
	done; \
	E="-e 's/\.F90/\.o/g' -e 's/\.F/\.o/g'" ; for f in ${EXCLUDELIST}; do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$l" | eval sed "$$E" >> ${OBJDIR}/LISTOBJ
	@l="B2F90OBJS ="; \
	for d in `echo "${FFPATH}" | tr : \ `; do \
		l="$$l `(cd $$d > /dev/null; echo *.F90)`"; \
	done; \
	E="-e 's/\.F90/\.o/g'" ; for f in ${EXCLUDELIST}; do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$l" | eval sed "$$E" >> ${OBJDIR}/LISTOBJ

${OBJDIR}/LISTOBJ: listobj

VERSION: ${SRCDIR}/include/git_version_B25.h

${SRCDIR}/include/git_version_B25.h: force
	@echo "      character*32 :: git_version_B25 = '`git describe --dirty --always`'" > ${SRCDIR}/include/git_version_new.h
	@if cmp -s ${SRCDIR}/include/git_version_new.h ${SRCDIR}/include/git_version_B25.h; then rm ${SRCDIR}/include/git_version_new.h; else mv ${SRCDIR}/include/git_version_new.h ${SRCDIR}/include/git_version_B25.h; fi

${OBJDIR}/dependencies: ${SRCDIR}/modules/.new_modules
ifeq ($(shell [ -d ${OBJDIR} ] && echo yes || echo no ),no)
	-mkdir -p ${OBJDIR}
endif
	touch ${OBJDIR}/dependencies
	${MAKE} tags
	${MAKE} VERSION
	${MAKE} local
	${MAKE} listobj
	${MAKE} depend

include ${OBJDIR}/dependencies

${OBJDIR}/process.o : process.F
	@- /bin/rm -f ${OBJDIR}/process.f ${OBJDIR}/process.o
ifeq ($(strip $(CPP)),)
	${FC} ${FCOPTS} ${DEFINES} ${EQUIVS} ${INCLUDE} -c $<
else
ifeq ($(strip $(SED)),)
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< ${OBJDIR}/process.f
else
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< | ${DBLSED} > ${OBJDIR}/process.f
endif
	${DBLFC} ${DBLOPTION} ${FCOPTS} -c ${INCMOD}${OBJDIR} ${OBJDIR}/process.f
endif
	@if [ -f process.o ] ; then /bin/mv process.o ${OBJDIR}/ ; fi

ifeq ($(COMPILER),g77)
${OBJDIR}/b2stbc.o : b2stbc.F
	@- /bin/rm -f ${OBJDIR}/b2stbc.f ${OBJDIR}/b2stbc.o
ifeq ($(strip $(CPP)),)
	${FC} ${FCOPTS} ${DEFINES} ${EQUIVS} ${INCLUDE} -c $<
else
ifeq ($(strip $(SED)),)
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< ${OBJDIR}/b2stbc.f
else
	-${CPP} ${DEFINES} ${EQUIVS} -P ${INCLUDE} $< | ${SED} > ${OBJDIR}/b2stbc.f
endif
	${FC} ${FCOPTS} -O1 -c ${INCMOD}${OBJDIR} ${OBJDIR}/b2stbc.f
endif
	@if [ -f b2stbc.o ] ; then /bin/mv b2stbc.o ${OBJDIR}/ ; fi
endif

echo:
	@echo LD_CATALYST=${LD_CATALYST}
	@echo INCLUDE=${INCLUDE}
	@echo LDLIBES=${LDLIBES}
	@echo DEFINES=${DEFINES}
	@echo EQUIVS=${EQUIVS}
	@echo VPATH=${VPATH}
	@echo FPATH=${FPATH}
	@echo OBJDIR=${OBJDIR}
	@echo OBJS=${OBJS}
	@echo SOLPS4OBJS=${SOLPS4OBJS}
	@echo EIR4MODS=${EIR4MODS}
	@echo EIR4OBJS=${EIR4OBJS}
	@echo MODLIST=${MODLIST}
	@echo MODLISTF=${MODLISTF}
	@echo MODLISTF90=${MODLISTF90}
	@echo MODULES=${MODULES}
	@echo MODOBJS=${MODOBJS}
	@echo MODMODS=${MODMODS}
	@echo MODLOCAL=${MODLOCAL}
	@echo $(filter-out ${MODOBJS},${ALLOBJS})
	@echo EXCLUDELIST=${EXCLUDELIST}
	@echo EXELIST=${EXELIST}
	@echo EX90LIST=${EX90LIST}
	@echo GREXE=${GREXE}
	@echo MNEXE=${MNEXE}
	@echo XDEXE=${XDEXE}
	@echo OTEXE=${OTEXE}
	@echo O9EXE=${O9EXE}
	@echo IDEXE=${IDEXE}

local: ${SRCLOCAL}/b2local.F ${MODLOCAL}/b2mod_local.F ${INCLOCAL}/b2local.h

${SRCLOCAL}/b2local.F:
	mkdir -p ${SRCLOCAL}
	echo "      subroutine b2local" > ${SRCLOCAL}/b2local.F
	echo "c" >> ${SRCLOCAL}/b2local.F
	echo "c store local or locally modified subroutines in this directory" >> ${SRCLOCAL}/b2local.F
	echo "c" >> ${SRCLOCAL}/b2local.F
	echo "      use b2mod_local" >> ${SRCLOCAL}/b2local.F
	echo '#include "b2local.h"' >> ${SRCLOCAL}/b2local.F
	echo "c" >> ${SRCLOCAL}/b2local.F
	echo "      end" >> ${SRCLOCAL}/b2local.F
ifdef SOLPS_CPP
ifeq ($(shell [ -d ../solps4-5/src/local_5 ] && echo yes || echo no ),no)
	ln -sf ${SRCLOCAL} ../solps4-5/src/local_5
endif
endif

${MODLOCAL}/b2mod_local.F:
	mkdir -p ${MODLOCAL}
	echo "      module b2mod_local" > ${MODLOCAL}/b2mod_local.F
	echo "c" >> ${MODLOCAL}/b2mod_local.F
	echo "c store local or locally modified modules in this directory" >> ${MODLOCAL}/b2mod_local.F
	echo "c" >> ${MODLOCAL}/b2mod_local.F
	echo "      end" >> ${MODLOCAL}/b2mod_local.F
ifdef SOLPS_CPP
ifeq ($(shell [ -d ../solps4-5/src/B2.5_modules.local ] && echo yes || echo no ),no)
	ln -sf ${MODLOCAL} ../solps4-5/src/B2.5_modules.local
endif
endif

${INCLOCAL}/b2local.h:
	mkdir -p ${INCLOCAL}
	echo "c" > ${INCLOCAL}/b2local.h
	echo "c store local or locally modified include files in this directory" >> ${INCLOCAL}/b2local.h
	echo "c" >> ${INCLOCAL}/b2local.h
ifdef SOLPS_CPP
ifeq ($(shell [ -d ../solps4-5/src/B2.5_include.local ] && echo yes || echo no ),no)
	ln -sf ${INCLOCAL} ../solps4-5/src/B2.5_include.local
endif
endif

