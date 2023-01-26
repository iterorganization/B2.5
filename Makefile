# Test whether necessary environment variables are defined; if not, exit
ifndef HOST_NAME
  $(error HOST_NAME not defined)
endif
ifndef COMPILER
  $(error COMPILER not defined)
endif

SRCB2    = ${PWD}
SRCDIR   = ${SRCB2}/src
DOCDIR   = ${SRCDIR}/documentation
PYTHON  ?= python
TAGSLIST =
SOLPSINCLUDE ?=

MAKES = ${SRCB2}/Makefile
DEFINES = ${B25_DEFINES} ${SOLPS_CPP}
# Include global SOLPS compiler settings
ifndef SOLPS_CPP
  NODENAME = $(shell echo `hostname`)
  ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER} ] && echo yes || echo no ),yes)
    include ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}
    MAKES += ${SOLPSTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER} ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}
  else
    $(warning ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER} not found.)
  endif
  ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/config.common.${COMPILER} ] && echo yes || echo no ),yes)
    include ${SOLPSTOP}/SETUP/config.common.${COMPILER}
    MAKES += ${SOLPSTOP}/SETUP/config.common.${COMPILER}
  endif
  ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local ] && echo yes || echo no ),yes)
    include ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local
    MAKES += ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local
  endif
else
  MAKES += ${SOLPSTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER} ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}
  ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/config.common.${COMPILER} ] && echo yes || echo no ),yes)
    MAKES += ${SOLPSTOP}/SETUP/config.common.${COMPILER}
  endif
  ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local ] && echo yes || echo no ),yes)
    MAKES += ${SOLPSTOP}/SETUP/config.${HOST_NAME}.${COMPILER}.local
  endif
endif
ifeq ($(shell [ -e ${SOLPSTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER}.local ] && echo yes || echo no ),yes)
  MAKES += ${SOLPSTOP}/SETUP/setup.csh.${HOST_NAME}.${COMPILER}.local
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

# Extensions for object directories when various options are used
ifdef USE_IMPGYRO
EXT_IMPGYRO = .ig
else
ifdef USE_MPI
EXT_MPI = .mpi
endif
endif
ifdef USE_OPENMP
EXT_OPENMP = .openmp
endif
ifdef SOLPS_DEBUG
EXT_DEBUG = .debug
IMAS_AMNS_DEBUG = yes
else
IMAS_AMNS_DEBUG = no
endif

# Directory where objectcode/binaries will be created
OBJDIR = ${SRCB2}/builds/${PREF_OBJDIR}.${HOST_NAME}.${COMPILER}${EXT_OPENMP}${EXT_MPI}${EXT_IMPGYRO}${EXT_DEBUG}

# If compiling with Eirene, look in default place for Eirene sources/lib
ifdef USE_EIRENE
  SRCEIR = ${SOLPSTOP}/modules/Eirene/src
  EIRDIR = ${SOLPSTOP}/modules/Eirene/builds/couple_SOLPS-ITER.${HOST_NAME}.${COMPILER}${EXT_MPI}${EXT_IMPGYRO}${EXT_DEBUG}
endif
ifdef SOLPSTOP
  NCSDIR = ${SOLPSTOP}/scripts/nc2text_simple
  NCODIR = ${SOLPSTOP}/scripts/${HOST_NAME}.${COMPILER}${EXT_DEBUG}
  NCXDIR = ${SOLPSTOP}/scripts/${HOST_NAME}.${COMPILER}
endif
ifdef USE_MPI
  include ${OBJDIR}/mpiversion.mk # defines MPI_VERSION, which is the MPI version number
endif

ifeq ($(shell [ -e ${OBJDIR}/LISTOBJ ] && echo yes || echo no ),yes)
  include ${OBJDIR}/LISTOBJ
endif
include ${SRCB2}/config/compile
MAKES += ${SRCB2}/config/compile ${SRCB2}/config/config.${HOST_NAME}.${COMPILER}
ifeq ($(shell [ -e ${SRCB2}/config/config.${HOST_NAME}.${COMPILER}.local ] && echo yes || echo no ),yes)
  include ${SRCB2}/config/config.${HOST_NAME}.${COMPILER}.local
  MAKES += ${SRCB2}/config/config.${HOST_NAME}.${COMPILER}.local
endif

# Add external includes first
ifdef NCDIR
SOLPSINCLUDE += -I${NCDIR}/include
endif

ifdef USE_IMPGYRO
SOLPSINCLUDE += ${MPI_CPP}
else
ifdef USE_MPI
SOLPSINCLUDE += ${MPI_CPP}
else
SOLPSINCLUDE += -I${SRCDIR}/mpi_dummy
endif
endif

ifdef MDSPLUS_DIR
ifndef LD_MDSPLUS
LD_MDSPLUS ?= -L${MDSPLUS_DIR}/lib -lMdsLib_client
endif
SOLPSINCLUDE += -I${MDSPLUS_DIR}/include
ifndef SOLPS_CPP
DEFINES += -DMDSPLUS
endif
endif

ifdef IMAS_PREFIX
ifeq ($(shell test ${IMAS_MINOR_VERSION} -ge 12; echo $$?),0)
ifeq ($(shell test ${GGD_MAJOR_VERSION} -eq 0; echo $$?),0)
$(warning Asking for an IMAS build but missing a GGD module: build may be incomplete)
endif
endif
endif

# If compiling with Paraview Catalyst
ifdef LD_CATALYST
SRCCAT = ${SRCDIR}/catalyst
PARAVIEW_MAJOR_VERSION ?= `paraview --version | cut -d '.' -f 1 | cut -d ' ' -f 3`
CMAKE_MAJOR_VERSION ?= `cmake --version | head -1 | cut -d ' ' -f 3 | cut -d '.' -f 1`
ifeq ($(shell test ${CMAKE_MAJOR_VERSION} -gt 2; echo $$?),0)
PARAVIEW_INCLUDE ?= $(shell paraview-config --cppflags)
else
PARAVIEW_INCLUDE ?= $(shell paraview-config --include)
endif
else
PARAVIEW_INCLUDE =
endif

# Local includes
BASEDIR = ${OBJDIR}
INCLOCAL = ${SRCDIR}/include.local
MODLOCAL = ${SRCDIR}/modules.local
SRCLOCAL = ${SRCB2}/src.local
SOLPS4 = ${SOLPSTOP}/modules/solps4-5/src/solps4_B2
EIR4 = ${SOLPSTOP}/modules/solps4-5/src/Eirene_modules
ifeq ($(shell [ -d ${MODLOCAL} ] && echo yes || echo no ),yes)
SOLPSINCLUDE += -I${MODLOCAL}
TAGSLIST += ${MODLOCAL}/*.F
endif
ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
SOLPSINCLUDE += -I${SRCLOCAL}
TAGSLIST += ${SRCLOCAL}/*.F
endif
ifeq ($(shell [ -d ${INCLOCAL} ] && echo yes || echo no ),yes)
SOLPSINCLUDE += -I${INCLOCAL}
TAGSLIST += ${INCLOCAL}/*.*
endif
SOLPSINCLUDE += -I${SRCDIR}/common -I${SRCDIR}/include
SOLPS4INCLUDE = -I${SOLPSTOP}/modules/solps4-5/src/B2_include
TAGSLIST += ${SRCDIR}/include/*.* ${SRCDIR}/common/*.* ${SRCDIR}/common/COUPLE/*.F ${SRCDIR}/*/*.F ${SRCDIR}/*/*.F90 ${DOCDIR}/*.xml ${DOCDIR}/*.py

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
# OpenMP parallelization, in order to do that, just compile without OpenMP
# compiler options (ifort -qopenmp or similar)

ifdef USE_OPENMP
#DEFINES += -DNO_OPENMP_B2XPFE
#DEFINES += -DNO_OPENMP_B2SIFRTF
#DEFINES += -DNO_OPENMP_B2SIHS
#DEFINES += -DNO_OPENMP_B2SPCX
#DEFINES += -DNO_OPENMP_B2SPEL
#DEFINES += -DNO_OPENMP_B2SQCX
#DEFINES += -DNO_OPENMP_B2SQEL
#DEFINES += -DNO_OPENMP_B2SRDT
#DEFINES += -DNO_OPENMP_B2SRST
#DEFINES += -DNO_OPENMP_B2STCX
#DEFINES += -DNO_OPENMP_B2STEL
#DEFINES += -DNO_OPENMP_MYBLAS
#DEFINES += -DNO_OPENMP_SFILL
#DEFINES += -DNO_OPENMP_B2NEWS_LOOP1
#DEFINES += -DNO_OPENMP_B2NEWS_LOOP2
#DEFINES += -DNO_OPENMP_B2NEWS_LOOP3
#DEFINES += -DNO_OPENMP_B2NEWS_LOOP4
#DEFINES += -DNO_OPENMP_B2NEWS_UNDERSCORE_LOOP1
#DEFINES += -DNO_OPENMP_B2NEWS_UNDERSCORE_LOOP2
#DEFINES += -DNO_OPENMP_B2NEWS_UNDERSCORE_LOOP3
#DEFINES += -DNO_OPENMP_B2NEWS_UNDERSCORE_LOOP4
#DEFINES += -DNO_OPENMP_B2NPMO
endif

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
FPATH := ${VPATH}
VPATH += :${SRCDIR}/test
VPATH += :${SRCDIR}/ids
VPATH += :${SRCDIR}/ids/archive
FFPATH += :${SRCDIR}/ids
FFPATH += :${SRCDIR}/modules
FFPATH += :${SRCDIR}/user

MODLIST =
MODLISTF =
MODLISTF90 =
ifeq ($(shell [ -d ${MODLOCAL} ] && echo yes || echo no ),yes)
MODLIST += ${MODLOCAL}/*.F
MODLISTF += ${MODLOCAL}/*.F
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
PROG_GR = b2yg.exe b2yi.exe b2ym.exe b2yn.exe b2yp.exe b2yq.exe b2yr.exe
PROG_MN = b2mn.exe b2mnastra.exe
PROG_AM = b2ar.exe
PROG_XD = b2xd.exe
PROG_OE = b2ag.exe b2fu.exe b2ts.exe b2uf.exe b2ye.exe b2yt.exe b2ymb.exe b2yrp.exe b2ydm.exe b2plasmastate_inspect.exe calc_atomic_data.exe
PROG_CO = b2co.exe
PROG_OT = b2ah.exe b2ai.exe b2yi_gnuplot.exe b2yh.exe b2yv.exe b2fgmtry_mod.exe
PROG_90 = check_b2_output.exe
PROG_OP = b2op.exe
PROG_OQ = b2mn_opt.exe
PROG_MD = b2md.exe b2rd.exe
PROG_ID = b2_ual_write.exe b2_ual_rewrite.exe b2_ual_write_b2mod.exe
PROG_TT = test_shrink_label.exe
PROG_NC = nc2text_simple.exe
PROG_NR = nc_reduce.exe

EXCLUDELIST = ${patsubst %.exe, %\\.o, ${PROG_GE} ${PROG_GR} ${PROG_MN} ${PROG_AM} ${PROG_XD} ${PROG_OE} ${PROG_CO} ${PROG_OT} ${PROG_90} ${PROG_MD} ${PROG_OP} ${PROG_OQ} ${PROG_ID} ${PROG_TT}}
EXELIST = ${patsubst %.exe, %.o, ${PROG_GE} ${PROG_GR} ${PROG_MN} ${PROG_AM} ${PROG_XD} ${PROG_OE} ${PROG_CO} ${PROG_OT} ${PROG_MD} ${PROG_OP} ${PROG_OQ}}
EX90LIST = ${patsubst %.exe, %.o, ${PROG_90} ${PROG_ID}}

GEEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_GE}}
GREXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_GR}}
XDEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_XD}}
MNEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_MN}}
AMEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_AM}}
OEEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_OE}}
COEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_CO}}
OTEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_OT}}
O9EXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_90}}
OPEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_OP}}
OQEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_OQ}}
MDEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_MD}}
IDEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_ID}}
TTEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_TT}}
NCEXE = ${patsubst %.exe, ${NCODIR}/%.exe, ${PROG_NC}}
NREXE = ${patsubst %.exe, ${NCODIR}/%.exe, ${PROG_NR}}

.PHONY: DEFAULT NOPLOT ALL VERSION mods clean depend listobj tags echo local force test nc2text_simple nc2text

DEFAULT: VERSION ${MNEXE} ${AMEXE} ${OEEXE} ${COEXE} ${OTEXE} ${O9EXE}
ALL: VERSION ${MNEXE} ${AMEXE} ${OEEXE} ${COEXE} ${OTEXE} ${O9EXE} ${XDEXE}
NOPLOT: VERSION ${MNEXE} ${AMEXE} ${OEEXE} ${OTEXE} ${O9EXE}
ifdef NCARG_ROOT
ifeq ($(strip ${GLI_HOME}),)
$(warning B2.5 graphical post-processing programs may not work because GLI_HOME is not defined.)
endif
DEFAULT: ${GEEXE} ${GREXE}
ALL: ${GEEXE} ${GREXE}
else
$(warning B2.5 graphical post-processing programs will not be compiled because NCARG_ROOT is not defined.)
endif
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
ALL: ${IDEXE}
ids: ${IDEXE}
NOPLOT: ${IDEXE}
endif
ifdef SOLPS_CPP
ifdef LD_NETCDF
DEFAULT: ${NCEXE} ${NREXE} nc2text
ALL: ${NCEXE} ${NREXE} nc2text
NOPLOT: ${NCEXE} ${NREXE} nc2text
endif
endif
MAIN: VERSION ${MNEXE}

ifdef USE_EIRENE
VPATH+=${SRCEIR}/modules:${SRCEIR}/interfaces/couple_SOLPS-ITER
MODLIST+=${SRCEIR}/modules/*.f ${SRCEIR}/modules/*.[fF]90 ${SRCEIR}/interfaces/couple_SOLPS-ITER/eirmod_*.f ${SRCEIR}/interfaces/couple_SOLPS-ITER/eirmod_*.F90
MODLISTF+=${SRCEIR}/modules/*.f ${SRCEIR}/interfaces/couple_SOLPS-ITER/eirmod_*.f
MODLISTF90+=${SRCEIR}/modules/*.[fF]90 ${SRCEIR}/interfaces/couple_SOLPS-ITER/eirmod_*.F90
MNEXTRA=${EIRDIR}/libeirene.a ${EIRDIR}/libgr_dummy.a ${EIRDIR}/ioflush.o
else
# MNEXTRA=${EIRDIR}/eirmod_balanced_strategy.o ${EIRDIR}/eirmod_braeir.o ${EIRDIR}/eirmod_brascl.o ${EIRDIR}/eirmod_braspoi.o ${EIRDIR}/eirmod_cadgeo.o ${EIRDIR}/eirmod_cai.o ${EIRDIR}/eirmod_calstr_buffered.o ${EIRDIR}/eirmod_ccona.o ${EIRDIR}/eirmod_ccoupl.o ${EIRDIR}/eirmod_ccrm.o ${EIRDIR}/eirmod_cestim.o ${EIRDIR}/eirmod_cfplk.o ${EIRDIR}/eirmod_cgeom.o ${EIRDIR}/eirmod_cgrid.o ${EIRDIR}/eirmod_cgrptl.o ${EIRDIR}/eirmod_cinit.o ${EIRDIR}/eirmod_clast.o ${EIRDIR}/eirmod_clgin.o ${EIRDIR}/eirmod_clogau.o ${EIRDIR}/eirmod_comnnl.o ${EIRDIR}/eirmod_comprt.o ${EIRDIR}/eirmod_comsig.o ${EIRDIR}/eirmod_comsou.o ${EIRDIR}/eirmod_comspl.o ${EIRDIR}/eirmod_comusr.o ${EIRDIR}/eirmod_comxs.o ${EIRDIR}/eirmod_coutau.o ${EIRDIR}/eirmod_cpes.o ${EIRDIR}/eirmod_cpl3d.o ${EIRDIR}/eirmod_cplmsk.o ${EIRDIR}/eirmod_cplot.o ${EIRDIR}/eirmod_cpolyg.o ${EIRDIR}/eirmod_crand.o ${EIRDIR}/eirmod_crech.o ${EIRDIR}/eirmod_cref.o ${EIRDIR}/eirmod_crefmod.o ${EIRDIR}/eirmod_csdvi.o ${EIRDIR}/eirmod_csdvi_bgk.o ${EIRDIR}/eirmod_csdvi_cop.o ${EIRDIR}/eirmod_cspei.o ${EIRDIR}/eirmod_cspez.o ${EIRDIR}/eirmod_cstep.o ${EIRDIR}/eirmod_ctetra.o ${EIRDIR}/eirmod_ctext.o ${EIRDIR}/eirmod_ctrcei.o ${EIRDIR}/eirmod_ctrig.o ${EIRDIR}/eirmod_ctsurf.o ${EIRDIR}/eirmod_cupd.o ${EIRDIR}/eirmod_cvarusr.o ${EIRDIR}/eirmod_czt1.o ${EIRDIR}/eirmod_eirbra.o ${EIRDIR}/eirmod_eirdiag.o ${EIRDIR}/eirmod_infcop.o ${EIRDIR}/eirmod_module_avltree.o ${EIRDIR}/eirmod_mpi.o ${EIRDIR}/eirmod_octree.o ${EIRDIR}/eirmod_parmmod.o ${EIRDIR}/eirmod_precision.o ${EIRDIR}/eirmod_refusr.o ${EIRDIR}/eirmod_solps.o
# EXCLUDELIST += ${patsubst ${OBJDIR}/%.o, %.o, ${MNEXTRA} }
endif
ifdef LD_CATALYST
VPATH += :${SRCDIR}/catalyst
FFPATH += :${SRCDIR}/catalyst
endif

IDSMODS = ${PROG_ID:%.exe=${OBJDIR}/%.${MOD}}
MODULES = ${patsubst %.f90,%.o,${patsubst %.F90,%.o,${patsubst %.f,%.o,${patsubst %.F,%.o,${shell basename -a ${MODLIST} } } } } }
MODOBJS = ${MODULES:%.o=${OBJDIR}/%.o}
MODMODS = $(filter-out ${IDSMODS},${MODOBJS:%.o=%.${MOD}})
SOLPS4OBJS = ${patsubst ${SOLPS4}/%.F,${OBJDIR}/%.o,${shell echo ${S4LIST} } }
EIR4MODS = ${patsubst ${EIR4}/%.F,${OBJDIR}/%.${MOD},${shell echo ${E4LIST} } }
EIR4OBJS = ${patsubst ${EIR4}/%.F,${OBJDIR}/%.o,${shell echo ${E4LIST} } }

ifeq (${MOD},o)
LIBOBJS = $(filter-out ${MODOBJS},${ALLOBJS})
else
LIBOBJS = ${ALLOBJS}
endif

mods : ${MODMODS}

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

${OBJDIR}/eirmod_wneutrals.${MOD}:
	ln -sf ${EIRDIR}/eirmod_wneutrals.${MOD} ${OBJDIR}

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

${OBJDIR}/eirmod_refusr.${MOD}:
	ln -sf ${EIRDIR}/eirmod_refusr.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_solps.${MOD}:
	ln -sf ${EIRDIR}/eirmod_solps.${MOD} ${OBJDIR}
endif

${OBJDIR}/eirmod_extrab25.o:
	ln -sf ${EIRDIR}/eirmod_extrab25.o ${OBJDIR}

${OBJDIR}/eirmod_wneutrals.o:
	ln -sf ${EIRDIR}/eirmod_wneutrals.o ${OBJDIR}

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

${OBJDIR}/eirmod_refusr.o:
	ln -sf ${EIRDIR}/eirmod_refusr.o ${OBJDIR}

${OBJDIR}/eirmod_solps.o:
	ln -sf ${EIRDIR}/eirmod_solps.o ${OBJDIR}
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

${OBJDIR}/eirmod_clgin.${MOD}:
	touch ${OBJDIR}/eirmod_clgin.${MOD}

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

${OBJDIR}/eirmod_refusr.${MOD}:
	touch ${OBJDIR}/eirmod_refusr.${MOD}

${OBJDIR}/eirmod_solps.${MOD}:
	touch ${OBJDIR}/eirmod_solps.${MOD}

${OBJDIR}/eirmod_wneutrals.${MOD}:
	touch ${OBJDIR}/eirmod_wneutrals.${MOD}
endif

ifeq ($(COMPILER),gfortran)
ifeq ($(shell test ${GFORTRAN_MAJOR_VERSION} -ge 10; echo $$?),0)
${OBJDIR}/b2mod_mdsplus.o : b2mod_mdsplus.F
	@- /bin/rm -f ${OBJDIR}/b2mod_mdsplus.f ${OBJDIR}/b2mod_mdsplus.o ${OBJDIR}/b2mod_mdsplus.${MOD}
ifeq ($(strip $(CPP)),)
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} ${DEFINES} ${EQUIVS} ${SOLPSINCLUDE} -c $<
else
ifeq ($(strip $(SED)),)
	-${CPP} ${DEFINES} ${EQUIVS} -P ${SOLPSINCLUDE} $< ${OBJDIR}/b2mod_mdsplus.f
else
	-${CPP} ${DEFINES} ${EQUIVS} -P ${SOLPSINCLUDE} $< | ${SED} > ${OBJDIR}/b2mod_mdsplus.f
endif
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${MODINCLUDE} ${INCMODS} -o ${OBJDIR}/b2mod_mdsplus.o ${OBJDIR}/b2mod_mdsplus.f
endif
	@if [ -f b2mod_mdsplus.o ] ; then /bin/mv b2mod_mdsplus.o ${OBJDIR}/ ; fi
	@if [ -f b2mod_mdsplus.${MOD} ] ; then /bin/mv b2mod_mdsplus.${MOD} ${OBJDIR}/ ; fi

ifneq (${MOD},o)
${OBJDIR}/b2mod_mdsplus.${MOD} : b2mod_mdsplus.F
	@- /bin/rm -f ${OBJDIR}/b2mod_mdsplus.f ${OBJDIR}/b2mod_mdsplus.o ${OBJDIR}/b2mod_mdsplus.${MOD}
ifeq ($(strip $(CPP)),)
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} ${DEFINES} ${EQUIVS} ${SOLPSINCLUDE} -c $<
else
ifeq ($(strip $(SED)),)
	-${CPP} ${DEFINES} ${EQUIVS} -P ${SOLPSINCLUDE} $< ${OBJDIR}/b2mod_mdsplus.f
else
	-${CPP} ${DEFINES} ${EQUIVS} -P ${SOLPSINCLUDE} $< | ${SED} > ${OBJDIR}/b2mod_mdsplus.f
endif
	${FC} ${FCOPTS} -fallow-argument-mismatch ${FFLAGSEXTRA} -c ${MODINCLUDE} ${INCMODS} -o ${OBJDIR}/b2mod_mdsplus.o ${OBJDIR}/b2mod_mdsplus.f
endif
	@if [ -f b2mod_mdsplus.o ] ; then /bin/mv b2mod_mdsplus.o ${OBJDIR}/ ; fi
ifeq ($(strip $(LINK_MOD)),)
	@if [ -f b2mod_mdsplus.${MOD} ] ; then /bin/mv b2mod_mdsplus.${MOD} ${OBJDIR}/ ; fi
else
	@[ -f b2mod_mdsplus.${MOD} ] && [ ! -L b2mod_mdsplus.${MOD} ] && /bin/mv b2mod_mdsplus.${MOD} ${OBJDIR}/ && /bin/ln -s ${OBJDIR}/b2mod_mdsplus.${MOD} .
endif
endif
endif
endif

${MNEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${IMASLIBS} ${PLLIBES} ${LDLIBES} ${LD_CATALYST} ${LDOPTSend}

${AMEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${IMASLIBS} ${LDLIBES} ${LDOPTSend}

${OEEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${LDLIBES} ${LDOPTSend}

${COEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${PLLIBES} ${LDLIBES} ${LDOPTSend}

${OPEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${LDLIBES} ${LDOPTSend}

${OQEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${LDLIBES} ${LD_CATALYST} ${LDOPTSend}

${OTEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${LDLIBES} ${LDOPTSend}

${O9EXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${LDLIBES}

${GEEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${PLLIBES} ${GRLIBES} ${LDLIBES} ${LDOPTSend}

${GREXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${GRLIBES} ${LDLIBES} ${LDOPTSend}

${XDEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${OBJDIR}/libsolps4.a ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${OBJDIR}/libsolps4.a ${GRLIBES} ${LDLIBES} ${LDEXTRA} ${LDOPTSend}

${MDEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${LDLIBES} ${LD_MDSPLUS} ${LDOPTSend}

${IDEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA} ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${MNEXTRA} ${IMASLIBS} ${PLLIBES} ${LDLIBES} ${LD_CATALYST} ${LDOPTSend}

${TTEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MAKES}
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${OBJDIR}/$*.o ${OBJDIR}/libb2.a ${LDLIBES} ${LDOPTSend}

${NCEXE}: ${NCODIR}/%.exe: ${NCODIR}/%.o ${MAKES}
ifdef LD_NETCDF
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${NCODIR}/$*.o ${LD_NETCDF}
	@-ln -sf $@ ${NCXDIR}/$*
	@-ln -sf ${NCXDIR}/nc2text_simple ${NCXDIR}/nc2text
else
	$(warning NETCDF library not present!)
endif

${NREXE}: ${NCODIR}/%.exe: ${NCODIR}/%.o ${OBJDIR}/cdf_routines.o ${OBJDIR}/chcase.o ${OBJDIR}/ifill.o ${OBJDIR}/isadigit.o ${OBJDIR}/lnblnk.o ${OBJDIR}/machsfr.o ${OBJDIR}/nagsubst.o ${OBJDIR}/open_file.o ${OBJDIR}/prgend.o ${OBJDIR}/prgini.o ${OBJDIR}/prvrt.o ${OBJDIR}/prvrti.o ${OBJDIR}/sfill.o ${OBJDIR}/streql.o ${OBJDIR}/sysend.o ${OBJDIR}/sysini.o ${OBJDIR}/xerrab.o ${OBJDIR}/xertst.o ${MAKES}
ifdef LD_NETCDF
	${LD} ${LDOPTS} ${FFLAGSEXTRA} -o $@ ${NCODIR}/$*.o ${OBJDIR}/b2mod_ipmain.o ${OBJDIR}/b2mod_lwimai.o ${OBJDIR}/b2mod_lwmain.o ${OBJDIR}/b2mod_math.o ${OBJDIR}/b2mod_openmp.o ${OBJDIR}/b2mod_subsys.o ${OBJDIR}/b2mod_xerset.o ${OBJDIR}/cdf_routines.o ${OBJDIR}/chcase.o ${OBJDIR}/ifill.o ${OBJDIR}/isadigit.o ${OBJDIR}/lnblnk.o ${OBJDIR}/machsfr.o ${OBJDIR}/nagsubst.o ${OBJDIR}/open_file.o ${OBJDIR}/prgend.o ${OBJDIR}/prgini.o ${OBJDIR}/prvrt.o ${OBJDIR}/prvrti.o ${OBJDIR}/sfill.o ${OBJDIR}/streql.o ${OBJDIR}/sysend.o ${OBJDIR}/sysini.o ${OBJDIR}/xerrab.o ${OBJDIR}/xertst.o ${LD_NETCDF} ${LD_NAG}
ifndef SOLPS_DEBUG
	@-ln -sf $@ ${NCXDIR}/$*
endif
else
	$(warning NETCDF library not present!)
endif

${OBJDIR}/libb2.a: ${LIBOBJS} ${SRCDIR}/include/git_version_B25.h ${DOCDIR}/b2cdci.F ${DOCDIR}/b2cdcn.F
	@${BLD} $@ ${LIBOBJS}

test:	${TTEXE}

nc2text: ${NCEXE}

nc2text_simple: ${NCEXE}

nc_reduce: ${NREXE}

${NCODIR}/nc2text_simple.o: ${NCSDIR}/nc2text_simple.F90
ifdef LD_NETCDF
	@-mkdir -p ${NCODIR}
	-${CPP} ${DEFINES} ${EQUIVS} -P ${SOLPSINCLUDE} $< $*.F90
	${FC} ${FCOPTS} ${FFLAGSEXTRA} -c -o $*.o $*.F90
else
	$(warning NETCDF library not present!)
endif

${NCODIR}/nc_reduce.o: ${NCSDIR}/nc_reduce.F90 ${OBJDIR}/b2mod_math.o ${OBJDIR}/b2mod_subsys.o ${OBJDIR}/b2mod_types.o
ifdef LD_NETCDF
	@-mkdir -p ${NCODIR}
	-${CPP} ${DEFINES} ${EQUIVS} -P ${SOLPSINCLUDE} $< $*.F90
	${FC} ${FCOPTS} ${FFLAGSEXTRA} -c ${MODINCLUDE} ${INCMODS} -o $*.o $*.F90
else
	$(warning NETCDF library not present!)
endif

${SOLPS4OBJS}: ${OBJDIR}/%.o: ${SOLPS4}/%.F
	@- /bin/rm -f ${OBJDIR}/$*.f ${OBJDIR}/$*.o
ifeq ($(strip $(CPP)),)
	${FC} ${FCOPTS} ${FFLAGSEXTRA} ${DEFINES} ${EQUIVS} ${SOLPSINCLUDE} ${SOLPS4INCLUDE} -c $<
else
ifeq ($(strip $(SED)),)
	-${CPP} ${DEFINES} -DSOLPS_ITER ${EQUIVS} -P ${SOLPSINCLUDE} ${SOLPS4INCLUDE} $< ${OBJDIR}/$*.f
else
	-${CPP} ${DEFINES} -DSOLPS_ITER ${EQUIVS} -P ${SOLPSINCLUDE} ${SOLPS4INCLUDE} $< | ${SED} > ${OBJDIR}/$*.f
endif
	${FC} ${FCOPTS} ${FFLAGSEXTRA} -c ${MODINCLUDE} ${INCMODS} ${OBJDEST} ${OBJDIR}/$*.f
endif
	@if [ -f $*.o ] ; then /bin/mv $*.o ${OBJDIR}/ ; fi

${EIR4OBJS}: ${OBJDIR}/%.o: ${EIR4}/%.F
	${FC} ${FCOPTS} ${FFLAGSEXTRA} ${DEFINES} ${EQUIVS} ${SOLPSINCLUDE} ${OBJDEST} -c $?

${EIR4MODS}: ${OBJDIR}/%.${MOD}: ${EIR4}/%.F
	@- /bin/rm -f ${OBJDIR}/$*.f ${OBJDIR}/$*.o
ifeq ($(strip $(CPP)),)
	${FC} ${FCOPTS} ${FFLAGSEXTRA} ${DEFINES} ${EQUIVS} ${SOLPSINCLUDE} ${SOLPS4INCLUDE} -c $<
else
ifeq ($(strip $(SED)),)
	-${CPP} ${DEFINES} ${EQUIVS} -P ${SOLPSINCLUDE} ${SOLPS4INCLUDE} $< ${OBJDIR}/$*.f
else
	-${CPP} ${DEFINES} ${EQUIVS} -P ${SOLPSINCLUDE} ${SOLPS4INCLUDE} $< | ${SED} > ${OBJDIR}/$*.f
endif
	${FC} ${FCOPTS} ${FFLAGSEXTRA} -c ${MODINCLUDE} ${INCMODS} ${SOLPS4INCLUDE} ${OBJDEST} ${OBJDIR}/$*.f
	@touch ${OBJDIR}/$*.${MOD}
endif

ifdef USE_EIRENE
${OBJDIR}/libsolps4.a: ${SOLPS4OBJS} ${EIR4MODS}
else
${OBJDIR}/libsolps4.a: ${SOLPS4OBJS}
endif
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
	-mv -i ${OBJDIR}/*.o ${OBJDIR}/*.f ${OBJDIR}/*.f90 ${OBJDIR}/*.a ${OBJDIR}/*.exe ${SRCDIR}/include/git_version_B25.h ${OBJDIR}/LISTOBJ ${OBJDIR}/dependencies ${OBJDIR}/mpiversion.mk ${OBJDIR}/.delete >& /dev/null
ifneq (${MOD},o)
	-mv -i ${OBJDIR}/*.${MOD} ${OBJDIR}/.delete >& /dev/null
endif
ifdef SOLPSTOP
ifdef LD_NETCDF
	-rm -f ${NCODIR}/*.f90 ${NCODIR}/*.o ${NCODIR}/*.exe
endif
endif
	-rm -rf ${OBJDIR}/.delete &

depend: ${OBJDIR}/LISTOBJ ${B2OBJS:.o=.F} ${B2F90OBJS:.o=.F90} ${EXELIST:.o=.F} ${EX90LIST:.o=.F90}
	@`which makedepend` -p'$${OBJDIR}/' ${DEFINES} -f- ${SOLPSINCLUDE} $^ | \
	sed 's,^$${OBJDIR}/[^ ][^ ]*/,\$${OBJDIR}/,' | \
        sed 's,: ${SOLPSTOP},: $${SOLPSTOP},' > ${OBJDIR}/dependencies
	@echo '# 1' >> ${OBJDIR}/dependencies
ifneq (${MOD},o)
	@`which makedepend` -p'$${OBJDIR}/' ${DEFINES} -f- ${SOLPSINCLUDE} ${MODLIST} -o.${MOD} | \
	sed 's,^$${OBJDIR}/[^ ][^ ]*/,\$${OBJDIR}/,' | \
        sed 's,: ${SOLPSTOP},: $${SOLPSTOP},' >> ${OBJDIR}/dependencies
	@echo '# 2' >> ${OBJDIR}/dependencies
endif
ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
	@egrep -aiH '^ {6,}use ' ${SRCLOCAL}/*.F | grep -v 'IGNORE' | tr , ' ' | awk '{sub("\\.F:",".o:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 3' >> ${OBJDIR}/dependencies
endif
	@egrep -aiH '^ {6,}use ' ${SRCDIR}/*/*.F | grep -v 'IGNORE' | tr , ' ' | awk '{sub("\\.F:",".o:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 4a' >> ${OBJDIR}/dependencies
	@egrep -aiH '^ {0,}use ' ${SRCDIR}/*/*.F90 | grep -v 'IGNORE' | tr , ' ' | awk '{sub("\\.F90:",".o:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 4b' >> ${OBJDIR}/dependencies
ifneq (${MOD},o)
	@egrep -aiH '^ {6,}use ' ${MODLISTF} | grep -v 'IGNORE' | tr , ' ' | awk '{sub("\\.F:",".${MOD}:",$$1);sub("\\.f:",".${MOD}:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 5a' >> ${OBJDIR}/dependencies
ifdef MODLISTF90
	@egrep -aiH '^ {0,}use ' ${MODLISTF90} | grep -v 'IGNORE' | tr , ' ' | awk '{sub("\\.F90:",".${MOD}:",$$1);sub("\\.f90:",".${MOD}:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"tolower($$3)".${MOD}"}' >> ${OBJDIR}/dependencies
	@echo '# 5b' >> ${OBJDIR}/dependencies
endif
endif

tags:
	rm -f ${SRCB2}/TAGS ; ctags -e -f ${SRCB2}/TAGS ${TAGSLIST} || touch ${SRCB2}/TAGS

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
	E="-e 's/ \*\.F90//g' -e 's/ \*\.F//g' -e 's/ eirmod_\*\.F90//g' -e 's/eirmod_\*\.f//g' -e 's/\.F90/\.o/g' -e 's/\.F/\.o/g' -e 's/\.f/\.o/g'" ; for f in ${EXCLUDELIST}; do \
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
	E="-e 's/ \*\.F90//g' -e 's/ \*\.F//g' -e 's/\.F90/\.o/g' -e 's/\.F/\.o/g'" ; for f in ${EXCLUDELIST}; do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$l" | eval sed "$$E" > ${OBJDIR}/LISTOBJ
endif
	@ll="B2OBJS ="; \
	for d in `echo "${FPATH}" | tr : \ `; do \
		ll="$$ll `(cd $$d > /dev/null; echo *.F)`"; \
	done; \
	E="-e 's/ \*\.F//g' -e 's/\.F/\.o/g'" ; for f in ${EXCLUDELIST}; do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$ll" | eval sed "$$E" >> ${OBJDIR}/LISTOBJ
	@lll="B2F90OBJS ="; \
	for d in `echo "${FFPATH}" | tr : \ `; do \
		lll="$$lll `(cd $$d > /dev/null; echo *.F90)`"; \
	done; \
	E="-e 's/ \*\.F90//g' -e 's/\.F90/\.o/g'" ; for f in ${EXCLUDELIST}; do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$lll" | eval sed "$$E" >> ${OBJDIR}/LISTOBJ

${OBJDIR}/LISTOBJ: listobj

VERSION: ${SRCDIR}/include/git_version_B25.h

${SRCDIR}/include/git_version_B25.h: force
	@echo "      character*32 :: git_version_B25 =" > ${SRCDIR}/include/git_version_new.h
	@echo "     . '`git describe --tags --dirty --always | cut -c 1-32`'" >> ${SRCDIR}/include/git_version_new.h
ifdef SOLPS_CPP
	@echo "      character*32 :: git_version_ADAS =" >> ${SRCDIR}/include/git_version_new.h
	@echo "     . '`( cd $${SOLPSTOP}/modules/adas ; git describe --tags --dirty --always | cut -c 1-32 )`'" >> ${SRCDIR}/include/git_version_new.h
	@echo "      character*32 :: git_version_SOLPS =" >> ${SRCDIR}/include/git_version_new.h
	@echo "     . '`( cd $${SOLPSTOP} ; git describe --tags --dirty --always | cut -c 1-32 )`'" >> ${SRCDIR}/include/git_version_new.h
else
	@echo "      character*32 :: git_version_ADAS = '0.0.0-0-g0000000'" >> ${SRCDIR}/include/git_version_new.h
	@echo "      character*32 :: git_version_SOLPS = '0.0.0-0-g0000000'" >> ${SRCDIR}/include/git_version_new.h
endif
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

ifeq ($(COMPILER),g77)
${OBJDIR}/b2stbc.o : b2stbc.F
	@- /bin/rm -f ${OBJDIR}/b2stbc.f ${OBJDIR}/b2stbc.o
ifeq ($(strip $(CPP)),)
	${FC} ${FCOPTS} ${FFLAGSEXTRA} ${DEFINES} ${EQUIVS} ${SOLPSINCLUDE} -c $<
else
ifeq ($(strip $(SED)),)
	-${CPP} ${DEFINES} ${EQUIVS} -P ${SOLPSINCLUDE} $< ${OBJDIR}/b2stbc.f
else
	-${CPP} ${DEFINES} ${EQUIVS} -P ${SOLPSINCLUDE} $< | ${SED} > ${OBJDIR}/b2stbc.f
endif
	${FC} ${FCOPTS} ${FFLAGSEXTRA} -O1 -c ${INCMOD}${OBJDIR} ${OBJDIR}/b2stbc.f
endif
	@if [ -f b2stbc.o ] ; then /bin/mv b2stbc.o ${OBJDIR}/ ; fi
endif

echo:
	@echo LD_CATALYST=${LD_CATALYST}
	@echo SOLPSINCLUDE=${SOLPSINCLUDE}
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
	@echo AMEXE=${AMEXE}
	@echo XDEXE=${XDEXE}
	@echo OTEXE=${OTEXE}
	@echo O9EXE=${O9EXE}
	@echo IDEXE=${IDEXE}
	@echo NCEXE=${NCEXE}
	@echo NREXE=${NREXE}

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

${MODLOCAL}/b2mod_local.F:
	mkdir -p ${MODLOCAL}
	echo "      module b2mod_local" > ${MODLOCAL}/b2mod_local.F
	echo "c" >> ${MODLOCAL}/b2mod_local.F
	echo "c store local or locally modified modules in this directory" >> ${MODLOCAL}/b2mod_local.F
	echo "c" >> ${MODLOCAL}/b2mod_local.F
	echo "      end" >> ${MODLOCAL}/b2mod_local.F

${INCLOCAL}/b2local.h:
	mkdir -p ${INCLOCAL}
	echo "c" > ${INCLOCAL}/b2local.h
	echo "c store local or locally modified include files in this directory" >> ${INCLOCAL}/b2local.h
	echo "c" >> ${INCLOCAL}/b2local.h

${OBJDIR}/mpiversion.mk: ${MAKES}
	printf "use mpi\nWRITE(*,fmt='(A12,I1)') 'MPI_VERSION=', MPI_VERSION\nWRITE(*,fmt='(A9)') 'MPI_MOD=1'\nEND" > ${OBJDIR}/mpi_version.f90
	(${FC} ${FCOPTS} ${SOLPSINCLUDE} -o ${OBJDIR}/mpi_version ${OBJDIR}/mpi_version.f90 ${LD_MPI} && ( ${OBJDIR}/mpi_version | tail -n2 ) || \
	( printf "include 'mpif.h'\nWRITE(*,fmt='(A12,I1)') 'MPI_VERSION=', MPI_VERSION\nWRITE(*,fmt='(A9)') 'MPI_MOD=0'\nEND" > ${OBJDIR}/mpi_version.f90 ; \
	${FC} ${FCOPTS} ${SOLPSINCLUDE} -o ${OBJDIR}/mpi_version ${OBJDIR}/mpi_version.f90 ${LD_MPI} && ( ${OBJDIR}/mpi_version | tail -n2 ) || ( echo MPI_VERSION=0 ; echo MPI_MOD=0 ) ) ) > ${OBJDIR}/mpiversion.mk
