SRCB2 = ${PWD}
SRCDIR = ${SRCB2}/src
ifeq ($(shell [ -d ${SOLPSTOP} ] && echo yes || echo no ),yes)
ifdef USE_MPI
OBJDIR=${SOLPSTOP}/bin/${OBJECTCODE}/B2.5.mpi
else
OBJDIR=${SOLPSTOP}/bin/${OBJECTCODE}/B2.5
endif
ifdef USE_EIRENE
SRCEIR = ${SOLPSTOP}/src/Eirene
ifdef USE_MPI
EIRDIR = ${SOLPSTOP}/bin/${OBJECTCODE}/B25eirene.mpi/Eirene
else
EIRDIR = ${SOLPSTOP}/bin/${OBJECTCODE}/B25eirene/Eirene
endif
endif
else
ifdef USE_MPI
OBJDIR=bin/${OBJECTCODE}.mpi
else
OBJDIR=bin/${OBJECTCODE}
endif
endif

BASEDIR = ${OBJDIR}
SRCLOCAL = ${SRCB2}/src.local
ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
INCLUDE = -I${SRCLOCAL}
else
INCLUDE =
endif
INCLUDE += -I${SRCDIR}/common -I${SRCDIR}/include.local -I${SRCDIR}/include

DEFINES = -DWANT_THIS ${SOLPS_CPP}
ifdef USE_MPI
DEFINES += ${USE_MPI}
endif
ifdef PERFMON
DEFINES += ${PERFMON}
endif
ifdef USE_EIRENE
DEFINES += ${USE_EIRENE}
endif

ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
VPATH+=${SRCDIR}/modules.local:${SRCLOCAL}:${SRCDIR}/modules:${SRCDIR}/b2aux:${SRCDIR}/convert:${SRCDIR}/documentation:${SRCDIR}/driver:${SRCDIR}/equations:${SRCDIR}/input:${SRCDIR}/output:${SRCDIR}/postprocessing:${SRCDIR}/preprocessing:${SRCDIR}/solvers:${SRCDIR}/sources:${SRCDIR}/transport:${SRCDIR}/utility:${SRCDIR}/b2plot:${SRCDIR}/user
else
VPATH+=${SRCDIR}/modules:${SRCDIR}/b2aux:${SRCDIR}/convert:${SRCDIR}/documentation:${SRCDIR}/driver:${SRCDIR}/equations:${SRCDIR}/input:${SRCDIR}/output:${SRCDIR}/postprocessing:${SRCDIR}/preprocessing:${SRCDIR}/solvers:${SRCDIR}/sources:${SRCDIR}/transport:${SRCDIR}/utility:${SRCDIR}/b2plot:${SRCDIR}/user
endif
FPATH:=${VPATH}

ifeq ($(shell [ -e ${OBJDIR}/LISTOBJ ] && echo yes || echo no ),yes)
include ${OBJDIR}/LISTOBJ
endif
include ${SRCB2}/config/compile
MAKES = ${SRCB2}/Makefile ${SRCB2}/config/compile ${SRCB2}/config/compiler.${OBJECTCODE}
ifeq ($(shell [ -e ${SRCB2}/config.local/compiler.${OBJECTCODE} ] && echo yes || echo no ),yes)
include ${SRCB2}/config.local/compiler.${OBJECTCODE}
MAKES+ = ${SRCB2}/config.local/compiler.${OBJECTCODE}
endif

ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
MODLIST = ${SRCDIR}/modules.local/*.F ${SRCDIR}/modules/*.F
else
MODLIST = ${SRCDIR}/modules/*.F
endif
ALLOBJS = ${OBJS:%.o=${OBJDIR}/%.o}
SOLPS4OBJS = ${SOLPS4}/readwrite.o ${SOLPS4}/b2rw.o ${SOLPS4}/calcalpha.o \
	${SOLPS4}/calcalphatrigger.o ${SOLPS4}/calcdifpr.o ${SOLPS4}/calcdifpr0.o \
	${SOLPS4}/calcdifni.o ${SOLPS4}/calcdifni0.o \
	${SOLPS4}/calcdt.o ${SOLPS4}/calcdrake.o ${SOLPS4}/calcfeqp.o ${SOLPS4}/calckxe.o \
	${SOLPS4}/calckxi.o ${SOLPS4}/calckye.o ${SOLPS4}/calckyi.o ${SOLPS4}/calclnlam.o \
	${SOLPS4}/calcparvis.o ${SOLPS4}/calcthc0.o ${SOLPS4}/calctravis.o ${SOLPS4}/calcvconv.o \
	${SOLPS4}/calcvis0.o ${SOLPS4}/calczeiler.o ${SOLPS4}/userdt.o ${SOLPS4}/userelm.o \
	${SOLPS4}/userfeqp.o ${SOLPS4}/userkxe.o ${SOLPS4}/userkxi.o ${SOLPS4}/userkye.o \
	${SOLPS4}/userkyi.o ${SOLPS4}/userni.o ${SOLPS4}/userparvis.o ${SOLPS4}/userpr.o \
	${SOLPS4}/usertravis.o ${SOLPS4}/uservconv.o ${SOLPS4}/xdr_logical.o 
SRCF = ${OBJS:%.o=%.F}

PROG_GR = b2yg.exe b2yi.exe b2ym.exe b2yn.exe b2yp.exe b2yq.exe b2yr.exe b2pl.exe
ifdef USE_MPI
PROG_MN = b2mn_mpi.exe
else
PROG_MN = b2mn_nompi.exe
endif
EXCL_MN = b2mn.exe b2mn_mpi.exe b2mn_nompi.exe
PROG_XD = b2xd.exe
PROG_OT = b2ag.exe b2ah.exe b2ai.exe b2ar.exe b2co.exe b2uf.exe b2fu.exe b2ts.exe b2yi_gnuplot.exe b2yh.exe b2yt.exe b2yv.exe b2fgmtry_mod.exe
PROG_OP = b2op.exe b2mn_opt.exe
PROG_MD = b2md.exe b2rd.exe

EXCLUDELIST = ${patsubst %.exe, %.o, ${PROG_GR} ${EXCL_MN} ${PROG_XD} ${PROG_OT} ${PROG_MD} ${PROG_OP} }
EXELIST = ${patsubst %.exe, %.o, ${PROG_GR} ${PROG_MN} ${PROG_XD} ${PROG_OT} ${PROG_MD} ${PROG_OP} }

GREXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_GR}}
XDEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_XD}}
MNEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_MN}}
OTEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_OT}}
OPEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_OP}}
MDEXE = ${patsubst %.exe, ${OBJDIR}/%.exe, ${PROG_MD}}

.PHONY: DEFAULT NOPLOT ALL VERSION clean depend listobj tags force

ifdef MDSPLUS_DIR
ifdef PAR_OPT
DEFAULT: VERSION ${MNEXE} ${OTEXE} ${GREXE} ${MDEXE} ${OPEXE}
ALL: VERSION ${MNEXE} ${OTEXE} ${GREXE} ${XDEXE} ${MDEXE} ${OPEXE}
NOPLOT: VERSION ${MNEXE} ${OTEXE} ${MDEXE} ${OPEXE}
else
DEFAULT: VERSION ${MNEXE} ${OTEXE} ${GREXE} ${MDEXE}
ALL: VERSION ${MNEXE} ${OTEXE} ${GREXE} ${XDEXE} ${MDEXE}
NOPLOT: VERSION ${MNEXE} ${OTEXE} ${MDEXE}
endif
else
ifdef PAR_OPT
DEFAULT: VERSION ${MNEXE} ${OTEXE} ${GREXE} ${OPEXE}
ALL: VERSION ${MNEXE} ${OTEXE} ${GREXE} ${XDEXE} ${OPEXE}
NOPLOT: VERSION ${MNEXE} ${OTEXE} ${OPEXE}
else
DEFAULT: VERSION ${MNEXE} ${OTEXE} ${GREXE}
ALL: VERSION ${MNEXE} ${OTEXE} ${GREXE} ${XDEXE}
NOPLOT: VERSION ${MNEXE} ${OTEXE}
endif
endif
MAIN: VERSION ${MNEXE}

ifdef USE_EIRENE
VPATH=${FPATH}:${SRCEIR}/modules:${SRCEIR}/interfaces/couple_SOLPS-ITER
MODLIST+=${SRCEIR}/modules/*.f ${SRCEIR}/interfaces/couple_SOLPS-ITER/eirmod_*.F90
MNEXTRA=${EIRDIR}/libeirene.a ${EIRDIR}/libgr_dummy.a ${EIRDIR}/ioflush.o
else
# MNEXTRA=${EIRDIR}/eirmod_braeir.o ${EIRDIR}/eirmod_brascl.o ${EIRDIR}/eirmod_braspoi.o ${EIRDIR}/eirmod_cadgeo.o ${EIRDIR}/eirmod_cai.o ${EIRDIR}/eirmod_ccona.o ${EIRDIR}/eirmod_ccoupl.o ${EIRDIR}/eirmod_cestim.o ${EIRDIR}/eirmod_cfplk.o ${EIRDIR}/eirmod_cgeom.o ${EIRDIR}/eirmod_cgrid.o ${EIRDIR}/eirmod_cgrptl.o ${EIRDIR}/eirmod_cinit.o ${EIRDIR}/eirmod_clast.o ${EIRDIR}/eirmod_clgin.o ${EIRDIR}/eirmod_clogau.o ${EIRDIR}/eirmod_comnnl.o ${EIRDIR}/eirmod_comprt.o ${EIRDIR}/eirmod_comsig.o ${EIRDIR}/eirmod_comsou.o ${EIRDIR}/eirmod_comspl.o ${EIRDIR}/eirmod_comusr.o ${EIRDIR}/eirmod_comxs.o ${EIRDIR}/eirmod_coutau.o ${EIRDIR}/eirmod_cpes.o ${EIRDIR}/eirmod_cpl3d.o ${EIRDIR}/eirmod_cplmsk.o ${EIRDIR}/eirmod_cplot.o ${EIRDIR}/eirmod_cpolyg.o ${EIRDIR}/eirmod_crand.o ${EIRDIR}/eirmod_crech.o ${EIRDIR}/eirmod_cref.o ${EIRDIR}/eirmod_crefmod.o ${EIRDIR}/eirmod_csdvi.o ${EIRDIR}/eirmod_csdvi_bgk.o ${EIRDIR}/eirmod_csdvi_cop.o ${EIRDIR}/eirmod_cspei.o ${EIRDIR}/eirmod_cspez.o ${EIRDIR}/eirmod_cstep.o ${EIRDIR}/eirmod_ctetra.o ${EIRDIR}/eirmod_ctext.o ${EIRDIR}/eirmod_ctrcei.o ${EIRDIR}/eirmod_ctrig.o ${EIRDIR}/eirmod_ctsurf.o ${EIRDIR}/eirmod_cupd.o ${EIRDIR}/eirmod_cvarusr.o ${EIRDIR}/eirmod_czt1.o ${EIRDIR}/eirmod_eirbra.o ${EIRDIR}/eirmod_module_avltree.o ${EIRDIR}/eirmod_octree.o ${EIRDIR}/eirmod_parmmod.o ${EIRDIR}/eirmod_precision.o 
# EXCLUDELIST += ${patsubst ${OBJDIR}/%.o, %.o, ${MNEXTRA} }
endif

MODULES = ${patsubst %.F %.f %.F90,%.o,${shell echo ${MODLIST} } }
MODMODS = ${MODULES:%.o=${OBJDIR}/%.${MOD}}
MODOBJS = ${MODULES:%.o=${OBJDIR}/%.o}

ifeq (${MOD},o)
LIBOBJS = $(filter-out ${MODOBJS},${ALLOBJS})
else
LIBOBJS = ${ALLOBJS}
endif

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

${OBJDIR}/eirmod_ccona.${MOD}:
	ln -sf ${EIRDIR}/eirmod_ccona.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_ccoupl.${MOD}:
	ln -sf ${EIRDIR}/eirmod_ccoupl.${MOD} ${OBJDIR}

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

${OBJDIR}/eirmod_cpl3ot.${MOD}:
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

${OBJDIR}/eirmod_module_avltree.${MOD}: 
	ln -sf ${EIRDIR}/eirmod_module_avltree.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_octree.${MOD}: 
	ln -sf ${EIRDIR}/eirmod_octree.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_parmmod.${MOD}:
	ln -sf ${EIRDIR}/eirmod_parmmod.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_precision.${MOD}:
	ln -sf ${EIRDIR}/eirmod_precision.${MOD} ${OBJDIR}
endif

${OBJDIR}/eirmod_extrab25.o:
	ln -sf ${EIRDIR}/eirmod_extrab25.o ${OBJDIR}

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

${OBJDIR}/eirmod_ccona.o:
	ln -sf ${EIRDIR}/eirmod_ccona.o ${OBJDIR}

${OBJDIR}/eirmod_ccoupl.o:
	ln -sf ${EIRDIR}/eirmod_ccoupl.o ${OBJDIR}

${OBJDIR}/eirmod_cestim.o:
	ln -sf ${EIRDIR}/eirmod_cestim.o ${OBJDIR}

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

${OBJDIR}/eirmod_module_avltree.o: 
	ln -sf ${EIRDIR}/eirmod_module_avltree.o ${OBJDIR}

${OBJDIR}/eirmod_octree.o: 
	ln -sf ${EIRDIR}/eirmod_octree.o ${OBJDIR}

${OBJDIR}/eirmod_parmmod.o:
	ln -sf ${EIRDIR}/eirmod_parmmod.o ${OBJDIR}

${OBJDIR}/eirmod_precision.o:
	ln -sf ${EIRDIR}/eirmod_precision.o ${OBJDIR}
else
${OBJDIR}/eirmod_braeir.${MOD}:
	touch ${OBJDIR}/eirmod_braeir.${MOD}

${OBJDIR}/eirmod_ccoupl.${MOD}:
	touch ${OBJDIR}/eirmod_ccoupl.${MOD}

${OBJDIR}/eirmod_comusr.${MOD}:
	touch ${OBJDIR}/eirmod_comusr.${MOD}

${OBJDIR}/eirmod_eirbra.${MOD}:
	touch ${OBJDIR}/eirmod_eirbra.${MOD}

${OBJDIR}/eirmod_extrab25.${MOD}:
	touch ${OBJDIR}/eirmod_extrab25.${MOD}

${OBJDIR}/eirmod_parmmod.${MOD}:
	touch ${OBJDIR}/eirmod_parmmod.${MOD}
endif

${MNEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA}
	${LD} ${LDOPTS} -o $@ $^ ${LDLIBES} ${LDOPTSend}

${OTEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA}
	${LD} ${LDOPTS} -o $@ $^ ${LDLIBES} ${LDOPTSend}

${GREXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA}
	${LD} ${LDOPTS} -o $@ $^ ${GRLIBES} ${LDLIBES} ${LDOPTSend}

${XDEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${SOLPS4OBJS} ${OBJDIR}/libb2.a ${MNEXTRA}
	${LD} ${LDOPTS} -o $@ $^ ${LCPP} ${GRLIBES} ${LDLIBES} ${LDEXTRA} ${LDOPTSend}

ifdef MDSPLUS_DIR
ifndef SOLPS_MDSPLUS_LIB
SOLPS_MDSPLUS_LIB=-L${MDSPLUS_DIR}/lib -lMdsLib_client 
endif
INCLUDE += -I${MDSPLUS_DIR}/include
DEFAULT: ${MDEXE}

${MDEXE}: ${OBJDIR}/%.exe: ${OBJDIR}/%.o ${OBJDIR}/libb2.a ${MNEXTRA}
	${LD} ${LDOPTS} -o $@ $^ ${LDLIBES} ${SOLPS_MDSPLUS_LIB} ${LDOPTSend}
endif

${OBJDIR}/libb2.a: ${LIBOBJS} ${SRCDIR}/include/git_version.h
	@${BLD} $@ ${LIBOBJS}

# target 'clean' cleans up the directory.
clean : 
	-mkdir ${OBJDIR}/.delete
	-mv -i ${OBJDIR}/*.o ${OBJDIR}/*.f ${OBJDIR}/*.a ${OBJDIR}/*.exe ${SRCDIR}/include/git_version.h ${OBJDIR}/.delete
ifneq (${MOD},o)
	-mv -i ${OBJDIR}/*.${MOD} ${OBJDIR}/.delete
endif
	-rm -rf ${OBJDIR}/.delete &

depend: ${OBJDIR}/LISTOBJ ${B2OBJS:.o=.F} ${EXELIST:.o=.F}
	@`which makedepend` -p'$${OBJDIR}/' ${DEFINES} -f- ${INCLUDE} $^ | \
	sed 's,^$${OBJDIR}/[^ ][^ ]*/,\$${OBJDIR}/,' > ${OBJDIR}/dependencies 
ifneq (${MOD},o)
ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
	@`which makedepend` -p'$${OBJDIR}/' ${DEFINES} -f- ${INCLUDE} ${SRCDIR}/modules.local/*.F ${SRCDIR}/modules/*.F -o.${MOD} | \
	sed 's,^$${OBJDIR}/[^ ][^ ]*/,\$${OBJDIR}/,' >> ${OBJDIR}/dependencies 
else
	@`which makedepend` -p'$${OBJDIR}/' ${DEFINES} -f- ${INCLUDE} ${SRCDIR}/modules/*.F -o.${MOD} | \
	sed 's,^$${OBJDIR}/[^ ][^ ]*/,\$${OBJDIR}/,' >> ${OBJDIR}/dependencies 
endif
endif
ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
	@egrep '^ {6,}use ' ${SRCLOCAL}/*.F ${SRCDIR}/*/*.F | grep -v 'IGNORE' | awk '{sub("\\.F:",".o:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"$$3".${MOD}"}' >> ${OBJDIR}/dependencies
ifneq (${MOD},o)
	@egrep '^ {6,}use ' ${SRCDIR}/modules.local/*.F ${SRCDIR}/modules/*.F | grep -v 'IGNORE' | awk '{sub("\\.F:",".${MOD}:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"$$3".${MOD}"}' >> ${OBJDIR}/dependencies
endif
else
	@egrep '^ {6,}use ' ${SRCDIR}/*/*.F | grep -v 'IGNORE' | awk '{sub("\\.F:",".o:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"$$3".${MOD}"}' >> ${OBJDIR}/dependencies
ifneq (${MOD},o)
	@egrep '^ {6,}use ' ${SRCDIR}/modules/*.F | grep -v 'IGNORE' | awk '{sub("\\.F:",".${MOD}:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"$$3".${MOD}"}' >> ${OBJDIR}/dependencies
endif
endif

tags:
ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
	rm -f ${SRCB2}/TAGS ; etags -o ${SRCB2}/TAGS ${SRCDIR}/modules.local/*.F ${SRCDIR}/include.local/*.* ${SRCDIR}/include/*.* ${SRCLOCAL}/*.F ${SRCDIR}/common/*.* ${SRCDIR}/common/COUPLE/*.F ${SRCDIR}/*/*.F
else
	rm -f ${SRCB2}/TAGS ; etags -o ${SRCB2}/TAGS ${SRCDIR}/include.local/*.* ${SRCDIR}/include/*.* ${SRCDIR}/common/*.* ${SRCDIR}/common/COUPLE/*.F ${SRCDIR}/*/*.F
endif

listobj:
ifdef USE_EIRENE
	@rm -f ${OBJDIR}/LISTOBJ; touch ${OBJDIR}/LISTOBJ; l="OBJS ="; \
	for d in `echo "${FPATH}" | tr : \ `; do \
		l="$$l `(cd $$d > /dev/null; echo *.F)`"; \
	done; \
	l="$$l `(cd ${SRCEIR}/modules > /dev/null; echo *.f)`"; \
	l="$$l `(cd ${SRCEIR}/interfaces/couple_SOLPS-ITER > /dev/null; echo *.F90)`"; \
	E="-e 's/\.F90/\.o/g' -e 's/\.F/\.o/g' -e 's/\.f/\.o/g'" ; for f in ${EXCLUDELIST}; do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$l" | eval sed "$$E" > ${OBJDIR}/LISTOBJ
else
	@rm -f ${OBJDIR}/LISTOBJ; touch ${OBJDIR}/LISTOBJ; l="OBJS ="; \
	for d in `echo "${FPATH}" | tr : \ `; do \
		l="$$l `(cd $$d > /dev/null; echo *.F)`"; \
	done; \
	E="-e 's/\.F/\.o/g'" ; for f in ${EXCLUDELIST}; do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$l" | eval sed "$$E" > ${OBJDIR}/LISTOBJ
endif
	@l="B2OBJS ="; \
	for d in `echo "${FPATH}" | tr : \ `; do \
		l="$$l `(cd $$d > /dev/null; echo *.F)`"; \
	done; \
	E="-e 's/\.F/\.o/g'" ; for f in ${EXCLUDELIST}; do \
		E="$$E -e 's/ $$f//'"; \
	done; \
	echo "$$l" | eval sed "$$E" >> ${OBJDIR}/LISTOBJ

${OBJDIR}/LISTOBJ: listobj

VERSION: ${SRCDIR}/include/git_version.h

${SRCDIR}/include/git_version.h: force
ifeq ($(shell [ -d ${SOLPSTOP} ] && echo yes || echo no ),yes)
	@echo "      character*15 :: gitversion ='`(cd ${SOLPSTOP}; git describe --dirty --always)`'" > ${SRCDIR}/include/git_version_new.h
else
	@echo "      character*15 :: gitversion ='`git describe --dirty --always`'" > ${SRCDIR}/include/git_version_new.h
endif
	@if cmp -s ${SRCDIR}/include/git_version_new.h ${SRCDIR}/include/git_version.h; then rm ${SRCDIR}/include/git_version_new.h; else mv ${SRCDIR}/include/git_version_new.h ${SRCDIR}/include/git_version.h; fi

${OBJDIR}/dependencies: ${SRCDIR}/modules/.new_modules
ifeq ($(shell [ -d ${OBJDIR} ] && echo yes || echo no ),no)
	-mkdir -p ${OBJDIR}
endif
	touch ${OBJDIR}/dependencies
	${MAKE} tags
	${MAKE} VERSION
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

ifeq ($(OBJECTCODE),g77)
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
	@echo VPATH=${VPATH}
	@echo FPATH=${FPATH}
	@echo OBJDIR=${OBJDIR}
#	@echo OBJS=${OBJS}
	@echo MODLIST=${MODLIST}
	@echo MODULES=${MODULES}
	@echo MODOBJS=${MODOBJS}
	@echo MODMODS=${MODMODS}
#	@echo $(filter-out ${MODOBJS},${ALLOBJS})
	@echo EXCLUDELIST=${EXCLUDELIST}
	@echo EXELIST=${EXELIST}
	@echo GREXE=${GREXE}
	@echo MNEXE=${MNEXE}
	@echo XDEXE=${XDEXE}
	@echo OTEXE=${OTEXE}
#	@echo ${SRCF}

local: 
	mkdir -p ${SRCLOCAL}
	echo "      subroutine b2local" > ${SRCLOCAL}/b2local.F
	echo "c" >> ${SRCLOCAL}/b2local.F
	echo "c store local or locally modified subroutines in this directory" >> ${SRCLOCAL}/b2local.F
	echo "c" >> ${SRCLOCAL}/b2local.F
	echo "      use b2mod_local" >> ${SRCLOCAL}/b2local.F
	echo '#include "b2local.h"' >> ${SRCLOCAL}/b2local.F
	echo "c" >> ${SRCLOCAL}/b2local.F
	echo "      end" >> ${SRCLOCAL}/b2local.F
	mkdir -p ${SRCDIR}/modules.local
	echo "      module b2mod_local" > ${SRCDIR}/modules.local/b2mod_local.F
	echo "c" >> ${SRCDIR}/modules.local/b2mod_local.F
	echo "c store local or locally modified modules in this directory" >> ${SRCDIR}/modules.local/b2mod_local.F
	echo "c" >> ${SRCDIR}/modules.local/b2mod_local.F
	echo "      end" >> ${SRCDIR}/modules.local/b2mod_local.F 
	mkdir -p ${SRCDIR}/include.local
	echo "c" > ${SRCDIR}/include.local/b2local.h
	echo "c store local or locally modified include files in this directory" >> ${SRCDIR}/include.local/b2local.h
	echo "c" >> ${SRCDIR}/include.local/b2local.h
