DEFINES = -DWANT_THIS ${SOLPS_CPP}
ifdef USE_EIRENE
DEFINES += ${USE_EIRENE}
SRCEIR = ${SOLPSTOP}/src/Eirene
ifdef USE_MPI
EIRDIR = ${SOLPSTOP}/bin/${OBJECTCODE}/Eirene
else
EIRDIR = ${SOLPSTOP}/bin/${OBJECTCODE}/Eirene.nompi
endif
endif
ifdef USE_MPI
DEFINES += ${USE_MPI}
endif
ifdef PERFMON
DEFINES += ${PERFMON}
endif

SRCB2 = ${PWD}
SRCDIR = ${SRCB2}/src
ifdef USE_MPI
OBJDIR = bin/${OBJECTCODE}
else
OBJDIR = bin/${OBJECTCODE}.nompi
endif
BASEDIR = ${OBJDIR}
SRCLOCAL = ${SRCB2}/src.local
ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
INCLUDE = -I${SRCLOCAL}
else
INCLUDE =
endif
INCLUDE += -I${SRCDIR}/common -I${SRCDIR}/include.local -I${SRCDIR}/include

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
PROG_OT = b2ag.exe b2ah.exe b2ai.exe b2ar.exe b2co.exe b2uf.exe b2fu.exe b2ts.exe b2yi_gnuplot.exe b2yh.exe b2yt.exe b2yv.exe
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

.PHONY: DEFAULT NOPLOT ALL VERSION clean realclean depend listobj tags

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
VPATH=${FPATH}:${SRCEIR}/modules:${SRCEIR}/extraB25
MODLIST+=${SRCEIR}/modules/*.f ${SRCEIR}/extraB25/eirmod_*.F90
MNEXTRA=${EIRDIR}/libeirene.a ${EIRDIR}/libgr_dummy.a
else
# MNEXTRA=${EIRDIR}/eirmod_precision.o ${EIRDIR}/eirmod_braeir.o ${EIRDIR}/eirmod_ccoupl.o ${EIRDIR}/eirmod_parmmod.o ${EIRDIR}/eirmod_comprt.o ${EIRDIR}/eirmod_comusr.o
# EXCLUDELIST += ${patsubst ${OBJDIR}/%.o, %.o, ${MNEXTRA} }
endif

MODULES = ${patsubst %.F %.f %.F90,%.o,${shell echo ${MODLIST} } }
MODMODS = ${MODULES:%.o=${OBJDIR}/%.${MOD}}
MODOBJS = ${MODULES:%.o=${OBJDIR}/%.o}

ifdef USE_EIRENE
${OBJDIR}/libgr_dummy.a:
	ln -sf ${EIRDIR}/libgr_dummy.a ${OBJDIR}

${OBJDIR}/libeirene.a:
	ln -sf ${EIRDIR}/libeirene.a ${OBJDIR}

ifneq (${MOD},o)
${OBJDIR}/eirmod_extrab25.${MOD}:
	ln -sf ${EIRDIR}/eirmod_extrab25.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_braeir.${MOD}: 
	ln -sf ${EIRDIR}/eirmod_braeir.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_ccoupl.${MOD}:
	ln -sf ${EIRDIR}/eirmod_ccoupl.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_comprt.${MOD}:
	ln -sf ${EIRDIR}/eirmod_comprt.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_comusr.${MOD}:
	ln -sf ${EIRDIR}/eirmod_comusr.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_eirbra.${MOD}: ${OBJDIR}/eirmod_precision.${MOD} ${OBJDIR}/eirmod_parmmod.${MOD} 
	ln -sf ${EIRDIR}/eirmod_eirbra.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_parmmod.${MOD}:
	ln -sf ${EIRDIR}/eirmod_parmmod.${MOD} ${OBJDIR}

${OBJDIR}/eirmod_precision.${MOD}:
	ln -sf ${EIRDIR}/eirmod_precision.${MOD} ${OBJDIR}
endif

${OBJDIR}/eirmod_extrab25.o:
	ln -sf ${EIRDIR}/eirmod_extrab25.o ${OBJDIR}

${OBJDIR}/eirmod_braeir.o: 
	ln -sf ${EIRDIR}/eirmod_braeir.o ${OBJDIR}

${OBJDIR}/eirmod_ccoupl.o:
	ln -sf ${EIRDIR}/eirmod_ccoupl.o ${OBJDIR}

${OBJDIR}/eirmod_comprt.o:
	ln -sf ${EIRDIR}/eirmod_comprt.o ${OBJDIR}

${OBJDIR}/eirmod_comusr.o:
	ln -sf ${EIRDIR}/eirmod_comusr.o ${OBJDIR}

${OBJDIR}/eirmod_eirbra.o: ${OBJDIR}/eirmod_precision.o ${OBJDIR}/eirmod_parmmod.o
	ln -sf ${EIRDIR}/eirmod_eirbra.o ${OBJDIR}

${OBJDIR}/eirmod_parmmod.o:
	ln -sf ${EIRDIR}/eirmod_parmmod.o ${OBJDIR}

${OBJDIR}/eirmod_precision.o:
	ln -sf ${EIRDIR}/eirmod_precision.o ${OBJDIR}
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

ifeq (${MOD},o)
${OBJDIR}/libb2.a: $(filter-out ${MODOBJS},${ALLOBJS})
else
${OBJDIR}/libb2.a: ${ALLOBJS}
endif
	${BLD} $@ $?

# targets 'clean' and 'realclean' clean up the directory.
clean : 
	-mkdir ${OBJDIR}/.delete
	-mv -i ${OBJDIR}/*.o ${OBJDIR}/*.f ${OBJDIR}/*.a ${OBJDIR}/*.exe ${OBJDIR}/.delete
ifneq (${MOD},o)
	-mv -i ${OBJDIR}/*.${MOD} ${OBJDIR}/.delete
endif
	-rm -rf ${OBJDIR}/.delete &

depend: ${OBJDIR}/LISTOBJ ${B2OBJS:.o=.F} ${EXELIST:.o=.F}
	@`which makedepend` -p'$${OBJDIR}/' ${DEFINES} -f- ${INCLUDE} $^ | \
	sed 's,^$${OBJDIR}/[^ ][^ ]*/,\$${OBJDIR}/,' > ${OBJDIR}/dependencies 
ifneq (${MOD},o)
	@`which makedepend` -p'$${OBJDIR}/' ${DEFINES} -f- ${INCLUDE} $^ -o.${MOD} | \
	sed 's,^$${OBJDIR}/[^ ][^ ]*/,\$${OBJDIR}/,' >> ${OBJDIR}/dependencies 
endif
ifdef USE_EIRENE
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
else
ifeq ($(shell [ -d ${SRCLOCAL} ] && echo yes || echo no ),yes)
	@egrep '^ {6,}use ' ${SRCLOCAL}/*.F ${SRCDIR}/*/*.F | grep -v 'IGNORE' | grep -v eirmod | awk '{sub("\\.F:",".o:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"$$3".${MOD}"}' >> ${OBJDIR}/dependencies
ifneq (${MOD},o)
	@egrep '^ {6,}use ' ${SRCDIR}/modules.local/*.F ${SRCDIR}/modules/*.F | grep -v eirmod | grep -v 'IGNORE' | awk '{sub("\\.F:",".${MOD}:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"$$3".${MOD}"}' >> ${OBJDIR}/dependencies
endif
else
	@egrep '^ {6,}use ' ${SRCDIR}/*/*.F | grep -v 'IGNORE' | grep -v eirmod | awk '{sub("\\.F:",".o:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"$$3".${MOD}"}' >> ${OBJDIR}/dependencies
ifneq (${MOD},o)
	@egrep '^ {6,}use ' ${SRCDIR}/modules/*.F | grep -v 'IGNORE' | grep -v eirmod | awk '{sub("\\.F:",".${MOD}:",$$1);sub("^.*/","$${OBJDIR}/",$$1); print $$1,"$${OBJDIR}/"$$3".${MOD}"}' >> ${OBJDIR}/dependencies
endif
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
	l="$$l `(cd ${SRCEIR}/extraB25 > /dev/null; echo *.F90)`"; \
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

VERSION: ${SRCLOCAL}/git_version.h

${SRCLOCAL}/git_version.h:
ifeq ($(shell [ -d ${SOLPSTOP} ] && echo yes || echo no ),yes)
	echo "      character*15 :: gitversion ='`(cd ${SOLPSTOP}; git describe --dirty --always)`'" > ${SRCLOCAL}/git_version.h
else
	echo "      character*15 :: gitversion ='`git describe --dirty --always`'" > ${SRCLOCAL}/git_version.h
endif

${OBJDIR}/dependencies: ${SRCDIR}/modules/.new_modules
ifeq ($(shell [ -d ${OBJDIR} ] && echo yes || echo no ),no)
	-mkdir -p ${OBJDIR}
endif
	touch ${OBJDIR}/dependencies
	${MAKE} tags
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
#	@echo OBJS=${OBJS}
#	@echo SOBJS=${SOBJS}
#	@echo DOBJS=${DOBJS}
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
#	@echo ${OBJDEST}

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
