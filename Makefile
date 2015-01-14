

OBJDIR    = bin/${OBJECTCODE}
OBJDIREIR = ../Eirene/bin/${OBJECTCODE}
SRCDIRIO  = ../ioflush

.SUFFIXES: .o .f .f90 .F

FF=mpif90
#FF=/private/libs/intel_64/mpich2-1.3.2-install/bin/mpif90
CC=gcc
DBG=1

ifeq ("$(DBG)","")
     FFLAGS= -c -O2 -r8
     CFLAGS= -c -O2
else
     FFLAGS= -c -g -r8 -traceback -check
     CFLAGS= -c -g -check
endif

FOPTS=
COPTS=
LIBS=

VPATH=./src/modules:./src/documentation:./src/equations:./src/input:./src/output:./src/sources:./src/transport:./src/solvers:./src/utility:./Add_Mod/Utility:./src/user:./src/b2aux:./src/driver:./src/b2plot_min

MODULES= b2mod_types.o  b2mod_layer.o b2mod_tallies.o b2mod_wall.o b2mod_constants.o b2mod_b2cmpa.o b2mod_user_namelist.o b2mod_elements.o b2mod_indirect.o \
	b2mod_boundary_namelist.o b2mod_version.o b2mod_b2cmpt.o b2mod_sources.o b2mod_time.o b2mod_geo.o b2mod_residuals.o b2mod_eirene_globals.o \
        b2mod_neutr_src_scaling.o \
	b2mod_neutrals_namelist.o  b2mod_transport.o \
	b2mod_anomalous_transport.o  b2mod_geo_corner.o         b2mod_numerics_namelist.o  b2mod_transport_models.o \
	b2mod_b2cmfs.o               b2mod_transport_namelist.o\
	b2mod_b2cmgs.o               b2mod_input_profile.o      \
	b2mod_plasma.o             \
	b2mod_b2cmpb.o               b2mod_diag.o               b2mod_ppout.o              \
	b2mod_ma28.o               b2mod_rates.o              \
	b2mod_b2cmrc.o               b2mod_ma28_for_7diag.o     \
	b2mod_ma28_for_9diag.o     b2mod_solpstop.o           b2mod_work.o\
	b2mod_ysmp_sdrv.o\
	b2mod_b2plot.o               b2mod_feedback.o           b2mod_sputter.o\
	b2mod_first_flight.o      \
	b2mod_astra_to_b2.o b2mod_b2_to_astra.o b2mod_nclass.o b2mod_neoclassical.o b2mod_ranges.o b2mod_b2plot_debug.o

OBJECTS= documentation.o equations.o input.o output.o sources.o transport.o solvers.o utility.o add_utility.o user.o b2aux.o driver.o b2plot.o

DESTM = $(MODULES:%.o=$(OBJDIR)/%.o)
DESTO = $(OBJECTS:%.o=$(OBJDIR)/%.o)

INCDIR= -I ./src/include -I ./src/common -I ./src/common/COUPLE -I ${OBJDIREIR}
LIBMOD= -L${OBJDIR} -lb25 -L${OBJDIREIR} -leirene
LIBGR=-L/usr/lib64 -L/b2-ext/ncarg/lib -Llib/grgli/grsoft -lgr -L../../lib/${OBJECTCODE} -Llib/grgli/gligks -lgks -lX11 -lXt

#default deps
$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) $(COPTS) -c -o $@ $<
$(OBJDIR)/%.o: %.f
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/%.o: %.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/%.o: %.f90
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/%.F90.o: %.F90
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -module $(OBJDIR) -c -o $@ $<


b2eirene: dirs $(OBJDIR)/b2eirene

$(OBJDIR)/b2eirene: ${OBJDIREIR}/libeirene.a ${OBJDIR}/libb25.a $(OBJDIR)/b2mn_mpi.o $(OBJDIR)/ioflush.o
	rm -f $(OBJDIR)/b2eirene
	$(FF) -o $(OBJDIR)/b2eirene $(INCDIR) $(OBJDIR)/b2mn_mpi.o $(OBJDIR)/ioflush.o $(LIBMOD) $(LIBGR) $(LIBMPI) 

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*.a $(OBJDIR)/b2eirene

dirs: $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

#deps
$(OBJDIR)/b2mn_mpi.o:
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -module $(OBJDIR) -c -o $@ src/driver/b2mn_mpi.F

# NOTE: is this ioflush needed?
$(OBJDIR)/ioflush.o:
	$(CC) $(CFLAGS) $(COPTS) -c -o $@ ${SRCDIRIO}/ioflush.c

#libinterface:

${OBJDIREIR}/libeirene.a:

###############################################################################

# NOTE: not all *.o -files needed for libb25 !?
$(OBJDIR)/libb25.a: $(DESTM) $(DESTO)
	ar vr $(OBJDIR)/libb25.a $(OBJDIR)/*.o


#deps
$(OBJDIR)/documentation.o: $(DESTM) $(addprefix $(OBJDIR)/, \
	b2cdca.o  b2cdci.o  b2cdcr.o  b2cdcv.o)

$(OBJDIR)/equations.o: $(DESTM) $(addprefix $(OBJDIR)/, \
	b2news.o   b2npc9.o  b2nph9.o  b2npmo.o  b2nppo.o  b2nxdp.o  b2nxdv.o  b2nxfm.o  b2nxfx.o  b2nxst.o \
	b2news_.o  b2npco.o  b2npht.o  b2npp7.o  b2nxcm.o  b2nxdu.o  b2nxfc.o  b2nxfv.o  )

$(OBJDIR)/input.o: $(DESTM) $(addprefix $(OBJDIR)/, \
	b2rfcp.o  b2rflb.o  b2rucp.o  b2rugm.o  b2rups.o  b2rurc.o  b2rusr.o  b2ruzd.o)

$(OBJDIR)/output.o: $(DESTM) $(addprefix $(OBJDIR)/, \
	b2wdat.o  b2wfpj.o  \
	b2wfcp.o  b2wucp.o  b2wuzd.o \
	b2wfgm.o  \
	b2wfpi.o  b2wups.o  tallies.o)

$(OBJDIR)/sources.o: $(DESTM) $(addprefix $(OBJDIR)/, \
	b2sifr.o   b2sihs_.o  b2sqel.o  b2srst.o   b2stbr_bas.o   eirene_f30f31.o   ggfill.o      setwrk0.o \
	average.o  b2sifr_.o  b2spcx.o   b2sral.o  b2stbc.o       b2stbc_spb.o     b2stbr_phys.o  eirene_mc.o       heatdiff1D.o\
	b2siav.o   b2sigp.o   b2spel.o   b2srdt.o  b2stbc_bas.o   b2stbm.o         b2stcx.o       heatdiff2D.o\
	b2sicf.o   b2sihs.o   b2sqcx.o   b2srsm.o  b2stbc_phys.o  b2stbr.o         b2stel.o       eseec0.o          integrate.o \
	eirene_mc_init.o)

$(OBJDIR)/transport.o: $(DESTM) $(addprefix $(OBJDIR)/, \
	b2tfcc.o   b2tfhi.o   b2tiner.o  b2tlhe.o  b2tqca.o  b2tr21.o  b2trno.o   b2tvspa.o  b2txfx.o\
	b2tfch.o   b2tfhi_.o  b2tinnt.o  b2tlhi.o  b2tqce.o  b2tral.o  b2trql.o   b2txfy.o\
	b2tanml.o  b2tfhe.o   b2tfnb.o   b2tlc0.o   b2tlmv.o  b2tqin.o  b2trcl.o  b2tstch.o  b2txcx.o   b2txsx.o\
	b2tcpa.o   b2tfhe_.o  b2tfrn.o   b2tlh0.o   b2tlnl.o  b2tqna.o  b2treq.o  b2ttia.o   b2txcy.o   b2txsy.o\
	bp.o class.o colxi.o fluxav.o menn.o neo_diagnostics.o neo_validity.o neoart.o neochi.o \
	neodv.o penq.o perr.o ps.o viscol.o viscos.o)

$(OBJDIR)/solvers.o: $(DESTM) $(addprefix $(OBJDIR)/, \
	b2upco.o  b2upht_.o  b2urmo.o  b2ursd.o  b2usco.o  b2usht.o  b2usp7.o   b2uspo.o  b2ux7p.o  b2uxm9.o\
	b2upht.o  b2uppo.o   b2ursc.o  b2usc9.o  b2ush9.o  b2usmo.o  b2usp7_.o  b2ux5p.o  b2ux9p.o)

$(OBJDIR)/utility.o: $(DESTM) $(addprefix $(OBJDIR)/, \
	cfwuin.o    ini1_7_solps.o     my_outi.o    slv5pt.o        xerrab.o\
	cfwure.o    ini1_solps.o    lefta.o          myblas.o     ratio.o     smax.o          xerset.o\
	chcase.o    int2d.o         len_of_digits.o  nagsubst.o   smin.o          xertst.o\
	hbout.o            intcor.o        lwimai.o         xxdata_11.o\
	bfill.o           daytim.o    ifill.o            intfaceh.o      lwmain.o         open_file.o  ssum.o          ysmp_d.o\
	dseval.o    illtern_7_solps.o  ma28.o           stabeq.o\
	dspline.o   ipgeti.o        ma28copy.o       streql.o\
	cfopen.o          ilu5g.o            ipgetr.o        samax.o     strip_spaces.o\
	cfruch.o          ipmain.o        machsfr.o        prgend.o     strmas.o\
	cfruin.o          get_jsep.o  iluter.o           ipos.o          prgini.o     subsys.o\
	cfrure.o          ilutern_7_solps.o  mstep.o          prvrt.o      sfill.o     sysend.o\
	cfvers.o          ilutern_9_solps.o  my_intcor.o      prvrti.o     sysini.o\
	cfwuch.o          ilutern_solps.o    jobnam.o        my_out.o         ratadas.o    sip5g.o     usrnam.o \
        intfacev.o\
	dfmin.o interp2d.o lubksb.o ludcmp.o sfluxav.o tfluxav.o)

$(OBJDIR)/add_utility.o: $(DESTM) $(addprefix $(OBJDIR)/, utility.o \
	dgbtf2.o  dgemm.o  dgetf2.o  dscal.o  dtrmm.o  dtrti2.o     idamax.o  xerbla.o\
	dgemv.o  dgetrf.o  dlamch.o  dswap.o  dtrmv.o  dtrtri.o     f01aaf_my.o  ilaenv.o  sdot_my.o\
	dcopy.o   dgbtrf.o   dgbtrs.o  dger.o   dgetri.o  dlaswp.o  dtbsv.o  dtrsm.o  lsame.o   )

$(OBJDIR)/user.o: $(DESTM) $(addprefix $(OBJDIR)/, \
	b2blnc.o  b2blne.o  b2blnm.o  b2file.o  b2trace.o  b2usrtrc.o  b2wrint.o  b2wrsep.o  )

$(OBJDIR)/b2aux.o: $(DESTM) $(addprefix $(OBJDIR)/, \
	b2xbzb.o  b2xpfe.o  b2xpne.o  b2xpnr.o  b2xppr.o  b2xvcp.o  b2xvfx.o  b2xvps.o  b2xxid.o\
	b2xgbs.o  b2xpfi.o  b2xpni.o  b2xppb.o  b2xppz.o  b2xvff.o  b2xvfy.o  b2xvsg.o  b2xxmm.o\
	b2xpnm.o  b2xppe.o  b2xpro.o  b2xvfv.o  b2xvgm.o  b2xxgs.o  b2xzdd.o)

$(OBJDIR)/driver.o: $(DESTM) $(addprefix $(OBJDIR)/, \
	b2mndr.o  b2mndt.o  b2mwmv.o  b2mwqt.o  b2mxac.o  b2mxnp.o  b2mxzr.o  \
	b2mnds.o  b2mwit.o  b2mwq0.o  b2mwti.o  b2mxar.o  b2mxnu.o  )

$(OBJDIR)/b2plot.o: $(DESTM) $(addprefix $(OBJDIR)/, \
	chord.o lower_case.o init_wall.o mapx.o mapy.o species.o)

#specific deps
$(OBJDIR)/b2mod_neutrals_namelist.o: src/modules/b2mod_neutrals_namelist.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DUSE_EIRENE -DNO_SAVE -module $(OBJDIR)  -c -o $@ $<
$(OBJDIR)/b2mod_diag.o: src/modules/b2mod_diag.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DUSE_EIRENE -module $(OBJDIR)  -c -o $@ $<
$(OBJDIR)/b2mod_work.o: src/modules/b2mod_work.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DUSE_EIRENE -module $(OBJDIR)  -c -o $@ $<
$(OBJDIR)/b2mod_geo.o: src/modules/b2mod_wall.F $(OBJDIR)/b2mod_geo_corner.o
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -module $(OBJDIR) -c -o $@ src/modules/b2mod_geo.F
$(OBJDIR)/b2mod_wall.o: src/modules/b2mod_wall.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_CDF -module $(OBJDIR)  -c -o $@ $<
$(OBJDIR)/b2mod_tallies.o: src/modules/b2mod_tallies.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_CDF -module $(OBJDIR)  -c -o $@ $<
$(OBJDIR)/b2mod_movies.o: src/modules/b2mod_movies.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_CDF -module $(OBJDIR)  -c -o $@ $<
$(OBJDIR)/b2stbr.o: src/sources/b2stbr.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DUSE_EIRENE -DNO_CDF -module $(OBJDIR)  -c -o $@ $<
$(OBJDIR)/eirene_f30f31.o: src/sources/eirene_f30f31.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_SAVE -DNO_SHORT_CYCLE -DNO_FORT10 -DUSE_EIRENE -DNO_CDF -DNO_NAG -module $(OBJDIR)  -c -o $@ $<
$(OBJDIR)/eirene_mc.o: src/sources/eirene_mc.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_SAVE -DNO_SHORT_CYCLE -DNO_FORT10 -DUSE_EIRENE -DNO_CDF -DNO_NAG -module $(OBJDIR)  -c -o $@ $<
$(OBJDIR)/eirene_mc_init.o: src/sources/eirene_mc_init.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_SAVE -DNO_SHORT_CYCLE -DNO_FORT10 -DUSE_EIRENE -DNO_CDF -DNO_NAG -module $(OBJDIR)  -c -o $@ $<
$(OBJDIR)/heatdiff1D.o: src/sources/heatdiff1D.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_NAG  -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/heatdiff2D.o: src/sources/heatdiff2D.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_NAG  -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/eseec0.o: src/sources/eseec0.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_NAG  -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/ggfill.o: src/sources/ggfill.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_NAG  -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/integrate.o: src/sources/integrate.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_NAG  -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/setwrk0.o: src/sources/setwrk0.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_NAG  -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/average.o: src/sources/average.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_NAG  -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/b2tr21.o: src/transport/b2tr21.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_NAG -Dsgetrf=dgetrf -Dsgetri=dgetri  -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/nagsubst.o: src/utility/nagsubst.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_NAG  -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/slv5pt.o: src/utility/slv5pt.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -Dsgbtrf=dgbtrf -Dsgbtrs=dgbtrs  -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/tallies.o: src/output/tallies.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS)  -DNO_CDF -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/b2mwti.o: src/driver/b2mwti.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DNO_CDF  -module $(OBJDIR) -c -o $@ $<
$(OBJDIR)/b2mndr.o: src/driver/b2mndr.F
	$(FF) $(FFLAGS) $(INCDIR) $(OPTS) -DUSE_EIRENE -DNO_CDF  -module $(OBJDIR) -c -o $@ $<



