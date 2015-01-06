#FC = gfortran 
FC=ifort 

OBJ0  =  eirmod_precision.o eirmod_parmmod.o eirmod_braeir.o eirmod_eirbra.o eirmod_extrab25.o
OBJ1  =  b2mod_types.o b2mod_constants.o b2mod_version.o b2mod_b2cmfs.o b2mod_b2cmgs.o b2mod_b2cmpa.o b2mod_b2cmpb.o b2mod_b2cmpt.o b2mod_b2cmrc.o b2mod_b2cmwg.o b2mod_anomalous_transport.o b2mod_eirene_globals.o 
OBJ2  =  b2mod_indirect.o b2mod_boundary_namelist.o b2mod_elements.o b2mod_eirene_globals.o b2mod_geo_corner.o b2mod_geo.o b2mod_layer.o b2mod_wall.o
OBJ3  =  b2mod_numerics_namelist.o b2mod_plasma.o b2mod_ppout.o b2mod_rates.o b2mod_residuals.o b2mod_movies.o 
OBJ4  =  b2mod_sources.o b2mod_solpstop.o b2mod_sputter.o b2mod_tallies.o b2mod_user_namelist.o b2mod_neutrals_namelist.o b2mod_neutr_src_scaling.o b2mod_time.o b2mod_diag.o
OBJ5  =  b2mod_transport.o b2mod_work.o  b2mod_b2plot.o b2mod_b2plot_debug.o b2mod_b2plot_help.o b2mod_b2plot_values.o b2mod_b2plot_xy.o b2mod_first_flight.o
OBJ6  =  b2tlnl.o b2spcx.o b2spel.o b2sqcx.o b2sqel.o eirene_f30f31.o eirene_mc_init.o heatdiff1D.o
OBJ7  =  cfopen.o cfruch.o cfruin.o cfrure.o cfvers.o cfwure.o cfwuin.o chcase.o daytim.o gfsub2.o get_jsep.o ifill.o ipgeti.o ipgetr.o
OBJ8  =  intfaceh.o intfacev.o ipmain.o jobnam.o len_of_digits.o lwimai.o lwmain.o
OBJ9  =  machsfr.o mstep.o myblas.o my_intcor.o my_out.o nagsubst.o open_file.o  
OBJ10 =  prgend.o prgini.o prvrt.o  prvrti.o rdneutrs.o samax.o sfill.o smax.o smin.o
OBJ11 =  sortda.o ssum.o streql.o subsys.o sysend.o sysini.o usrnam.o xerrab.o xerset.o xertst.o
OBJ12 =  b2rfcp.o b2rflb.o b2rucp.o b2rugm.o b2rups.o b2rurc.o b2ruzd.o
OBJ13 =  b2xgbs.o b2xpne.o b2xpni.o b2xvcp.o b2xvff.o b2xvfx.o b2xvfy.o b2xvgm.o b2xvps.o b2xvsg.o b2xxgs.o b2xxid.o
OBJ14 =  additional_elements.o ank_interface.o avey.o caltraj.o ceil.o chord.o clear_extralabel.o clip_box.o clip_line.o clvd.o clvl.o clvn.o clvp.o clvs.o color.o cont.o contourplot.o copy.o  coronal.o
OBJ15 =  dhdata.o dropstackbottom.o end_page.o extend.o getint.o gridon.o gslwsc_d.o init_page.o init_wall.o initialize_plotting.o ipolate.o ipop.o fill.o ipush.o isadigit.o isnegative.o iswap.o isplusorminus.o floor.o flux.o flx.o fly.o font.o fx.o fy.o
OBJ16 =  labelon.o lblalongcl.o lblbartraj.o line.o lined_d.o linewidth.o lower_case.o map.o map_chords.o mapx.o mapy.o
OBJ17 =  mesh.o meshit.o ngread.o nice.o nofill.o outlineon.o page.o process_chord_slice.o process_chords.o process_chords_boxcar.o pswa.o psx.o psy.o pvol.o readline.o rlcomp.o rlshdw.o rpop.o rpush.o rswap.o separatrixon.o set_a4l.o set_a4p.o set_extralabel.o
OBJ18 =  set_globalheader.o set_midplane_map.o set_new_page.o set_same_page.o set_square.o shrink_label.o 
OBJ19 =  skipblanks.o species.o spect.o spop.o spush.o sumx.o sumy.o strswap.o surf.o symbol.o targeton.o terminate_plotting.o trgsumw.o trgsurf.o trgw.o trgz.o vesselon.o wallon.o
OBJ20 =  b2pinp.o b2paux.o b2pldr.o b2plot.o b2pnwl.o b2prdl.o b2psh.o b2ptrdl.o b2ptrg.o b2ptrgl.o b2pwinm.o b2pwinp.o b2pwld.o b2pwlprp.o b2pwrdld.o b2pwsort.o b2mwti.o mesh.o meshit.o plot_chords.o plot_wall_elements.o plottarget.o plotcontour.o plotpatch.o pltoutline.o pltsptrix.o pne.o plttarget.o pltvessel.o pltwall.o
OBJ21 =  b2pl.o

NCARG=/b2-ext/ncarg/lib
#NCARG=../../libs/ncarg

MSCL=lib/mscl/linux.ifort64
#MSCL=../../libs/mscl/linux.intel

b2pl.exe: $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9) $(OBJ10) $(OBJ11) $(OBJ12) $(OBJ13) $(OBJ14) $(OBJ15) $(OBJ16) $(OBJ17) $(OBJ18) $(OBJ19) $(OBJ20) $(OBJ21)
#	$(FC) *.o -o Release/b2pl.exe  $(MSCL)/libmscl.a -L$(NCARG)/lib -lncarg -lncarg_gks -lncarg_c -lngmath -lcgm -lncarg_ras -lcairo -lfreetype -L/usr/lib64/ -lX11
	$(FC) *.o -o Release/b2pl.exe  $(MSCL)/libmscl.a $(NCARG)/lib/libncarg.a $(NCARG)/lib/libncarg_gks.a $(NCARG)/lib/libncarg_c.a $(NCARG)/lib/libngmath.a $(NCARG)/lib/libcgm.a $(NCARG)/lib/libncarg_ras.a -lcairo -lfreetype -L/usr/lib64/ -lX11
#	$(FC) *.o -o Release/b2pl.exe /afs/ipp-garching.mpg.de/common/soft/ncar.4.0/@sys/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm -lsocket -lf77compat
eirmod_precision.o: ../../Eirene_dummy/modules/eirmod_precision.f
	$(FC) -c ../../Eirene_dummy/modules/eirmod_precision.f
eirmod_parmmod.o: ../../Eirene_dummy/modules/eirmod_parmmod.f
	$(FC) -c ../../Eirene_dummy/modules/eirmod_parmmod.f
eirmod_braeir.o: ../../Eirene_dummy/modules/eirmod_braeir.f
	$(FC) -c ../../Eirene_dummy/modules/eirmod_braeir.f
eirmod_eirbra.o: ../../Eirene_dummy/modules/eirmod_eirbra.f
	$(FC) -c ../../Eirene_dummy/modules/eirmod_eirbra.f
eirmod_extrab25.o: ../../Eirene_dummy/extraB25/eirmod_extraB25.F90
	$(FC) -c ../../Eirene_dummy/extraB25/eirmod_extraB25.F90
b2mod_types.o: ../src/modules/b2mod_types.F
	$(FC) -c ../src/modules/b2mod_types.F
b2mod_constants.o: ../src/modules/b2mod_constants.F
	$(FC) -c ../src/modules/b2mod_constants.F
b2mod_version.o: ../src/modules/b2mod_version.F
	$(FC) -c ../src/modules/b2mod_version.F
b2mod_b2cmfs.o: ../src/modules/b2mod_b2cmfs.F
	$(FC) -c ../src/modules/b2mod_b2cmfs.F
b2mod_b2cmgs.o: ../src/modules/b2mod_b2cmgs.F
	$(FC) -c ../src/modules/b2mod_b2cmgs.F
b2mod_b2cmpa.o: ../src/modules/b2mod_b2cmpa.F
	$(FC) -c ../src/modules/b2mod_b2cmpa.F
b2mod_b2cmpb.o: ../src/modules/b2mod_b2cmpb.F
	$(FC) -c ../src/modules/b2mod_b2cmpb.F
b2mod_b2cmpt.o: ../src/modules/b2mod_b2cmpt.F
	$(FC) -c ../src/modules/b2mod_b2cmpt.F
b2mod_b2cmrc.o: ../src/modules/b2mod_b2cmrc.F
	$(FC) -c ../src/modules/b2mod_b2cmrc.F
b2mod_b2cmwg.o: ../src/modules/b2mod_b2cmwg.F
	$(FC) -c ../src/modules/b2mod_b2cmwg.F
b2mod_anomalous_transport.o: ../src/modules/b2mod_anomalous_transport.F
	$(FC) -c ../src/modules/b2mod_anomalous_transport.F
b2mod_indirect.o: ../src/modules/b2mod_indirect.F
	$(FC) -c ../src/modules/b2mod_indirect.F
b2mod_boundary_namelist.o: ../src/modules/b2mod_boundary_namelist.F
	$(FC) -c ../src/modules/b2mod_boundary_namelist.F -I../src/include/
b2mod_elements.o: ../src/modules/b2mod_elements.F
	$(FC) -c ../src/modules/b2mod_elements.F
b2mod_eirene_globals.o: ../src/modules/b2mod_eirene_globals.F
	$(FC) -c ../src/modules/b2mod_eirene_globals.F
b2mod_geo_corner.o: ../src/modules/b2mod_geo_corner.F
	$(FC) -c ../src/modules/b2mod_geo_corner.F
b2mod_geo.o: ../src/modules/b2mod_geo.F
	$(FC) -c ../src/modules/b2mod_geo.F
b2mod_layer.o: ../src/modules/b2mod_layer.F
	$(FC) -c ../src/modules/b2mod_layer.F -I../src/include/
b2mod_movies.o: ../src/modules/b2mod_movies.F
	$(FC) -c -DNO_CDF ../src/modules/b2mod_movies.F
b2mod_numerics_namelist.o: ../src/modules/b2mod_numerics_namelist.F
	$(FC) -c ../src/modules/b2mod_numerics_namelist.F -I../src/include/
b2mod_plasma.o: ../src/modules/b2mod_plasma.F
	$(FC) -c ../src/modules/b2mod_plasma.F
b2mod_ppout.o: ../src/modules/b2mod_ppout.F
	$(FC) -c ../src/modules/b2mod_ppout.F
b2mod_rates.o: ../src/modules/b2mod_rates.F
	$(FC) -c ../src/modules/b2mod_rates.F
b2mod_first_flight.o: ../src/modules/b2mod_first_flight.F
	$(FC) -c ../src/modules/b2mod_first_flight.F -I../src/common/
b2mod_residuals.o: ../src/modules/b2mod_residuals.F
	$(FC) -c ../src/modules/b2mod_residuals.F
b2mod_sources.o: ../src/modules/b2mod_sources.F
	$(FC) -c ../src/modules/b2mod_sources.F
b2mod_solpstop.o: ../src/modules/b2mod_solpstop.F
	$(FC) -c ../src/modules/b2mod_solpstop.F
b2mod_sputter.o: ../src/modules/b2mod_sputter.F
	$(FC) -c ../src/modules/b2mod_sputter.F -I../src/include/
b2mod_tallies.o: ../src/modules/b2mod_tallies.F
	$(FC) -c -DNO_CDF ../src/modules/b2mod_tallies.F -I../src/common/
b2mod_time.o: ../src/modules/b2mod_time.F
	$(FC) -c ../src/modules/b2mod_time.F
b2mod_neutrals_namelist.o: ../src/modules/b2mod_neutrals_namelist.F
	$(FC) -c ../src/modules/b2mod_neutrals_namelist.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
b2mod_neutr_src_scaling.o: ../src/modules/b2mod_neutr_src_scaling.F
	$(FC) -c ../src/modules/b2mod_neutr_src_scaling.F -I../src/include/
b2mod_transport.o: ../src/modules/b2mod_transport.F
	$(FC) -c ../src/modules/b2mod_transport.F
b2mod_user_namelist.o: ../src/modules/b2mod_user_namelist.F
	$(FC) -c ../src/modules/b2mod_user_namelist.F -I../src/include/
b2mod_work.o: ../src/modules/b2mod_work.F
	$(FC) -c ../src/modules/b2mod_work.F
b2mod_b2plot.o: ../src/modules/b2mod_b2plot.F
	$(FC) -c ../src/modules/b2mod_b2plot.F -I../src/include/
b2mod_b2plot_debug.o: ../src/modules/b2mod_b2plot_debug.F
	$(FC) -c ../src/modules/b2mod_b2plot_debug.F
b2mod_b2plot_help.o: ../src/modules/b2mod_b2plot_help.F
	$(FC) -c ../src/modules/b2mod_b2plot_help.F
b2mod_b2plot_values.o: ../src/modules/b2mod_b2plot_values.F
	$(FC) -c ../src/modules/b2mod_b2plot_values.F
b2mod_b2plot_xy.o: ../src/modules/b2mod_b2plot_xy.F
	$(FC) -c ../src/modules/b2mod_b2plot_xy.F
b2mod_wall.o: ../src/modules/b2mod_wall.F
	$(FC) -c -DNO_CDF ../src/modules/b2mod_wall.F -I../src/include/
b2mod_diag.o: ../src/modules/b2mod_diag.F
	$(FC) -c ../src/modules/b2mod_diag.F


b2rfcp.o: ../src/input/b2rfcp.F
	$(FC) -c ../src/input/b2rfcp.F -I../src/common/
b2rflb.o: ../src/input/b2rflb.F
	$(FC) -c ../src/input/b2rflb.F
b2rucp.o: ../src/input/b2rucp.F
	$(FC) -c ../src/input/b2rucp.F -I../src/common/
b2rugm.o: ../src/input/b2rugm.F
	$(FC) -c ../src/input/b2rugm.F -I../src/common/
b2rups.o: ../src/input/b2rups.F
	$(FC) -c ../src/input/b2rups.F
b2rurc.o: ../src/input/b2rurc.F
	$(FC) -c ../src/input/b2rurc.F -I../src/common/
b2ruzd.o: ../src/input/b2rusr.F
	$(FC) -c ../src/input/b2ruzd.F -I../src/common/

b2spcx.o: ../src/sources/b2spcx.F
	$(FC) -c ../src/sources/b2spcx.F -I../src/common/
b2spel.o: ../src/sources/b2spel.F
	$(FC) -c ../src/sources/b2spel.F -I../src/common/
b2sqcx.o: ../src/sources/b2sqcx.F
	$(FC) -c ../src/sources/b2sqcx.F
b2sqel.o: ../src/sources/b2sqel.F
	$(FC) -c ../src/sources/b2sqel.F
eirene_f30f31.o: ../src/sources/eirene_f30f31.F
	$(FC) -c ../src/sources/eirene_f30f31.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
eirene_mc_init.o: ../src/sources/eirene_mc_init.F
	$(FC) -c ../src/sources/eirene_mc_init.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
heatdiff1D.o: ../src/sources/heatdiff1D.F
	$(FC) -c -DNO_NAG ../src/sources/heatdiff1D.F -I../src/common/
 
b2tlnl.o: ../src/transport/b2tlnl.F
	$(FC) -c ../src/transport/b2tlnl.F -I../src/common/

cfopen.o: ../src/utility/cfopen.F
	$(FC) -c ../src/utility/cfopen.F
cfruch.o: ../src/utility/cfruch.F
	$(FC) -c ../src/utility/cfruch.F
cfruin.o: ../src/utility/cfruin.F
	$(FC) -c ../src/utility/cfruin.F
cfrure.o: ../src/utility/cfrure.F
	$(FC) -c ../src/utility/cfrure.F
cfvers.o: ../src/utility/cfvers.F
	$(FC) -c ../src/utility/cfvers.F
cfwuin.o: ../src/utility/cfwuin.F
	$(FC) -c ../src/utility/cfwuin.F
cfwure.o: ../src/utility/cfwure.F
	$(FC) -c ../src/utility/cfwure.F
chcase.o: ../src/utility/chcase.F
	$(FC) -c ../src/utility/chcase.F
daytim.o: ../src/utility/daytim.F
	$(FC) -c ../src/utility/daytim.F
gfsub2.o: ../src/utility/gfsub2.F
	$(FC) -c ../src/utility/gfsub2.F
get_jsep.o: ../src/utility/get_jsep.F
	$(FC) -c ../src/utility/get_jsep.F
ifill.o: ../src/utility/ifill.F
	$(FC) -c ../src/utility/ifill.F
intfaceh.o: ../src/utility/intfaceh.F
	$(FC) -c ../src/utility/intfaceh.F
intfacev.o: ../src/utility/intfacev.F
	$(FC) -c ../src/utility/intfacev.F
ipgeti.o: ../src/utility/ipgeti.F
	$(FC) -c ../src/utility/ipgeti.F
ipgetr.o: ../src/utility/ipgetr.F
	$(FC) -c ../src/utility/ipgetr.F
ipmain.o: ../src/utility/ipmain.F
	$(FC) -c ../src/utility/ipmain.F
jobnam.o: ../src/utility/jobnam.F
	$(FC) -c ../src/utility/jobnam.F
len_of_digits.o: ../src/utility/len_of_digits.F
	$(FC) -c ../src/utility/len_of_digits.F
lwimai.o: ../src/utility/lwimai.F
	$(FC) -c ../src/utility/lwimai.F -I../src/common/
lwmain.o: ../src/utility/lwmain.F
	$(FC) -c ../src/utility/lwmain.F -I../src/common/
machsfr.o: ../src/utility/machsfr.F
	$(FC) -c ../src/utility/machsfr.F
mstep.o: ../src/utility/mstep.F
	$(FC) -c ../src/utility/mstep.F
my_intcor.o: ../src/utility/my_intcor.F
	$(FC) -c ../src/utility/my_intcor.F
my_out.o: ../src/utility/my_out.F
	$(FC) -c ../src/utility/my_out.F
myblas.o: ../src/utility/myblas.F
	$(FC) -c ../src/utility/myblas.F
nagsubst.o: ../src/utility/nagsubst.F
	$(FC) -c -DNO_NAG ../src/utility/nagsubst.F
open_file.o: ../src/utility/open_file.F
	$(FC) -c ../src/utility/open_file.F
prgend.o: ../src/utility/prgend.F
	$(FC) -c ../src/utility/prgend.F
prgini.o: ../src/utility/prgini.F
	$(FC) -c ../src/utility/prgini.F -I../src/common/
prvrt.o: ../src/utility/prvrt.F
	$(FC) -c ../src/utility/prvrt.F
prvrti.o: ../src/utility/prvrti.F
	$(FC) -c ../src/utility/prvrti.F
rdneutrs.o: ../src/utility/rdneutrs.F
	$(FC) -c ../src/utility/rdneutrs.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
samax.o: ../src/utility/samax.F
	$(FC) -c ../src/utility/samax.F
sfill.o: ../src/utility/sfill.F
	$(FC) -c ../src/utility/sfill.F
smax.o: ../src/utility/smax.F
	$(FC) -c ../src/utility/smax.F
smin.o: ../src/utility/smin.F
	$(FC) -c ../src/utility/smin.F
sortda.o: ../src/utility/sortda.F
	$(FC) -c ../src/utility/sortda.F
ssum.o: ../src/utility/ssum.F
	$(FC) -c ../src/utility/ssum.F
streql.o: ../src/utility/streql.F
	$(FC) -c ../src/utility/streql.F
subsys.o: ../src/utility/subsys.F
	$(FC) -c ../src/utility/subsys.F
sysend.o: ../src/utility/sysend.F
	$(FC) -c ../src/utility/sysend.F
sysini.o: ../src/utility/sysini.F
	$(FC) -c ../src/utility/sysini.F
usrnam.o: ../src/utility/usrnam.F
	$(FC) -c ../src/utility/usrnam.F
xerrab.o: ../src/utility/xerrab.F
	$(FC) -c ../src/utility/xerrab.F
xerset.o: ../src/utility/xerset.F
	$(FC) -c ../src/utility/xerset.F
xertst.o: ../src/utility/xertst.F
	$(FC) -c ../src/utility/xertst.F

b2xgbs.o: ../src/b2aux/b2xgbs.F
	$(FC) -c ../src/b2aux/b2xgbs.F
b2xpne.o: ../src/b2aux/b2xpne.F
	$(FC) -c ../src/b2aux/b2xpne.F
b2xpni.o: ../src/b2aux/b2xpni.F
	$(FC) -c ../src/b2aux/b2xpni.F
b2xvcp.o: ../src/b2aux/b2xvcp.F
	$(FC) -c ../src/b2aux/b2xvcp.F -I../src/common/
b2xvff.o: ../src/b2aux/b2xvff.F
	$(FC) -c ../src/b2aux/b2xvff.F
b2xvfx.o: ../src/b2aux/b2xvfx.F
	$(FC) -c ../src/b2aux/b2xvfx.F
b2xvfy.o: ../src/b2aux/b2xvfy.F
	$(FC) -c ../src/b2aux/b2xvfy.F
b2xvgm.o: ../src/b2aux/b2xvgm.F
	$(FC) -c ../src/b2aux/b2xvgm.F -I../src/common/
b2xvps.o: ../src/b2aux/b2xvps.F
	$(FC) -c ../src/b2aux/b2xvps.F
b2xvsg.o: ../src/b2aux/b2xvsg.F
	$(FC) -c ../src/b2aux/b2xvsg.F
b2xxgs.o: ../src/b2aux/b2xxgs.F
	$(FC) -c ../src/b2aux/b2xxgs.F
b2xxid.o: ../src/b2aux/b2xxid.F
	$(FC) -c ../src/b2aux/b2xxid.F

additional_elements.o: ../src/b2plot/additional_elements.F
	$(FC) -c ../src/b2plot/additional_elements.F -I../src/common/
ank_interface.o: ../src/b2plot/ank_interface.F
	$(FC) -c ../src/b2plot/ank_interface.F -I../src/common/
avey.o: ../src/b2plot/avey.F
	$(FC) -c ../src/b2plot/avey.F
b2paux.o: ../src/b2plot/b2paux.F
	$(FC) -c ../src/b2plot/b2paux.F
b2pinp.o: ../src/b2plot/b2pinp.F
	$(FC) -c ../src/b2plot/b2pinp.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
b2pldr.o: ../src/b2plot/b2pldr.F
	$(FC) -c ../src/b2plot/b2pldr.F
b2plot.o: ../src/b2plot/b2plot.F
	$(FC) -c -DNO_CDF ../src/b2plot/b2plot.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
b2pnwl.o: ../src/b2plot/b2pnwl.F
	$(FC) -c ../src/b2plot/b2pnwl.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
b2prdl.o: ../src/b2plot/b2prdl.F
	$(FC) -c ../src/b2plot/b2prdl.F
b2psh.o: ../src/b2plot/b2psh.F
	$(FC) -c ../src/b2plot/b2psh.F
b2ptrdl.o: ../src/b2plot/b2ptrdl.F
	$(FC) -c ../src/b2plot/b2ptrdl.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
b2ptrg.o: ../src/b2plot/b2ptrg.F
	$(FC) -c ../src/b2plot/b2ptrg.F
b2ptrgl.o: ../src/b2plot/b2ptrgl.F
	$(FC) -c ../src/b2plot/b2ptrgl.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
b2pwinm.o: ../src/b2plot/b2pwinm.F
	$(FC) -c ../src/b2plot/b2pwinm.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
b2pwinp.o: ../src/b2plot/b2pwinp.F
	$(FC) -c ../src/b2plot/b2pwinp.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
b2pwld.o: ../src/b2plot/b2pwld.F
	$(FC) -c ../src/b2plot/b2pwld.F
b2pwlprp.o: ../src/b2plot/b2pwlprp.F
	$(FC) -c ../src/b2plot/b2pwlprp.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
b2pwrdld.o: ../src/b2plot/b2pwrdld.F
	$(FC) -c ../src/b2plot/b2pwrdld.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
b2pwsort.o: ../src/b2plot/b2pwsort.F
	$(FC) -c ../src/b2plot/b2pwsort.F -I../src/common/ -I../src/common/COUPLE/ -I../src/include/
caltraj.o: ../src/b2plot/caltraj.F
	$(FC) -c -DNO_NAG ../src/b2plot/caltraj.F
ceil.o: ../src/b2plot/ceil.F
	$(FC) -c ../src/b2plot/ceil.F
chord.o: ../src/b2plot/chord.F
	$(FC) -c ../src/b2plot/chord.F
clear_extralabel.o: ../src/b2plot/clear_extralabel.F
	$(FC) -c ../src/b2plot/clear_extralabel.F -I../src/common/
clip_box.o: ../src/b2plot/clip_box.F
	$(FC) -c ../src/b2plot/clip_box.F
clip_line.o: ../src/b2plot/clip_line.F
	$(FC) -c ../src/b2plot/clip_line.F
clvd.o: ../src/b2plot/clvd.F
	$(FC) -c ../src/b2plot/clvd.F -I../src/common/
clvl.o: ../src/b2plot/clvl.F
	$(FC) -c ../src/b2plot/clvl.F -I../src/common/
clvn.o: ../src/b2plot/clvn.F
	$(FC) -c ../src/b2plot/clvn.F -I../src/common/
clvp.o: ../src/b2plot/clvp.F
	$(FC) -c ../src/b2plot/clvp.F -I../src/common/
clvs.o: ../src/b2plot/clvs.F
	$(FC) -c ../src/b2plot/clvs.F -I../src/common/
color.o: ../src/b2plot/color.F
	$(FC) -c ../src/b2plot/color.F -I../src/common/
cont.o: ../src/b2plot/cont.F
	$(FC) -c ../src/b2plot/cont.F
contourplot.o: ../src/b2plot/contourplot.F
	$(FC) -c ../src/b2plot/contourplot.F -I../src/common/
copy.o: ../src/b2plot/copy.F
	$(FC) -c ../src/b2plot/copy.F -I../src/common/
coronal.o: ../src/b2plot/coronal.F
	$(FC) -c ../src/b2plot/coronal.F
dhdata.o: ../src/b2plot/dhdata.F
	$(FC) -c ../src/b2plot/dhdata.F
dropstackbottom.o: ../src/b2plot/dropstackbottom.F
	$(FC) -c ../src/b2plot/dropstackbottom.F
end_page.o: ../src/b2plot/end_page.F
	$(FC) -c -DNCAR4 ../src/b2plot/end_page.F -I../src/common/
extend.o: ../src/b2plot/extend.F
	$(FC) -c ../src/b2plot/extend.F
fill.o: ../src/b2plot/fill.F
	$(FC) -c ../src/b2plot/fill.F -I../src/common/
floor.o: ../src/b2plot/floor.F
	$(FC) -c ../src/b2plot/floor.F
flux.o: ../src/b2plot/flux.F
	$(FC) -c ../src/b2plot/flux.F -I../src/common/
flx.o: ../src/b2plot/flx.F
	$(FC) -c ../src/b2plot/flx.F
fly.o: ../src/b2plot/fly.F
	$(FC) -c ../src/b2plot/fly.F
font.o: ../src/b2plot/font.F
	$(FC) -c ../src/b2plot/font.F
fx.o: ../src/b2plot/fx.F
	$(FC) -c ../src/b2plot/fx.F
fy.o: ../src/b2plot/fy.F
	$(FC) -c ../src/b2plot/fy.F
getint.o: ../src/b2plot/getint.F
	$(FC) -c ../src/b2plot/getint.F
gridon.o: ../src/b2plot/gridon.F
	$(FC) -c ../src/b2plot/gridon.F -I../src/common/
gslwsc_d.o: ../src/b2plot/gslwsc_d.F
	$(FC) -c ../src/b2plot/gslwsc_d.F
init_page.o: ../src/b2plot/init_page.F
	$(FC) -c ../src/b2plot/init_page.F -I../src/common/
init_wall.o: ../src/b2plot/init_wall.F
	$(FC) -c ../src/b2plot/init_wall.F -I../src/include/
initialize_plotting.o: ../src/b2plot/initialize_plotting.F
	$(FC) -c -DNCAR4 ../src/b2plot/initialize_plotting.F -I../src/common/
ipolate.o: ../src/b2plot/ipolate.F
	$(FC) -c ../src/b2plot/ipolate.F -I../src/common/
ipop.o: ../src/b2plot/ipop.F
	$(FC) -c ../src/b2plot/ipop.F
ipush.o: ../src/b2plot/ipush.F
	$(FC) -c ../src/b2plot/ipush.F
isadigit.o: ../src/b2plot/isadigit.F
	$(FC) -c ../src/b2plot/isadigit.F
isnegative.o: ../src/b2plot/isnegative.F
	$(FC) -c ../src/b2plot/isnegative.F
isplusorminus.o: ../src/b2plot/isplusorminus.F
	$(FC) -c ../src/b2plot/isplusorminus.F
iswap.o: ../src/b2plot/iswap.F
	$(FC) -c ../src/b2plot/iswap.F
labelon.o: ../src/b2plot/labelon.F
	$(FC) -c ../src/b2plot/labelon.F -I../src/common/
lblalongcl.o: ../src/b2plot/lblalongcl.F
	$(FC) -c ../src/b2plot/lblalongcl.F -I../src/common/
lblbartraj.o: ../src/b2plot/lblbartraj.F
	$(FC) -c ../src/b2plot/lblbartraj.F -I../src/common/
line.o: ../src/b2plot/line.F
	$(FC) -c ../src/b2plot/line.F
lined_d.o: ../src/b2plot/lined_d.F
	$(FC) -c ../src/b2plot/lined_d.F -I../src/common/
linewidth.o: ../src/b2plot/linewidth.F
	$(FC) -c ../src/b2plot/linewidth.F -I../src/common/
lower_case.o: ../src/b2plot/lower_case.F
	$(FC) -c ../src/b2plot/lower_case.F
map.o: ../src/b2plot/map.F
	$(FC) -c ../src/b2plot/map.F
map_chords.o: ../src/b2plot/map_chords.F
	$(FC) -c ../src/b2plot/map_chords.F
mapx.o: ../src/b2plot/mapx.F
	$(FC) -c ../src/b2plot/mapx.F
mapy.o: ../src/b2plot/mapy.F
	$(FC) -c ../src/b2plot/mapy.F
mesh.o: ../src/b2plot/mesh.F
	$(FC) -c ../src/b2plot/mesh.F
meshit.o: ../src/b2plot/meshit.F
	$(FC) -c ../src/b2plot/meshit.F
ngread.o: ../src/b2plot/ngread.F
	$(FC) -c ../src/b2plot/ngread.F -I../src/common/COUPLE/ -I../src/common/ -I../src/include/
nice.o: ../src/b2plot/nice.F
	$(FC) -c ../src/b2plot/nice.F
nofill.o: ../src/b2plot/nofill.F
	$(FC) -c ../src/b2plot/nofill.F -I../src/common/
outlineon.o: ../src/b2plot/outlineon.F
	$(FC) -c ../src/b2plot/outlineon.F -I../src/common/
page.o: ../src/b2plot/page.F
	$(FC) -c ../src/b2plot/page.F -I../src/common/
plot_chords.o: ../src/b2plot/plot_chords.F
	$(FC) -c ../src/b2plot/plot_chords.F
plot_wall_elements.o: ../src/b2plot/plot_wall_elements.F
	$(FC) -c ../src/b2plot/plot_wall_elements.F -I../src/common/
plotcontour.o: ../src/b2plot/plotcontour.F
	$(FC) -c ../src/b2plot/plotcontour.F -I../src/common/
plotpatch.o: ../src/b2plot/plotpatch.F
	$(FC) -c ../src/b2plot/plotpatch.F -I../src/common/
plottarget.o: ../src/b2plot/plottarget.F
	$(FC) -c ../src/b2plot/plottarget.F -I../src/common/
pltoutline.o: ../src/b2plot/pltoutline.F
	$(FC) -c ../src/b2plot/pltoutline.F
pltsptrix.o: ../src/b2plot/pltsptrix.F
	$(FC) -c ../src/b2plot/pltsptrix.F
plttarget.o: ../src/b2plot/plttarget.F
	$(FC) -c ../src/b2plot/plttarget.F -I../src/common/
pltvessel.o: ../src/b2plot/pltvessel.F
	$(FC) -c ../src/b2plot/pltvessel.F
pltwall.o: ../src/b2plot/pltwall.F
	$(FC) -c ../src/b2plot/pltwall.F -I../src/common/
pne.o: ../src/b2plot/pne.F
	$(FC) -c ../src/b2plot/pne.F
process_chord_slice.o: ../src/b2plot/process_chord_slice.F
	$(FC) -c ../src/b2plot/process_chord_slice.F
process_chords.o: ../src/b2plot/process_chords.F
	$(FC) -c ../src/b2plot/process_chords.F
process_chords_boxcar.o: ../src/b2plot/process_chords_boxcar.F
	$(FC) -c ../src/b2plot/process_chords_boxcar.F
pswa.o: ../src/b2plot/pswa.F
	$(FC) -c ../src/b2plot/pswa.F
psx.o: ../src/b2plot/psx.F
	$(FC) -c ../src/b2plot/psx.F
psy.o: ../src/b2plot/psy.F
	$(FC) -c ../src/b2plot/psy.F
pvol.o: ../src/b2plot/pvol.F
	$(FC) -c ../src/b2plot/pvol.F
readline.o: ../src/b2plot/readline.F
	$(FC) -c ../src/b2plot/readline.F
rlcomp.o: ../src/b2plot/rlcomp.F
	$(FC) -c ../src/b2plot/rlcomp.F -I../src/common/
rlshdw.o: ../src/b2plot/rlshdw.F
	$(FC) -c ../src/b2plot/rlshdw.F -I../src/common/
rpop.o: ../src/b2plot/rpop.F
	$(FC) -c ../src/b2plot/rpop.F
rpush.o: ../src/b2plot/rpush.F
	$(FC) -c ../src/b2plot/rpush.F
rswap.o: ../src/b2plot/rswap.F
	$(FC) -c ../src/b2plot/rswap.F
separatrixon.o: ../src/b2plot/separatrixon.F
	$(FC) -c ../src/b2plot/separatrixon.F -I../src/common/
set_a4l.o: ../src/b2plot/set_a4l.F
	$(FC) -c ../src/b2plot/set_a4l.F -I../src/common/
set_a4p.o: ../src/b2plot/set_a4p.F
	$(FC) -c ../src/b2plot/set_a4p.F -I../src/common/
set_extralabel.o: ../src/b2plot/set_extralabel.F
	$(FC) -c ../src/b2plot/set_extralabel.F -I../src/common/COUPLE/ -I../src/common/ -I../src/include/
set_globalheader.o: ../src/b2plot/set_globalheader.F
	$(FC) -c ../src/b2plot/set_globalheader.F -I../src/common/
set_midplane_map.o: ../src/b2plot/set_midplane_map.F
	$(FC) -c ../src/b2plot/set_midplane_map.F
set_new_page.o: ../src/b2plot/set_new_page.F
	$(FC) -c ../src/b2plot/set_new_page.F -I../src/common/
set_same_page.o: ../src/b2plot/set_same_page.F
	$(FC) -c ../src/b2plot/set_same_page.F -I../src/common/
set_square.o: ../src/b2plot/set_square.F
	$(FC) -c ../src/b2plot/set_square.F -I../src/common/
shrink_label.o: ../src/b2plot/shrink_label.F
	$(FC) -c ../src/b2plot/shrink_label.F
skipblanks.o: ../src/b2plot/skipblanks.F
	$(FC) -c ../src/b2plot/skipblanks.F -I../src/common/COUPLE/ -I../src/common/ -I../src/include/
species.o: ../src/b2plot/species.F
	$(FC) -c ../src/b2plot/species.F
spect.o: ../src/b2plot/spect.F
	$(FC) -c ../src/b2plot/spect.F -I../src/common/
spop.o: ../src/b2plot/spop.F
	$(FC) -c ../src/b2plot/spop.F
spush.o: ../src/b2plot/spush.F
	$(FC) -c ../src/b2plot/spush.F
sumx.o: ../src/b2plot/sumx.F
	$(FC) -c ../src/b2plot/sumx.F
sumy.o: ../src/b2plot/sumy.F
	$(FC) -c ../src/b2plot/sumy.F
strswap.o: ../src/b2plot/strswap.F
	$(FC) -c ../src/b2plot/strswap.F
surf.o: ../src/b2plot/surf.F
	$(FC) -c ../src/b2plot/surf.F
symbol.o: ../src/b2plot/symbol.F
	$(FC) -c ../src/b2plot/symbol.F
targeton.o: ../src/b2plot/targeton.F
	$(FC) -c ../src/b2plot/targeton.F -I../src/common/
terminate_plotting.o: ../src/b2plot/terminate_plotting.F
	$(FC) -c -DNCAR4 ../src/b2plot/terminate_plotting.F -I../src/common/
trgsumw.o: ../src/b2plot/trgsumw.F
	$(FC) -c ../src/b2plot/trgsumw.F
trgsurf.o: ../src/b2plot/trgsurf.F
	$(FC) -c ../src/b2plot/trgsurf.F
trgw.o: ../src/b2plot/trgw.F
	$(FC) -c ../src/b2plot/trgw.F
trgz.o: ../src/b2plot/trgz.F
	$(FC) -c ../src/b2plot/trgz.F
vesselon.o: ../src/b2plot/vesselon.F
	$(FC) -c ../src/b2plot/vesselon.F -I../src/common/
wallon.o: ../src/b2plot/wallon.F
	$(FC) -c ../src/b2plot/wallon.F -I../src/common/

b2mwti.o: ../src/driver/b2mwti.F
	$(FC) -c -DNO_CDF ../src/driver/b2mwti.F  -I../src/common/ -I../src/include/
#b2pl.o: ../src/driver/b2pl.F
#	$(FC) -c ../src/driver/b2pl.F
b2pl.o: ../../rest/B2.5/src/driver/b2pl.F
	$(FC) -c ../../rest/B2.5/src/driver/b2pl.F

clean:
	rm -f *.o *.mod





