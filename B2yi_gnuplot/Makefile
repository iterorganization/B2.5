FC = gfortran -g
NCARG=~/usr/local/ncarg
MSCL=~/usr/local/mscl/linux.intel

.SUFFIXES: .o .f .f90 .F



OBJ1 =  b2mod_types.o b2mod_rates.o b2mod_version.o b2mod_geo_corner.o b2mod_geo.o b2mod_work.o  b2mod_constants.o b2mod_b2cmpa.o b2mod_indirect.o b2mod_residuals.o b2mod_sources.o b2mod_transport.o b2mod_plasma.o b2mod_ppout.o b2mod_b2cmrc.o b2mod_b2cmfs.o b2mod_b2cmgs.o b2mod_time.o b2mod_b2cmpb.o b2mod_b2cmpt.o
OBJ2 =  cfendf.o cfopen.o cfruch.o cfruin.o cfrure.o chcase.o  ifill.o  myblas.o
OBJ3 =  ipmain.o lwimai.o lwmain.o machsfr.o nagsubst.o  ipgetr.o ipgeti.o smax.o smin.o samax.o ssum.o
OBJ4 =  prvrt.o prvrti.o 
OBJ5 =  scal10.o sfill.o streql.o subsys.o sysend.o sysini.o  xerrab.o xerset.o xertst.o open_file.o cfwure.o cfvers.o
OBJ6 =  b2rugm.o b2ruzd.o b2rups.o b2xvsg.o b2xvgm.o b2xpne.o b2xpni.o b2rurc.o b2spel.o b2sqel.o b2xvfx.o b2xvfy.o b2xvps.o b2xvff.o b2xpfe.o
OBJ7 =  b2yi_gnuplot.o b2wgpi_gnuplot.o b2wgpj_gnuplot.o
OBJ8 =  b2yidr_gnuplot_2.o
OBJ9 =  intfaceh.o intfacev.o b2xgbs.o b2xxid.o b2xxgs.o b2rfcp.o b2rflb.o jobnam.o prgini.o prgend.o daytim.o usrnam.o
INCLUDE = ../src/common

.F.o:
	$(FC) $(FFLAGS) $(INCDIR) $(OPTS) $?
.f.o:
	$(FC) $(FFLAGS) $(INCDIR) $(OPTS) $?



b2yi.exe: $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9)
#	$(FC) *.o -o Release/b2yi.exe  $(MSCL)/libmscl.a $(NCARG)/lib/libncarg.a $(NCARG)/lib/libncarg_gks.a $(NCARG)/lib/libncarg_c.a $(NCARG)/lib/libngmath.a $(NCARG)/lib/libcgm.a $(NCARG)/lib/libncarg_ras.a -lX11 -lm 
	$(FC) *.o -o Release/b2yi.exe  



b2mod_types.o: ../src/modules/b2mod_types.F
	$(FC) -c ../src/modules/b2mod_types.F
b2mod_time.o: ../src/modules/b2mod_time.F
	$(FC) -c ../src/modules/b2mod_time.F
b2mod_b2cmfs.o: ../src/modules/b2mod_b2cmfs.F
	$(FC) -c ../src/modules/b2mod_b2cmfs.F
b2mod_b2cmrc.o: ../src/modules/b2mod_b2cmrc.F
	$(FC) -c ../src/modules/b2mod_b2cmrc.F
b2mod_b2cmgs.o: ../src/modules/b2mod_b2cmgs.F
	$(FC) -c ../src/modules/b2mod_b2cmgs.F

b2mod_rates.o: ../src/modules/b2mod_rates.F
	$(FC) -c ../src/modules/b2mod_rates.F

b2mod_version.o: ../src/modules/b2mod_version.F
	$(FC) -c ../src/modules/b2mod_version.F

b2mod_geo_corner.o: ../src/modules/b2mod_geo_corner.F
	$(FC) -c ../src/modules/b2mod_geo_corner.F

b2mod_geo.o: ../src/modules/b2mod_geo.F
	$(FC) -c ../src/modules/b2mod_geo.F
b2mod_work.o: ../src/modules/b2mod_work.F
	$(FC) -c ../src/modules/b2mod_work.F
b2mod_plasma.o: ../src/modules/b2mod_plasma.F
	$(FC) -c ../src/modules/b2mod_plasma.F
b2mod_constants.o: ../src/modules/b2mod_constants.F
	$(FC) -c ../src/modules/b2mod_constants.F

b2mod_b2cmpa.o: ../src/modules/b2mod_b2cmpa.F
	$(FC) -c ../src/modules/b2mod_b2cmpa.F
b2mod_b2cmpb.o: ../src/modules/b2mod_b2cmpb.F
	$(FC) -c ../src/modules/b2mod_b2cmpb.F
b2mod_b2cmpt.o: ../src/modules/b2mod_b2cmpt.F
	$(FC) -c ../src/modules/b2mod_b2cmpt.F

b2mod_indirect.o: ../src/modules/b2mod_indirect.F
	$(FC) -c ../src/modules/b2mod_indirect.F
b2mod_residuals.o: ../src/modules/b2mod_residuals.F
	$(FC) -c ../src/modules/b2mod_residuals.F
b2mod_sources.o: ../src/modules/b2mod_sources.F
	$(FC) -c ../src/modules/b2mod_sources.F
b2mod_transport.o: ../src/modules/b2mod_transport.F
	$(FC) -c ../src/modules/b2mod_transport.F
b2mod_ppout.o: ../src/modules/b2mod_ppout.F
	$(FC) -c ../src/modules/b2mod_ppout.F

myblas.o: ../src/utility/myblas.F
	$(FC) -c ../src/utility/myblas.F
cfvers.o: ../src/utility/cfvers.F
	$(FC) -c ../src/utility/cfvers.F
cfwure.o: ../src/utility/cfwure.F
	$(FC) -c ../src/utility/cfwure.F
cfendf.o: ../src/utility/cfendf.F
	$(FC) -c ../src/utility/cfendf.F
cfopen.o: ../src/utility/cfopen.F
	$(FC) -c ../src/utility/cfopen.F
cfruch.o: ../src/utility/cfruch.F
	$(FC) -c ../src/utility/cfruch.F
cfruin.o: ../src/utility/cfruin.F
	$(FC) -c ../src/utility/cfruin.F
cfrure.o: ../src/utility/cfrure.F
	$(FC) -c ../src/utility/cfrure.F
chcase.o: ../src/utility/chcase.F
	$(FC) -c ../src/utility/chcase.F
gnopen.o: ../src/utility/gnopen.F
	$(FC) -c ../src/utility/gnopen.F
gwsrep.o: ../src/utility/gwsrep.F
	$(FC) -c ../src/utility/gwsrep.F
ifill.o: ../src/utility/ifill.F
	$(FC) -c ../src/utility/ifill.F
ipmain.o: ../src/utility/ipmain.F
	$(FC) -c ../src/utility/ipmain.F
lwimai.o: ../src/utility/lwimai.F
	$(FC) -c ../src/utility/lwimai.F -I../src/common/
lwmain.o: ../src/utility/lwmain.F
	$(FC) -c ../src/utility/lwmain.F -I../src/common/
machsfr.o: ../src/utility/machsfr.F
	$(FC) -c ../src/utility/machsfr.F
nagsubst.o: ../src/utility/nagsubst.F
	$(FC) -c -DNO_NAG ../src/utility/nagsubst.F
ncardummy.o: ../src/utility/ncardummy.F
	$(FC) -c ../src/utility/ncardummy.F
pagdef.o: ../src/utility/pagdef.F
	$(FC) -c ../src/utility/pagdef.F
pagend.o: ../src/utility/pagend.F
	$(FC) -c ../src/utility/pagend.F
pagini.o: ../src/utility/pagini.F
	$(FC) -c ../src/utility/pagini.F
prgend.o: ../src/utility/prgend.F
	$(FC) -c ../src/utility/prgend.F
prgini.o: ../src/utility/prgini.F
	$(FC) -c ../src/utility/prgini.F -I../src/common/
prvrt.o: ../src/utility/prvrt.F
	$(FC) -c ../src/utility/prvrt.F
prvrti.o: ../src/utility/prvrti.F
	$(FC) -c ../src/utility/prvrti.F
scal10.o: ../src/utility/scal10.F
	$(FC) -c ../src/utility/scal10.F
sfill.o: ../src/utility/sfill.F
	$(FC) -c ../src/utility/sfill.F
streql.o: ../src/utility/streql.F
	$(FC) -c ../src/utility/streql.F
subsys.o: ../src/utility/subsys.F
	$(FC) -c ../src/utility/subsys.F
sysend.o: ../src/utility/sysend.F
	$(FC) -c ../src/utility/sysend.F
sysini.o: ../src/utility/sysini.F
	$(FC) -c ../src/utility/sysini.F
xerrab.o: ../src/utility/xerrab.F
	$(FC) -c ../src/utility/xerrab.F
xerset.o: ../src/utility/xerset.F
	$(FC) -c ../src/utility/xerset.F
xertst.o: ../src/utility/xertst.F
	$(FC) -c ../src/utility/xertst.F
open_file.o: ../src/utility/open_file.F
	$(FC) -c ../src/utility/open_file.F
samax.o: ../src/utility/samax.F
	$(FC) -c ../src/utility/samax.F
smax.o: ../src/utility/smax.F
	$(FC) -c ../src/utility/smax.F
ssum.o: ../src/utility/ssum.F
	$(FC) -c ../src/utility/ssum.F
smin.o: ../src/utility/smin.F
	$(FC) -c ../src/utility/smin.F
ipgetr.o: ../src/utility/ipgetr.F
	$(FC) -c ../src/utility/ipgetr.F
ipgeti.o: ../src/utility/ipgeti.F
	$(FC) -c ../src/utility/ipgeti.F
gslwsc_d.o: ../src/b2plot/gslwsc_d.F
	$(FC) -c ../src/b2plot/gslwsc_d.F


b2rugm.o: ../src/input/b2rugm.F
	$(FC) -c ../src/input/b2rugm.F
b2ruzd.o: ../src/input/b2ruzd.F
	$(FC) -c ../src/input/b2ruzd.F
b2rups.o: ../src/input/b2rups.F
	$(FC) -c ../src/input/b2rups.F
b2rurc.o: ../src/input/b2rurc.F
	$(FC) -c ../src/input/b2rurc.F
b2rfcp.o: ../src/input/b2rfcp.F
	$(FC) -c ../src/input/b2rfcp.F
b2rflb.o: ../src/input/b2rflb.F
	$(FC) -c ../src/input/b2rflb.F


b2spel.o: ../src/sources/b2spel.F
	$(FC) -c ../src/sources/b2spel.F
b2sqel.o: ../src/sources/b2sqel.F
	$(FC) -c ../src/sources/b2sqel.F

b2xvsg.o: ../src/b2aux/b2xvsg.F
	$(FC) -c ../src/b2aux/b2xvsg.F
b2xvgm.o: ../src/b2aux/b2xvgm.F
	$(FC) -c ../src/b2aux/b2xvgm.F
b2xpne.o: ../src/b2aux/b2xpne.F
	$(FC) -c ../src/b2aux/b2xpne.F
b2xpni.o: ../src/b2aux/b2xpni.F
	$(FC) -c ../src/b2aux/b2xpni.F
b2xpfe.o: ../src/b2aux/b2xpfe.F
	$(FC) -c ../src/b2aux/b2xpfe.F
b2xvfx.o: ../src/b2aux/b2xvfx.F
	$(FC) -c ../src/b2aux/b2xvfx.F
b2xvfy.o: ../src/b2aux/b2xvfy.F
	$(FC) -c ../src/b2aux/b2xvfy.F
b2xvps.o: ../src/b2aux/b2xvps.F
	$(FC) -c ../src/b2aux/b2xvps.F
b2xvff.o: ../src/b2aux/b2xvff.F
	$(FC) -c ../src/b2aux/b2xvff.F
b2xgbs.o: ../src/b2aux/b2xgbs.F
	$(FC) -c ../src/b2aux/b2xgbs.F
b2xxid.o: ../src/b2aux/b2xxid.F
	$(FC) -c ../src/b2aux/b2xxid.F
b2xxgs.o: ../src/b2aux/b2xxgs.F
	$(FC) -c ../src/b2aux/b2xxgs.F

b2wgpi_gnuplot.o: ../src/output/b2wgpi_gnuplot.F
	$(FC) -c ../src/output/b2wgpi_gnuplot.F
b2wgpj_gnuplot.o: ../src/output/b2wgpj_gnuplot.F
	$(FC) -c ../src/output/b2wgpj_gnuplot.F



b2yi_gnuplot.o: ../src/postprocessing/b2yi_gnuplot.F
	$(FC) -c ../src/postprocessing/b2yi_gnuplot.F
b2yidr_gnuplot_2.o: ../src/postprocessing/b2yidr_gnuplot_2.F
	$(FC) -c ../src/postprocessing/b2yidr_gnuplot_2.F

intfaceh.o: ../src/utility/intfaceh.F
	$(FC) -c ../src/utility/intfaceh.F
intfacev.o: ../src/utility/intfacev.F
	$(FC) -c ../src/utility/intfacev.F

jobnam.o: ../src/utility/jobnam.F
	$(FC) -c ../src/utility/jobnam.F
usrnam.o: ../src/utility/usrnam.F
	$(FC) -c ../src/utility/usrnam.F
daytim.o: ../src/utility/daytim.F
	$(FC) -c ../src/utility/daytim.F

prgini.o: ../src/utility/prgini.F
	$(FC) -c ../src/utility/prgini.F -I$(INCLUDE) 
prgend.o: ../src/utility/prgend.F
	$(FC) -c ../src/utility/prgend.F


clean:
	rm -f *.o *.mod





