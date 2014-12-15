FC = ifort -g

OBJ1 =  b2mod_types.o  b2mod_constants.o b2mod_version.o b2mod_b2cmfs.o b2mod_b2cmgs.o b2mod_b2cmpa.o b2mod_geo_corner.o b2mod_geo.o b2mod_geo2.o      b2mod_indirect.o 
OBJ2 =  b2mod_plasma.o b2mod_rates.o     b2mod_residuals.o b2mod_sources.o b2mod_transport.o b2mod_work.o
OBJ3 =  cfopen.o cfrure.o cfvers.o cfwuch.o cfwuin.o cfwure.o chcase.o daytim.o    esseq.o    eval.o ifill.o cfruin.o cfruch.o
OBJ4 =  intfaceh.o intfacev.o ipgeti.o ipgetr.o ipmain.o jobnam.o lwimai.o lwmain.o nagsubst.o open_file.o  
OBJ5 =  prgend.o prgini.o prvrt.o  prvrti.o ratio.o  samax.o  sfill.o     smax.o     smin.o
OBJ6 =  ssum.o   streql.o strmas.o subsys.o sysend.o sysini.o xerrab.o  xerset.o   xertst.o 
OBJ7 =  usrnam.o 
OBJ8 =  b2xgbs.o b2xvff.o b2xvfx.o b2xvfy.o b2xvgm.o b2xvsg.o b2xxgs.o b2xxid.o
OBJ9 =  b2agbb.o b2agcf.o b2agfs.o b2aggs.o b2agmt.o b2agmx.o b2agsp.o b2agx0.o b2fgmtry_moddr.o 
OBJ10 =  b2wfgm.o b2wugm.o b2rugm.o
OBJ11 =  b2fgmtry_mod.o

b2fgmtry_mod.exe: $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9) $(OBJ10) $(OBJ11)
	$(FC) *.o -o Release/b2fgmtry_mod.exe
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
b2mod_geo_corner.o: ../src/modules/b2mod_geo_corner.F
	$(FC) -c ../src/modules/b2mod_geo_corner.F
b2mod_geo.o: ../src/modules/b2mod_geo.F
	$(FC) -c ../src/modules/b2mod_geo.F
b2mod_geo2.o: ../src/modules/b2mod_geo2.F
	$(FC) -c ../src/modules/b2mod_geo2.F
b2mod_indirect.o: ../src/modules/b2mod_indirect.F
	$(FC) -c ../src/modules/b2mod_indirect.F
b2mod_plasma.o: ../src/modules/b2mod_plasma.F
	$(FC) -c ../src/modules/b2mod_plasma.F
b2mod_rates.o: ../src/modules/b2mod_rates.F
	$(FC) -c ../src/modules/b2mod_rates.F
b2mod_residuals.o: ../src/modules/b2mod_residuals.F
	$(FC) -c ../src/modules/b2mod_residuals.F
b2mod_sources.o: ../src/modules/b2mod_sources.F
	$(FC) -c ../src/modules/b2mod_sources.F
b2mod_transport.o: ../src/modules/b2mod_transport.F
	$(FC) -c ../src/modules/b2mod_transport.F
b2mod_work.o: ../src/modules/b2mod_work.F
	$(FC) -c ../src/modules/b2mod_work.F

cfopen.o: ../src/utility/cfopen.F
	$(FC) -c ../src/utility/cfopen.F
cfrure.o: ../src/utility/cfrure.F
	$(FC) -c ../src/utility/cfrure.F
cfvers.o: ../src/utility/cfvers.F
	$(FC) -c ../src/utility/cfvers.F
cfwuch.o: ../src/utility/cfwuch.F
	$(FC) -c ../src/utility/cfwuch.F
cfruch.o: ../src/utility/cfruch.F
	$(FC) -c ../src/utility/cfruch.F
cfwuin.o: ../src/utility/cfwuin.F
	$(FC) -c ../src/utility/cfwuin.F
cfruin.o: ../src/utility/cfruin.F
	$(FC) -c ../src/utility/cfruin.F
cfwure.o: ../src/utility/cfwure.F
	$(FC) -c ../src/utility/cfwure.F
chcase.o: ../src/utility/chcase.F
	$(FC) -c ../src/utility/chcase.F
daytim.o: ../src/utility/daytim.F
	$(FC) -c ../src/utility/daytim.F
esseq.o: ../src/utility/esseq.F
	$(FC) -c ../src/utility/esseq.F
eval.o: ../src/utility/eval.F
	$(FC) -c ../src/utility/eval.F
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
lwimai.o: ../src/utility/lwimai.F
	$(FC) -c ../src/utility/lwimai.F -I../src/common/
lwmain.o: ../src/utility/lwmain.F
	$(FC) -c ../src/utility/lwmain.F -I../src/common/
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
ratio.o: ../src/utility/ratio.F
	$(FC) -c ../src/utility/ratio.F
samax.o: ../src/utility/samax.F
	$(FC) -c ../src/utility/samax.F
sfill.o: ../src/utility/sfill.F
	$(FC) -c ../src/utility/sfill.F
smax.o: ../src/utility/smax.F
	$(FC) -c ../src/utility/smax.F
smin.o: ../src/utility/smin.F
	$(FC) -c ../src/utility/smin.F
ssum.o: ../src/utility/ssum.F
	$(FC) -c ../src/utility/ssum.F
streql.o: ../src/utility/streql.F
	$(FC) -c ../src/utility/streql.F
strmas.o: ../src/utility/strmas.F
	$(FC) -c ../src/utility/strmas.F
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

usrnam.o: ../src/utility/usrnam.F
	$(FC) -c ../src/utility/usrnam.F

b2xgbs.o: ../src/b2aux/b2xgbs.F
	$(FC) -c ../src/b2aux/b2xgbs.F
b2xvff.o: ../src/b2aux/b2xvff.F
	$(FC) -c ../src/b2aux/b2xvff.F
b2xvfx.o: ../src/b2aux/b2xvfx.F
	$(FC) -c ../src/b2aux/b2xvfx.F
b2xvfy.o: ../src/b2aux/b2xvfy.F
	$(FC) -c ../src/b2aux/b2xvfy.F
b2xvgm.o: ../src/b2aux/b2xvgm.F
	$(FC) -c ../src/b2aux/b2xvgm.F -I../src/common/
b2xvsg.o: ../src/b2aux/b2xvsg.F
	$(FC) -c ../src/b2aux/b2xvsg.F
b2xxgs.o: ../src/b2aux/b2xxgs.F
	$(FC) -c ../src/b2aux/b2xxgs.F
b2xxid.o: ../src/b2aux/b2xxid.F
	$(FC) -c ../src/b2aux/b2xxid.F

b2wfgm.o: ../src/output/b2wfgm.F
	$(FC) -c ../src/output/b2wfgm.F -I../src/common/
b2wugm.o: ../src/output/b2wugm.F
	$(FC) -c ../src/output/b2wugm.F -I../src/common/
b2rugm.o: ../src/input/b2rugm.F
	$(FC) -c ../src/input/b2rugm.F -I../src/common/

b2agbb.o: ../src/preprocessing/b2agbb.F
	$(FC) -c ../src/preprocessing/b2agbb.F
b2agcf.o: ../src/preprocessing/b2agcf.F
	$(FC) -c ../src/preprocessing/b2agcf.F
b2agfs.o: ../src/preprocessing/b2agfs.F
	$(FC) -c ../src/preprocessing/b2agfs.F -I../src/common/
b2aggs.o: ../src/preprocessing/b2aggs.F
	$(FC) -c ../src/preprocessing/b2aggs.F -I../src/common/
b2agmt.o: ../src/preprocessing/b2agmt.F
	$(FC) -c ../src/preprocessing/b2agmt.F
b2agmx.o: ../src/preprocessing/b2agmx.F
	$(FC) -c ../src/preprocessing/b2agmx.F
b2agsp.o: ../src/preprocessing/b2agmx.F
	$(FC) -c ../src/preprocessing/b2agsp.F -I../src/common/
b2agx0.o: ../src/preprocessing/b2agsp.F
	$(FC) -c ../src/preprocessing/b2agx0.F
b2fgmtry_moddr.o: b2fgmtry_moddr.F
	$(FC) -c b2fgmtry_moddr.F -I../src/common/
b2fgmtry_mod.o: b2fgmtry_mod.F
	$(FC) -c b2fgmtry_mod.F

clean:
	rm -f *.o *.mod





