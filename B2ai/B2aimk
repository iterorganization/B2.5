FC = ifort

OBJ1  =  b2mod_types.o  b2mod_constants.o b2mod_b2cmfs.o b2mod_b2cmgs.o b2mod_b2cmpa.o b2mod_version.o b2mod_geo_corner.o b2mod_geo.o b2mod_indirect.o b2mod_plasma.o b2mod_rates.o b2mod_time.o
OBJ2  =  cfopen.o cfruin.o cfruch.o cfrure.o cfvers.o cfwuch.o cfwuin.o cfwure.o chcase.o daytim.o ifill.o  
OBJ3  =  intfaceh.o intfacev.o ipgeti.o ipmain.o jobnam.o lwimai.o lwmain.o nagsubst.o  
OBJ4  =  open_file.o prgend.o prgini.o prvrt.o  prvrti.o samax.o  sfill.o smax.o smin.o
OBJ5  =  streql.o strmas.o subsys.o sysend.o sysini.o xerrab.o  xerset.o   xertst.o 
OBJ6  =  usrnam.o 
OBJ7  =  b2xgbs.o b2xpne.o b2xvff.o b2xvfx.o b2xvfy.o b2xvgm.o b2xvps.o b2xvsg.o b2xxgs.o b2xxid.o
OBJ8  =  b2aix0.o b2aidr.o 
OBJ9  =  b2rflb.o b2rugm.o b2rups.o b2ruzd.o
OBJ10 =  b2wups.o b2wuzd.o
OBJ11 =  b2ai.o

b2ai.exe: $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9) $(OBJ10) $(OBJ11)
	ifort *.o -o Release/b2ai.exe
#	ifort *.o -o b2ai.exe  -lU77
b2mod_types.o: ../src/modules/b2mod_types.F
	ifort -c ../src/modules/b2mod_types.F
b2mod_constants.o: ../src/modules/b2mod_constants.F
	ifort -c ../src/modules/b2mod_constants.F
b2mod_b2cmfs.o: ../src/modules/b2mod_b2cmfs.F
	ifort -c ../src/modules/b2mod_b2cmfs.F
b2mod_b2cmgs.o: ../src/modules/b2mod_b2cmgs.F
	ifort -c ../src/modules/b2mod_b2cmgs.F
b2mod_b2cmpa.o: ../src/modules/b2mod_b2cmpa.F
	ifort -c ../src/modules/b2mod_b2cmpa.F
b2mod_version.o: ../src/modules/b2mod_version.F
	ifort -c ../src/modules/b2mod_version.F
b2mod_geo_corner.o: ../src/modules/b2mod_geo_corner.F
	ifort -c ../src/modules/b2mod_geo_corner.F
b2mod_geo.o: ../src/modules/b2mod_geo.F
	ifort -c ../src/modules/b2mod_geo.F
b2mod_indirect.o: ../src/modules/b2mod_indirect.F
	ifort -c ../src/modules/b2mod_indirect.F
b2mod_plasma.o: ../src/modules/b2mod_plasma.F
	ifort -c ../src/modules/b2mod_plasma.F
b2mod_rates.o: ../src/modules/b2mod_rates.F
	ifort -c ../src/modules/b2mod_rates.F
b2mod_time.o: ../src/modules/b2mod_time.F
	ifort -c ../src/modules/b2mod_time.F

cfopen.o: ../src/utility/cfopen.F
	ifort -c ../src/utility/cfopen.F
cfruch.o: ../src/utility/cfruch.F
	ifort -c ../src/utility/cfruch.F
cfruin.o: ../src/utility/cfruin.F
	ifort -c ../src/utility/cfruin.F
cfrure.o: ../src/utility/cfrure.F
	ifort -c ../src/utility/cfrure.F
cfvers.o: ../src/utility/cfvers.F
	ifort -c ../src/utility/cfvers.F
cfwuch.o: ../src/utility/cfwuch.F
	ifort -c ../src/utility/cfwuch.F
cfwuin.o: ../src/utility/cfwuin.F
	ifort -c ../src/utility/cfwuin.F
cfwure.o: ../src/utility/cfwure.F
	ifort -c ../src/utility/cfwure.F
chcase.o: ../src/utility/chcase.F
	ifort -c ../src/utility/chcase.F
daytim.o: ../src/utility/daytim.F
	ifort -c ../src/utility/daytim.F
ifill.o: ../src/utility/ifill.F
	ifort -c ../src/utility/ifill.F
intfaceh.o: ../src/utility/intfaceh.F
	ifort -c ../src/utility/intfaceh.F
intfacev.o: ../src/utility/intfacev.F
	ifort -c ../src/utility/intfacev.F
ipgeti.o: ../src/utility/ipgeti.F
	ifort -c ../src/utility/ipgeti.F
ipmain.o: ../src/utility/ipmain.F
	ifort -c ../src/utility/ipmain.F
jobnam.o: ../src/utility/jobnam.F
	ifort -c ../src/utility/jobnam.F
lwimai.o: ../src/utility/lwimai.F
	ifort -c ../src/utility/lwimai.F -I../src/common/
lwmain.o: ../src/utility/lwmain.F
	ifort -c ../src/utility/lwmain.F -I../src/common/
nagsubst.o: ../src/utility/nagsubst.F
	ifort -c -DNO_NAG ../src/utility/nagsubst.F
open_file.o: ../src/utility/open_file.F
	ifort -c ../src/utility/open_file.F
prgend.o: ../src/utility/prgend.F
	ifort -c ../src/utility/prgend.F
prgini.o: ../src/utility/prgini.F
	ifort -c ../src/utility/prgini.F -I../src/common/
prvrt.o: ../src/utility/prvrt.F
	ifort -c ../src/utility/prvrt.F
prvrti.o: ../src/utility/prvrti.F
	ifort -c ../src/utility/prvrti.F
samax.o: ../src/utility/samax.F
	ifort -c ../src/utility/samax.F
sfill.o: ../src/utility/sfill.F
	ifort -c ../src/utility/sfill.F
smax.o: ../src/utility/smax.F
	ifort -c ../src/utility/smax.F
smin.o: ../src/utility/smin.F
	ifort -c ../src/utility/smin.F
ssum.o: ../src/utility/ssum.F
	ifort -c ../src/utility/ssum.F
streql.o: ../src/utility/streql.F
	ifort -c ../src/utility/streql.F
strmas.o: ../src/utility/strmas.F
	ifort -c ../src/utility/strmas.F
subsys.o: ../src/utility/subsys.F
	ifort -c ../src/utility/subsys.F
sysend.o: ../src/utility/sysend.F
	ifort -c ../src/utility/sysend.F
sysini.o: ../src/utility/sysini.F
	ifort -c ../src/utility/sysini.F
xerrab.o: ../src/utility/xerrab.F
	ifort -c ../src/utility/xerrab.F
xerset.o: ../src/utility/xerset.F
	ifort -c ../src/utility/xerset.F
xertst.o: ../src/utility/xertst.F
	ifort -c ../src/utility/xertst.F

usrnam.o: ../src/utility/usrnam.F
	ifort -c ../src/utility/usrnam.F

b2xgbs.o: ../src/b2aux/b2xgbs.F
	ifort -c ../src/b2aux/b2xgbs.F
b2xpne.o: ../src/b2aux/b2xpne.F
	ifort -c ../src/b2aux/b2xpne.F
b2xvff.o: ../src/b2aux/b2xvff.F
	ifort -c ../src/b2aux/b2xvff.F
b2xvfx.o: ../src/b2aux/b2xvfx.F
	ifort -c ../src/b2aux/b2xvfx.F
b2xvfy.o: ../src/b2aux/b2xvfy.F
	ifort -c ../src/b2aux/b2xvfy.F
b2xvgm.o: ../src/b2aux/b2xvgm.F
	ifort -c ../src/b2aux/b2xvgm.F -I../src/common/
b2xvps.o: ../src/b2aux/b2xvps.F
	ifort -c ../src/b2aux/b2xvps.F
b2xvsg.o: ../src/b2aux/b2xvsg.F
	ifort -c ../src/b2aux/b2xvsg.F
b2xxgs.o: ../src/b2aux/b2xxgs.F
	ifort -c ../src/b2aux/b2xxgs.F
b2xxid.o: ../src/b2aux/b2xxid.F
	ifort -c ../src/b2aux/b2xxid.F

b2rflb.o: ../src/input/b2rflb.F
	ifort -c ../src/input/b2rflb.F
b2rugm.o: ../src/input/b2rugm.F
	ifort -c ../src/input/b2rugm.F -I../src/common/
b2rups.o: ../src/input/b2rups.F
	ifort -c ../src/input/b2rups.F -I../src/common/
b2ruzd.o: ../src/input/b2ruzd.F
	ifort -c ../src/input/b2ruzd.F -I../src/common/

b2wups.o: ../src/output/b2wups.F
	ifort -c ../src/output/b2wups.F
b2wuzd.o: ../src/output/b2wuzd.F
	ifort -c ../src/output/b2wuzd.F -I../src/common/

b2aix0.o: ../src/preprocessing/b2aix0.F
	ifort -c ../src/preprocessing/b2aix0.F
b2aidr.o: ../src/preprocessing/b2aidr.F
	ifort -c ../src/preprocessing/b2aidr.F -I../src/common/
b2ai.o: ../src/preprocessing/b2ai.F
	ifort -c ../src/preprocessing/b2ai.F

clean:
	rm -f *.o *.mod





