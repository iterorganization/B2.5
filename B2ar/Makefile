FC = ifort

OBJ1  =  b2mod_types.o b2mod_constants.o b2mod_version.o b2mod_b2cmpa.o b2mod_b2cmpb.o b2mod_b2cmpt.o b2mod_b2cmrc.o b2mod_indirect.o b2mod_solpstop.o
OBJ2  =  adpak.o cfopen.o cfruch.o cfruin.o cfrure.o cfvers.o cfwuch.o cfwuin.o cfwure.o chcase.o daytim.o ifill.o int2d.o  
OBJ3  =  ipgeti.o ipmain.o ipos.o jobnam.o lwimai.o lwmain.o machsfr.o nagsubst.o open_file.o  
OBJ4  =  prgend.o prgini.o prvrt.o prvrti.o ratadas.o ratadp.o ratprt.o ratstr.o rattai.o ratwei.o 
OBJ5  =  sfill.o smax.o smin.o streql.o strip_spaces.o subsys.o sysend.o sysini.o 
OBJ6  =  xerrab.o xerset.o xertst.o xxdata_11.o
OBJ7  =  usrnam.o
OBJ8  =  b2xvcp.o b2xvsg.o b2xxid.o
OBJ9  =  b2arin.o b2ardr.o
OBJ10 =  b2rucp.o b2ruzd.o
OBJ11 =  b2wfr0.o b2wfr1.o b2wfr2.o b2wurc.o
OBJ12 =  b2ar.o

b2ar.exe: $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9) $(OBJ10) $(OBJ11) $(OBJ12)
	ifort *.o -o Release/b2ar.exe
#	ifort *.o -o b2ar.exe  -lU77

b2mod_types.o: ../src/modules/b2mod_types.F
	ifort -c ../src/modules/b2mod_types.F
b2mod_constants.o: ../src/modules/b2mod_constants.F
	ifort -c ../src/modules/b2mod_constants.F
b2mod_version.o: ../src/modules/b2mod_version.F
	ifort -c ../src/modules/b2mod_version.F
b2mod_b2cmpa.o: ../src/modules/b2mod_b2cmpa.F
	ifort -c ../src/modules/b2mod_b2cmpa.F
b2mod_b2cmpb.o: ../src/modules/b2mod_b2cmpb.F
	ifort -c ../src/modules/b2mod_b2cmpb.F
b2mod_b2cmpt.o: ../src/modules/b2mod_b2cmpt.F
	ifort -c ../src/modules/b2mod_b2cmpt.F
b2mod_b2cmrc.o: ../src/modules/b2mod_b2cmrc.F
	ifort -c ../src/modules/b2mod_b2cmrc.F
b2mod_indirect.o: ../src/modules/b2mod_indirect.F
	ifort -c ../src/modules/b2mod_indirect.F
b2mod_solpstop.o: ../src/modules/b2mod_solpstop.F
	ifort -c ../src/modules/b2mod_solpstop.F

adpak.o: ../src/utility/adpak.F
	ifort -c ../src/utility/adpak.F -I../src/common/
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
int2d.o: ../src/utility/int2d.F
	ifort -c ../src/utility/int2d.F
ipgeti.o: ../src/utility/ipgeti.F
	ifort -c ../src/utility/ipgeti.F
ipmain.o: ../src/utility/ipmain.F
	ifort -c ../src/utility/ipmain.F
ipos.o: ../src/utility/ipos.F
	ifort -c ../src/utility/ipos.F
jobnam.o: ../src/utility/jobnam.F
	ifort -c ../src/utility/jobnam.F
lwimai.o: ../src/utility/lwimai.F
	ifort -c ../src/utility/lwimai.F -I../src/common/
lwmain.o: ../src/utility/lwmain.F
	ifort -c ../src/utility/lwmain.F -I../src/common/
machsfr.o: ../src/utility/machsfr.F
	ifort -c ../src/utility/machsfr.F
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
ratadas.o: ../src/utility/ratadas.F
	ifort -c ../src/utility/ratadas.F
ratadp.o: ../src/utility/ratadp.F
	ifort -c ../src/utility/ratadp.F
ratprt.o: ../src/utility/ratprt.F
	ifort -c ../src/utility/ratprt.F
ratstr.o: ../src/utility/ratstr.F
	ifort -c ../src/utility/ratstr.F -I../src/common/
rattai.o: ../src/utility/rattai.F
	ifort -c ../src/utility/rattai.F
ratwei.o: ../src/utility/ratwei.F
	ifort -c ../src/utility/ratwei.F
sfill.o: ../src/utility/sfill.F
	ifort -c ../src/utility/sfill.F
smax.o: ../src/utility/smax.F
	ifort -c ../src/utility/smax.F
smin.o: ../src/utility/smin.F
	ifort -c ../src/utility/smin.F
streql.o: ../src/utility/streql.F
	ifort -c ../src/utility/streql.F
strip_spaces.o: ../src/utility/strip_spaces.F
	ifort -c ../src/utility/strip_spaces.F
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
xxdata_11.o: ../src/utility/xxdata_11.F
	ifort -c ../src/utility/xxdata_11.F

usrnam.o: ../src/utility/usrnam.F
	ifort -c ../src/utility/usrnam.F

b2xvcp.o: ../src/b2aux/b2xvcp.F
	ifort -c ../src/b2aux/b2xvcp.F -I../src/common/
b2xvsg.o: ../src/b2aux/b2xvsg.F
	ifort -c ../src/b2aux/b2xvsg.F
b2xxid.o: ../src/b2aux/b2xxid.F
	ifort -c ../src/b2aux/b2xxid.F

b2rucp.o: ../src/input/b2rucp.F
	ifort -c ../src/input/b2rucp.F -I../src/common/
b2ruzd.o: ../src/input/b2rusr.F
	ifort -c ../src/input/b2ruzd.F -I../src/common/

b2wfr0.o: ../src/output/b2wfr0.F
	ifort -c ../src/output/b2wfr0.F
b2wfr1.o: ../src/output/b2wfr1.F
	ifort -c ../src/output/b2wfr1.F
b2wfr2.o: ../src/output/b2wfr2.F
	ifort -c ../src/output/b2wfr2.F
b2wurc.o: ../src/output/b2wurc.F
	ifort -c ../src/output/b2wurc.F -I../src/common/

b2arin.o: ../src/preprocessing/b2arin.F
	ifort -c ../src/preprocessing/b2arin.F
b2ardr.o: ../src/preprocessing/b2ardr.F
	ifort -c ../src/preprocessing/b2ardr.F -I../src/common/
b2ar.o: ../src/preprocessing/b2ar.F
	ifort -c ../src/preprocessing/b2ar.F

clean:
	rm -f *.o *.mod





