#! /bin/csh -f

set num = `egrep -c -ai "ii.=.,SIZE\(stated.*, [^1]\)" b2mod_driver_diffv.F90`

set ii = 1
while ( $ii <= $num )
  set var1 = `egrep -ai -m 1 "ii.=[^0],SIZE\(stated.*, [^1]\)" b2mod_driver_diffv.F90`
  set var2 = `echo $var1 | awk '{sub("\\=1,SIZE","=0,SIZE");sub("\\)",")-1");print}' `
  sed -i -e "s/$var1/$var2/g" b2mod_driver_diffv.F90
  @ ii ++
end

set num = `egrep -c -ai "ii.=.,SIZE\(.*" b2mod_driver_diffv.F90`

set ii = 1
while ( $ii <= $num )
   set record = `egrep -ai "ii.=.,SIZE\(.*" b2mod_driver_diffv.F90 | awk "NR==$ii{print}"`
   set field = `echo $record | awk '{print $3}' | awk 'BEGIN {FS = ")" };{print $1}'`
   @ field = $field + 1
   set dum = '{sub(".)","'$field')");print}'
   set var = `echo $record | awk "$dum" `
   sed -i -e "s/$record/$var/g" b2mod_driver_diffv.F90
  @ ii ++
end
