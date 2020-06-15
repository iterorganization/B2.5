#! /bin/csh -f

set num = `egrep -c -ai "ii.=.,SIZE\(stated.*, [^1]\)" b2mod_driver_diff.F90`

set ii = 1
while ( $ii <= $num )
	set var1 = `egrep -ai -m 1 "ii.=[^0],SIZE\(stated.*, [^1]\)" b2mod_driver_diff.F90`

        set var2 = `echo $var1 | awk '{sub("\\=1,SIZE","=0,SIZE");sub("\\)",")-1");print}' `

	sed -i -e "s/$var1/$var2/g" b2mod_driver_diff.F90
	@ ii ++
end
