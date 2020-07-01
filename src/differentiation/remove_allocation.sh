#! /bin/csh -f
# this scripts removes the initialization to 0.0_R8 of differentiated variables in b2mod_driver_diff.F90. 
# If this is not done then run continuation with AD b2fdiff* files is not possible, i.e. each time the diff'ed quantities will start from zero.
set num = `egrep -c -ai "cfrure" read_plasma_state_diff.F` #number of times a variable is read from the diff'ed state file
set ii = 1
while ( $ii <= $num )
   # read lines with cfrure, take only the ii-th one and record the variable dimension, the field (pl or dv) and the name
   set record = `grep -i "cfrure" read_plasma_state_diff.F | grep -i -v "time" | awk "NR==$ii{print}" | awk '{sub(/,/, "", $9);print $4, $7, $9}'`
   # count variable dimension using * as separator as in nCv*2 -> fields = 2
   set fields = `echo "$record" | awk 'BEGIN {FS = "*" };{print NF}'`
   if ($fields > 0) then
     # build variable name concatenating the field and the name
     set v1 = `echo "$record" | awk '{print $2}'`
     set v2 = `echo "$record" | awk '{print $3}'`
     set varname = "stated%$v1%$v2"
     # find the line of first occurence of the initialization loop e.g. DO ii1=1,SIZE(varname)
     set l1 = `egrep -ai -m 1 -n "size\($varname" b2mod_driver_diff.F90 | awk '{sub(/:/, "", $1);print $1}'`
     if ($l1 > 0) then
       # assign to $fields the line number after the first occurence where the last END DO is closed
       # DO ii1=..     line l1
       #   DO ii2=..
       #     var = 0.0
       #   END DO
       # END DO        line l1+2*fields+1
       @ l2 = $fields * 2 + $l1
       sed -i "$l1,$l2"d b2mod_driver_diff.F90
     endif
   endif
  @ ii ++
end

# remove a wrong character allocation in b2us_plasm_diff.F90
set l1 = `egrep -ai -m 1 -n "ALLOCATE\(state_extd%text" b2us_plasma_diff.F90 | awk '{sub(/:/, "", $1);print $1}'`
@ l1 = $l1 + 1
@ l2 = $l1 + 4
sed -i "$l1,$l2"d b2us_plasma_diff.F90
sed -i -e "/ ALLOCATE(state_extd%text/a\      state_extd%text = ''" b2us_plasma_diff.F90
