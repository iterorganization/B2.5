string=''

$SOLPSTOP/modules/B2.5/src/differentiation/exclude_files.sh

find $SOLPSTOP/modules/B2.5/builds/standalone.LEUVEN.ifort64 -name \*.f -exec basename \{} .f \;| ( while read filename 
do
if grep -q -w "$filename" excluded.txt
then
string+=""
else
string+=" $SOLPSTOP/modules/B2.5/builds/standalone.LEUVEN.ifort64/"$filename".f"
fi
done

echo $string > testfile)


find $SOLPSTOP/modules/B2.5/builds/standalone.LEUVEN.ifort64 -name \*.f90 -exec basename \{} .f90 \;| ( while read filename 
do
if grep -q -w "$filename" excluded.txt
then
string+=""
else
string+="  $SOLPSTOP/modules/B2.5/builds/standalone.LEUVEN.ifort64/"$filename".f90"
fi
done

echo $string >> testfile)

cat testfile

