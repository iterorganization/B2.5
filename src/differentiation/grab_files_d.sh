MODLIST=`ls b2mod*_d.F90 b2us_*_d.F90 | sed -e 's/_d.F90//g'`
cat $SOLPSTOP/modules/B2.5/src/differentiation/files_to_exclude.txt > tmp
ls b2mod*_d.F90 b2us_*_d.F90 | sed -e 's/_d.F90//g' >> tmp
for d in $MODLIST; do
mv $d"_d.F90" $d"_diff.F90"
done;
ls *_d.F90 | sed -e 's/_d.F90//g' >> tmp
#check that _d routines have their counterpart _nodiff, otherwise they need to be recovered (only for TAPENADE 3.14)
#DLIST=`ls *_d.F90 | sed -e 's/_d.F90//g'`
#for d in $DLIST; do
#f=$d"_nodiff"
#if (grep -q -w -i $f $d"_d.F90") then
#echo $d >> tmp
#fi
#done;
SRCPATH="modules b2aux convert documentation driver equations input output postprocessing preprocessing solvers sources transport utility b2plot user ids"
# grab the files that have been excluded from differentiation but not the MAIN programs
rm -r temp
mkdir temp
cd temp
mv ../tmp .
for d in $SRCPATH; do
find $SOLPSTOP/modules/B2.5/src/$d -name \*.F -exec basename \{} .F \; | ( while read filename
do
if ! (grep -q -w -i $filename tmp) then
file=`find $SOLPSTOP/modules/B2.5/src/$d -name \$filename.F`
cp $file .
echo "Copied filed "$filename".F to differentiation"
fi
done)

find $SOLPSTOP/modules/B2.5/src/$d -name \*.F90 -exec basename \{} .F90 \; | ( while read filename
do
if ! (grep -q -w -i $filename tmp) then
file=`find $SOLPSTOP/modules/B2.5/src/$d -name \$filename.F90`
cp $file .
echo "Copied filed "$filename".F90 to differentiation"
fi
done)
done;

rm tmp

cp $SOLPSTOP/modules/B2.5/src/utility/sfill.F .
cp $SOLPSTOP/modules/B2.5/src/utility/intface.F .
cp $SOLPSTOP/modules/B2.5/src/utility/myblas.F .
cp $SOLPSTOP/modules/B2.5/src/utility/calc_err.F .
cp $SOLPSTOP/modules/B2.5/src/utility/intp_2dtable.F .
cp $SOLPSTOP/modules/B2.5/src/utility/cpeir_bilinear_int.F .
cp $SOLPSTOP/modules/B2.5/src/utility/mstep.F .
cp $SOLPSTOP/modules/B2.5/src/utility/expu.F .
cp $SOLPSTOP/modules/B2.5/src/utility/ma28copy.F .
cp $SOLPSTOP/modules/B2.5/src/transport/b2tlnl.F .
cp $SOLPSTOP/modules/B2.5/src/b2aux/b2xpne.F .
cp $SOLPSTOP/modules/B2.5/src/b2aux/b2xpnm.F .
cp $SOLPSTOP/modules/B2.5/src/b2aux/b2xpnr.F .
cp $SOLPSTOP/modules/B2.5/src/equations/b2nxst.F .
cp $SOLPSTOP/modules/B2.5/src/equations/b2nxcm.F .
cp $SOLPSTOP/modules/B2.5/src/solvers/b2uppo.F .
cp $SOLPSTOP/modules/B2.5/src/sources/b2sqcx.F .
cp $SOLPSTOP/modules/B2.5/src/differentiation/myerf_d.F .
cp $SOLPSTOP/modules/B2.5/src/differentiation/dim_d.F90 .
cp $SOLPSTOP/modules/B2.5/src/differentiation/b2mwqt_diff.F .
cp $SOLPSTOP/modules/B2.5/src/differentiation/b2mxar_diff.F .
cp $SOLPSTOP/modules/B2.5/src/differentiation/b2mxac_diff.F .
cp $SOLPSTOP/modules/B2.5/src/differentiation/b2uxus_d.F .
#cp $SOLPSTOP/modules/B2.5/src/differentiation/invert_a_d.F .
#cp $SOLPSTOP/modules/B2.5/src/differentiation/slv5pt_d.F .
#cp $SOLPSTOP/modules/B2.5/src/differentiation/b2rups_diff.F .


# and now modify the 'use modules' which have been differentiated
files=`ls *.F*`
for d in $MODLIST; do
f=$d"_diff"
echo "Now modifying the use of "$d" into "$f
sed -i -e "s/\<use $d\>/use $f/g" $files
done;

echo "Files that have been excluded from differentiation have been copied to differentiation directory for compiling"

#sed -i -e "s/call get_jsep(/call get_jsep_nodiff(/g" ./*
#sed -i -e "s/call b2usr_loads(/call b2usr_loads_nodiff(/g" ./*
#sed -i -e "s/call species(/call species_nodiff(/g" ./*
#sed -i -e "s/\<call b2xppr\>/call b2xppr_nodiff/g" ./*
#sed -i -e "s/call b2xzef(/call b2xzef_nodiff(/g" ./*
#sed -i -e "s/call b2xpni (/call b2xpni_nodiff (/g" ./*
#sed -i -e "s/call b2xpro(/call b2xpro_nodiff(/g" ./*
#sed -i -e "s/call b2xpfe (/call b2xpfe_nodiff (/g" ./*
#sed -i -e "s/call b2xpfi (/call b2xpfi_nodiff (/g" ./*
#sed -i -e "s/call parsehdr(/call parsehdr_nodiff(/g" ./*
#sed -i -e "s/\<fix_user\>/fix_user_nodiff/g" ./*
#sed -i -e "s/\<call mapx\>/call mapx_nodiff/g" ./*
#sed -i -e "s/\<call mapy\>/call mapy_nodiff/g" ./*

mv ./*.F ../
mv ./*.F90 ../
cd ../
rm -r temp


