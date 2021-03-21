#!/bin/sh

infile=$1
tag=$2

infilesym=`echo $1 | sed 's/_input/_symm/'`
cstfile='state'$2'.cst'
wt_exp=40
wt_cst=2
maxdis=6
if ["$2" = "1"] || ["$2" = "3"]; then
  orig=2013
else
  orig=1973

mkdir state$2

cat $cstfile | sed 's/SCALARWEIGHTEDFUNC 1.0/SCALARWEIGHTEDFUNC '$wt_exp'/' > combo$2.$tag.cst
./make_coord_csts.py --pdb $infile --residues 1-193 --atoms CA,CB,C,H,O,N --weight 0.1 --vrt $orig >> combo$2.$tag.cst
./make_csts.py --pdb $infilesym --seqsep 3 --residues 196-345 --maxdis $maxdis --atoms CA,CB,C,H,O,N --weight $wt_cst >> combo$2.$tag.cst
./make_csts.py --pdb $infilesym --seqsep 3 --residues 349-491 --toall --maxdis $maxdis --atoms CA,CB,C,H,O,N --weight $wt_cst >> combo$2.$tag.cst

~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
    -database ~/Rosetta/main/database \
  	-in:file:s $infile \
  	-extra_res_fa CMP.params \
  	-parser::protocol relax.xml \
  	-parser::script_vars cstfile=combo$2.$tag.cst \
 	  -default_max_cycles 200 \
  	-beta_cart \
  	-overwrite \
  	-out:suffix _$tag \
  	-out:prefix state$2/ \
  	-nstruct 1
