#!/bin/sh

~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
	-in:file:fasta hybrid.fasta \
	-parser::protocol hybrid.xml \
	-parser::script_vars tmpl=$1'_thread.pdb' \
	-constrain_relax_to_start_coords \
	-relax::cartesian \
	-relax:default_repeats 1 \
	-default_max_cycles 200 \
	-score:weights beta_cart \
	-beta \
	-nstruct 1 \
	-out::prefix $1'_' \
	-out::suffix _$2
