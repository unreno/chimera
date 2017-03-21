#!/usr/bin/env bash


#
#	When called like ...
#	chimera_positions_within_10bp.bash $base.*.bowtie2.$human.$mapq.insertion_points ....
#		$1 will be $base.post.bowtie2.$human.$mapq.insertion_points
#		$2 will be $base.pre.bowtie2.$human.$mapq.insertion_points
#
#	Both can be of the format ...
#		chr8|99063786
#		chr9|105365407
#		chrX|74554211
#		chrY|14564844
#		... OR ...
#		chrX:154433612-155433612|68228
#		chrX:154433612-155433612|790825
#		chrX:93085499-94085499|110644
#		chrX:93085499-94085499|112146
#		... OR EVEN ...
#		chr8|99063786|F
#		chr9|105365407|F
#		chrX|74554211|R
#		chrY|14564844|R
#
for line in `cat $1` ; do
	#	echo "line:" $line
	chr=${line%%|*}	#	remove everything after first pipe (including pipe)
	#	echo "chr: " $chr    #	chrX or chrX:154433612-155433612
	line=${line#*|}	#	remove everything before first pipe (including pipe)
	#	echo "line:" $line   #	74554211 or 790825
	pos=${line%%|*}	#	remove everything after first pipe (including pipe)
	#	echo "pos: " $pos    #	74554211 or 790825

	#	Print out the lines with the same reference chromosome
	#		and a position within 10bp in either direction.
	awk -F\| -v chr="$chr" -v pos="$pos" '
		( ( $1 == chr ) && ( (pos-10) < $2 ) && ( (pos+10) > $2 ) ){
			print
		}' $2

done

