#
#	When called like ...
#	awk -f chimera_positions_within_10bp.bash $base.*.bowtie2.$human.$mapq.insertion_points ....
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

BEGIN {
	FS="|"
}
( NR == FNR ){
	positions[$1][$2]++
}
( NR != FNR && $1 in positions ){
	for( pos in positions[$1] )
		if( ( (pos-10) < $2 ) && ( (pos+10) > $2 ) )
			for( i=0; i<positions[$1][pos]; i++ )
				print
}
