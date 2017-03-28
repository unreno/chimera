#!/usr/bin/env bash


#
#	Can't pass options when using env
#
#	with the -x the commands are sent to STDERR before execution
#
#	bowtie output goes to stderr for some reason
#	probably because the SAM file usually goes to stdout
#	so, wrap everything in curly braces and direct both
#	to files.
#
#	Explicit redirection within the block will override this
#


#	This really needs to be run in the data directory.
#
#	Eventually, may want to pass number of cpus or threads so
#	execs can use the same number.

script=`basename $0`
basedir=`dirname $0`
human='hg19'
viral='herv_k113'

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "Searches for reads that are part viral, part human."
	echo "Completely ignores any laning."
	echo
	echo "$script [--human STRING] [--viral STRING] <fastq file(s)>"
	echo
	echo "Defaults:"
	echo "  human ..... : $human"
	echo "  viral ..... : $viral"
	echo
	echo "Note: all output files will be based on the working directory's name"
	echo
	echo "Note: fastq files can be gzipped (.gz) or bzipped (.bz2),"
#	Not true, I now join them. This doesn't work in other scripts where they are laned.
#	I would have to add a special parameter to determine which is which.
#	echo " but MUST be separated by commas without white space."
#	echo "They will be passed directly to bowtie2 as such."
#	echo "For example:"
#	echo "  a.fastq.gz,b.fastq.gz"
	echo
	exit
}

while [ $# -ne 0 ] ; do
	case $1 in
		-h|--h*)
			shift; human=$1; shift ;;
		-v|--v*)
			shift; viral=$1; shift ;;
		-t|--t*)
			shift; threads=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done


#       Basically, this is TRUE AND DO ...
#[ $# -eq 0 ] && usage


base=`basename $PWD`

#	Print a lot more stuff
set -x

{
	echo "Starting at ..."
	date

	base="$base.bowtie2.$viral.very_sensitive_local"

	aligned="$base.paired.aligned"

	echo "Seeking insertion points and overlaps"

	for q in 20 10 00 ; do
		echo $q
		mapq="Q${q}"

		samtools view -q $q -F 20 $aligned.pre.bowtie2.$human.bam \
			| awk '( $6 != "*" && $6 != "101M" ){print $3"|"$4+length($10)}' \
			| sort > $aligned.pre.bowtie2.$human.$mapq.insertion_points
		samtools view -q $q -F 20 $aligned.post.bowtie2.$human.bam \
			| awk '( $6 != "*" && $6 != "101M" ){print $3"|"$4}' \
			| sort > $aligned.post.bowtie2.$human.$mapq.insertion_points
		awk -f $basedir/chimera_positions_within.awk $aligned.*.bowtie2.$human.$mapq.insertion_points \
			| sort | uniq -c > $aligned.both.bowtie2.$human.$mapq.insertion_points.overlappers

		samtools view -q $q -F 4 -f 16 $aligned.pre.bowtie2.$human.bam \
			| awk '( $6 != "*" && $6 != "101M" ){print $3"|"$4}' \
			| sort > $aligned.pre.bowtie2.$human.$mapq.rc_insertion_points
		samtools view -q $q -F 4 -f 16 $aligned.post.bowtie2.$human.bam \
			| awk '( $6 != "*" && $6 != "101M" ){print $3"|"$4+length($10)}' \
			| sort > $aligned.post.bowtie2.$human.$mapq.rc_insertion_points
		awk -f $basedir/chimera_positions_within.awk $aligned.*.bowtie2.$human.$mapq.rc_insertion_points \
			| sort | uniq -c > $aligned.both.bowtie2.$human.$mapq.rc_insertion_points.rc_overlappers

	done

	echo
	echo "Finished at ..."
	date

} 1>>$base.$script.out 2>&1 &


#		REFERENCES
#
#	Flags:
#		1    0x1   PAIRED        .. paired-end (or multiple-segment) sequencing technology
#		2    0x2   PROPER_PAIR   .. each segment properly aligned according to the aligner
#		4    0x4   UNMAP         .. segment unmapped
#		8    0x8   MUNMAP        .. next segment in the template unmapped
#		16   0x10  REVERSE       .. SEQ is reverse complemented
#		32   0x20  MREVERSE      .. SEQ of the next segment in the template is reversed
#		64   0x40  READ1         .. the first segment in the template
#		128  0x80  READ2         .. the last segment in the template
#		256  0x100 SECONDARY     .. secondary alignment
#		512  0x200 QCFAIL        .. not passing quality controls
#		1024 0x400 DUP           .. PCR or optical duplicate
#		2048 0x800 SUPPLEMENTARY .. supplementary alignment
#
#
#	Sam file columns
#	1 QNAME String Query template NAME
#	2 FLAG Int bitwise FLAG
#	3 RNAME String Reference sequence NAME
#	4 POS Int 1-based leftmost mapping POSition
#	5 MAPQ Int MAPping Quality
#	6 CIGAR String CIGAR string
#	7 RNEXT String Ref.  name of the mate/next read
#	8 PNEXT Int Position of the mate/next read
#	9 TLEN Int observed Template LENgth
#	10 SEQ String segment SEQuence
#	11 QUAL String
#
