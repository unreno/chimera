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
human='hg19'
viral='herv_k113'
threads=2

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "Searches for pairs of reads, where one is viral and the other is human."
	echo "Requires pairs of properly sorted, laned fasta or fastq files."
	echo
	echo "$script [--human STRING] [--viral STRING] [--threads INTEGER] <lane 1 fastq file(s)> <lane 2 fastq file(s)>"
	echo
	echo "Defaults:"
	echo "  human ..... : $human"
	echo "  viral ..... : $viral"
	echo "  threads ... : $threads (for bowtie2)"
	echo
	echo "Note: all output files will be based on the working directory's name"
	echo
	echo "Note: fastq files can be gzipped (.gz) or bzipped (.bz2),"
	echo " but MUST be listed in the same order to preserve synchronization"
	echo " and separated by commas without white space."
	echo "They will be passed directly to bowtie2 as such."
	echo "For example:"
	echo "  a_1.fastq.gz,b_1.fastq.gz a_2.fastq.gz,b_2.fastq.gz"
	echo
	echo "WARNING:"
	echo "  bowtie2 DOES NOT base pairs on NAMES, it bases them on POSITION in file."
	echo "  ie. The first read in each file are considered a pair."
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
[ $# -ne 2 ] && usage



base=`basename $PWD`

#	Print a lot more stuff
set -x

{
	echo "Starting at ..."
	date

	#indexes=/Volumes/cube/working/indexes
	#	leading with the ": " stops execution
	#	just ${BOWTIE2_INDEXES:"/Volumes/cube/working/indexes"}
	#	would try to execute the result.  I just want the OR/EQUALS feature
	#: ${BOWTIE2_INDEXES:="/Volumes/cube/working/indexes"}
	: ${BOWTIE2_INDEXES:="$HOME/BOWTIE2_INDEXES"}

	#	they MUST be exported, to be picked up by bowtie2
	export BOWTIE2_INDEXES

	#	One line "if-then-else" to determine filetype by last character of first file name.
	[ "${1:(-1)}" == "q" -o "${1:(-5)}" == "q.bz2" -o "${1:(-4)}" == "q.gz" ] \
		&& filetype='-q' || filetype='-f'

	base="$base.bowtie2.$viral.__very_sensitive"

	bowtie2 --very-sensitive --threads $threads -x $viral \
		$filetype -1 $1 -2 $2 -S $base.sam
	status=$?
	if [ $status -ne 0 ] ; then
		date
		echo "bowtie failed with $status"
		exit $status
	fi

	#	Convert to bam and remove the sam.
	samtools view -b -o $base.bam $base.sam
	rm $base.sam

	#	BEGIN extraction of JUST UNALIGNED READ and alignment to HUMAN
	umabase="$base.unaligned_mate_aligned"

	samtools view -f 4 -F 8 $base.bam | gawk '
		{	l=(and($2,64))?"1":"2";
			print ">"$1"/"l;
			print $10;}' > $umabase.fasta

	#	Align the reads to the human reference.
	bowtie2 -x $human --threads $threads -f -U $umabase.fasta \
		-S $umabase.bowtie2.$human.sam
	status=$?
	if [ $status -ne 0 ] ; then
		date
		echo "bowtie failed with $status"
		exit $status
	fi

	#	Convert to bam and remove the sam.
	samtools view -b -o $umabase.bowtie2.$human.bam \
		$umabase.bowtie2.$human.sam
	rm $umabase.bowtie2.$human.sam
	#	END extraction of JUST UNALIGNED READ and alignment to HUMAN



	#	BEGIN extraction of ALIGNED READ AND UNALIGNED MATE and alignment to HUMAN

	#	If one read aligns and the mate doesn't, 
	#	the mate will still have a non-"*" value in the reference column $3 so ...

	#	READ1 IS NOT ALWAYS FIRST!!!!!!  ERRRREEREREREREREWR!

	#	bitwise math requires gawk, or perhaps a very new version awk

	habase=$base.half_aligned
	samtools view $base.bam | awk '( $3 != "*" )' | gawk  -v base=$habase '
		function print_to_fasta(a){
			lane=(and(a[2],64))?"1":"2";
			print ">"a[1]"/"lane >> base"_"lane".fasta"
			print a[10]          >> base"_"lane".fasta"
		}
		BEGIN{
			split("",b);split("",l);
		}
		( b[1] == $1 ){
			for(i=0;i<=NF;i++)l[i]=$i;

			#	many ways to check this
			if( ( and(l[2],4) && !and(l[2],8) ) || ( !and(l[2],4) && and(l[2],8) ) ){
				print_to_fasta( b )
				print_to_fasta( l )
			}

			delete b; delete l;
			next; #	do not buffer this line
		}
		( b[1] != $1 ){ 
			for(i=0;i<=NF;i++)b[i]=$i;
		}'  

	#	Align the reads to the human reference.
	bowtie2 -x $human --threads $threads -f -1 ${habase}_1.fasta -2 ${habase}_2.fasta \
		-S $habase.bowtie2.$human.sam
	status=$?
	if [ $status -ne 0 ] ; then
		date
		echo "bowtie or samtools failed with $status"
		exit $status
	fi

	#	Convert to bam and remove the sam.
	samtools view -b -o $habase.bowtie2.$human.bam \
		$habase.bowtie2.$human.sam
	rm $habase.bowtie2.$human.sam
	#	END extraction of ALIGNED READ AND UNALIGNED MATE and alignment to HUMAN

#	rm $base.bam

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
#    f4 = unmapped
#    F4 = NOT unmapped = mapped
#    F8 = mate NOT unmapped = mate mapped
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
