#!/usr/bin/env bash

#
#	Can't pass options when using env.
#	With the -x the commands are sent to STDERR before execution.
#	Using "set -x" has the same effect.
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

human='hg19'
viral='herv_k113'
threads=2

function usage(){
	echo
	echo
	echo "THIS SCRIPT IS STILL IN DEVELOPMENT. DO NOT USE."
	echo "THE VIRAL DATABASE MUST BE ALONE."
	echo
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "Searches for pairs of reads, where one is viral and the other is human,"
	echo "as well as individual chimeric reads."
	echo "Requires pairs of properly sorted, laned fasta or fastq files."
	echo
	echo "`basename $0` [--human STRING] [--viral STRING] [--threads INTEGER] <lane 1 fastq file(s)> <lane 2 fastq file(s)>"
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
	echo "  bowtie2, very sadly, DOES NOT base pairs on NAMES, it bases them on POSITION in file."
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

	base="$base.bowtie2.$viral.__very_sensitive_local"


#	Very Sensitive Local. I want paired chimeras as well as individual chimeras.
#	Will it work?

#	Sadly, can't force paired alignment to be concordant (ie on the same reference sequence)
#	--no-discordant doesn't seem to do anything at all

#	With --all, the pairing can be across different references, which is not helpful.

#	Looks like I will have to separate all sequences into different databases,
#	then loop over and align to each.

#	What about the soft clipping? Will I have to separate? Answers to come
#	on the next episode of "what is jake up to now"

#	split $viral .each do ...

#	IFS=',' read -r -a viral_dbs <<< $viral


	bowtie2 --very-sensitive-local --threads $threads -x $viral \
		$filetype -1 $1 -2 $2 -S $base.sam

	status=$?
	if [ $status -ne 0 ] ; then
		date
		echo "bowtie failed with $status"
		exit $status
	fi

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

#	newbase="$base.unaligned_mate_aligned"
#	samtools view -b -f 4 -F 8 -o $newbase.bam $base.sam

#	Use awk script something like ...
#	Buffer first read's info.
#	On second read, if first and second read meet your demands, print it.
#	Make sure that you skip the header!
#BEGIN {
#	if( !base ){
#		print "Requires gawk and "
#		print "please call with '-v base=YOUR_BASE'"
#		print "samtools view $flag $1 | gawk -v base=YOUR_BASE -f "
#		exit
#	}
#}
#( and( $2 , 64 ) ){
#	b1=$1
#	b10=$10
#	b11=$11
#}
#( and( $2 , 128 ) ){
#	if ( $1 == b1 ){
#		if ( length(b10) != length($10) ){
#			print $1 >> base".diff_length_reads"
#		}
#		if ( length(b11) != length($11) ){
#			print $1 >> base".diff_length_quality"
#		}
#		print "@"b1"/1" >> base"_R1.fastq"
#		print b10 >> base"_R1.fastq"
#		print "+" >> base"_R1.fastq"
#		print b11 >> base"_R1.fastq"
#
#		print "@"$1"/2" >> base"_R2.fastq"
#		print $10 >> base"_R2.fastq"
#		print "+" >> base"_R2.fastq"
#		print $11 >> base"_R2.fastq"
#		b1=b10=b11=""
#	}
}


	samtools view -h -b -f 4 -F 8 -o $base.unaligned_mate_aligned.bam $base.sam
	samtools view -h -b -F 4 -f 8 -o $base.aligned_mate_unaligned.bam $base.sam
	samtools view -h -b -F 12     -o $base.both_aligned.bam           $base.sam

#Usage:   samtools merge [-nurlf] [-h inh.sam] [-b <bamlist.fofn>] <out.bam> <in1.bam> <in2.bam> [<in3.bam> ... <inN.bam>]
#
#Options: -n       sort by read names
#         -r       attach RG tag (inferred from file names)
#         -u       uncompressed BAM output
#         -f       overwrite the output BAM if exist
#         -1       compress level 1
#         -l INT   compression level, from 0 to 9 [-1]
#         -@ INT   number of BAM compression threads [0]
#         -R STR   merge file in the specified region STR [all]
#         -h FILE  copy the header in FILE to <out.bam> [in1.bam]
#         -c       combine RG tags with colliding IDs rather than amending them
#         -p       combine PG tags with colliding IDs rather than amending them
#         -s VALUE override random seed
#         -b FILE  list of input BAM filenames, one per line [null]

	samtools merge -p -n -h $base.both_aligned.bam $base.align_mix.bam \
		$base.unaligned_mate_aligned.bam \
		$base.aligned_mate_unaligned.bam \
		$base.both_aligned.bam 


#	rm $base.sam
#	base=$newbase

#	#samtools view $base.bam | awk '{print "@"$1;print $10;print "+";print $11;}' > $base.fastq
#	samtools view $base.bam | gawk '{l=(and($2,64))?"1":"2";print ">"$1"/"l;print $10;}' > $base.fasta


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

#	#	Align the reads to the human reference.
#	bowtie2 -x $human --threads $threads -f -U $base.fasta \
#		-S $base.bowtie2.$human.sam
#	status=$?
#	if [ $status -ne 0 ] ; then
#		date
#		echo "bowtie failed with $status"
#		exit $status
#	fi
#
#	#	Convert to bam and remove the sam.
#	samtools view -b -o $base.bowtie2.$human.bam \
#		$base.bowtie2.$human.sam
#	rm $base.bowtie2.$human.sam

	echo
	echo "Finished at ..."
	date

} 1>>$base.`basename $0`.out 2>&1 &
