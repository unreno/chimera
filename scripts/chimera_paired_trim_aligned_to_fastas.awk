#	samtools view -h $base.sam | gawk -v base=$aligned \

#	This script requires that the input bam file be sorted in such a fashion
#	that the reads are grouped together in pairs. (sort by name)
#	As the output fasta files will be used with bowtie2,
#	they too need to be in sync by actual pair.
function lg(s){ if( logging ) print s >> log_file }
function reverse_complement(s){
	x=""
	for(i=length(s);i!=0;i--)
		x=x comp[substr(s,i,1)];
	return x;
}
function print_to_fasta(a,dir){
	lane=(and(a[2],64))?"1":"2";
	print ">"a[1]"/"lane >> base"."dir"."pre_or_post"_"lane".fasta"
	print a[10]          >> base"."dir"."pre_or_post"_"lane".fasta"
}
function trim(r){
	#	Ensure at least 2-digit soft clip and ensure matches near the beginning of the reference.
	#	"near the beginning" means starts at position <= 5
	#	20S50M
	if( ( r[6] ~ /^[0-9]{2,}S[0-9IDMX=]*$/ ) && ( r[4] <= 5 ) ){
		lg( "Soft clipped just at beginning and near start of ref")
		lg( r[6] )
		lg( r[4] )
		split(r[6],a,"S");
#	base clip on just CIGAR's soft clipped
		clip=a[1];
#	or on where the reference SHOULD start (was this way)
#		clip=a[1]-r[4]+1;
		lg( "Keeping first "clip" bases" )
		r[10]=substr(r[10],1,clip);
		lg( length(r[10]) )
		pre_or_post="pre";
	}

	#	Ensure at least 2-digit soft clip and ensure matches near the end of the reference.
	#	"near the end" means starts at position >= 5 more than the reference minus the length of the read
	#	50M20S
# why just "if" and not "else if"? Testing "else if"
	else if( ( r[6] ~ /^[0-9IDMX=]*[0-9]{2,}S$/ ) && ( r[4] >= ( ref[r[3]] - length(r[10]) + 5 ) ) ){
		lg( "Soft clipped just at end and near end of ref" )
		lg( r[6] )
		lg( r[4] )
		lg( ( ref[r[3]] - length(r[10]) ) )

#	clip off just the non-S CIGAR? 
		split(r[6],a,/[IDMX=]/);
		clip=sprintf("%d",a[length(a)]);
		lg( "Keeping last "clip" bases" )
		r[10]=substr(r[10],length(r[10])-clip+1);
#	or clip off from where end of reference should be (if didn't match to the end) (was this way)
#		clip=ref[r[3]]-r[4]+2;
#		r[10]=substr(r[10],clip);
		lg( length(r[10]) )
		pre_or_post="post";
	}

	else {
		lg( "Not soft clipped on just one end or not matched at end of reference." )
		pre_or_post="neither";
	}

#	substr(s,a,b) -> from string s, from position a, return b chars. No b = all to end.

#			return r;	#	cannot return arrays. Arrays passed as reference so mods made are actual.
	return pre_or_post;
}
BEGIN {
	FS="\t";	#	20171108 - wasn't here, but likely not needed
	comp["A"]="T";
	comp["T"]="A";
	comp["C"]="G";
	comp["G"]="C";
	comp["N"]="N";
	split("",b);split("",l);
	log_file="chimera_paired_trim_aligned_to_fastas.awk"
#	if( logging ) print "\nNew file\n\n" >> log_file
	lg( "\nNew file\n\n" )
}
#	Store all of the reference lengths
#	Ex. @SQ	SN:chr1	LN:249250621
#	... ref[chr1] = 249250621
( /^@SQ/ ){ ref[substr($2,4)] = substr($3,4); next; }

#	Simply for progress
( ( !/^@/ ) && ( ( NR % 10000000 ) == 0 ) ){ print "Read "NR" records" }

#	Non-sequence reference lines, with a mapped reference, matching previous sequence name
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] == $1 ) ){
	lg( "Matched buffer : "$0 )
	for(i=0;i<=NF;i++)l[i]=$i;

	#	1 and only 1 read aligned. Many ways to check this.
	if( ( and(l[2],4) && !and(l[2],8) ) || ( !and(l[2],4) && and(l[2],8) ) ){
		lg( "Appears 1 read in pair aligned" )
		#	If this read unmapped and mate not unmapped (mate mapped) ...
		if( and(l[2],4) && !and(l[2],8) ){
			before_trim=b[10]
			lg( "Trimming" )
			lg( before_trim )
			pre_or_post=trim(b);
			lg( pre_or_post )
			lg( "Trimmed" )
			lg( b[10] )
			if( b[10] != before_trim ){
#				if( and(b[2],16) ) l[10]=reverse_complement(l[10])
#	Reversing mate to match makes bowtie infer discordant pair alignment.
#	Expects opposite orientation so reverse read back instead.
#	Reverse this read rather than its mate to match.
				dir="F"
				if( and(b[2],16) ){
					b[10]=reverse_complement(b[10])
					dir="R"
				}
				print_to_fasta( b, dir )
				print_to_fasta( l, dir )
			} else {
				lg( "No change in trimming so not written to fasta files" )
			}
		#	If mate read unmapped and this not unmapped (this mapped) ...
		}else if( !and(l[2],4) && and(l[2],8) ){
			before_trim=l[10]
			lg( "Trimming" )
			lg( before_trim )
			pre_or_post=trim(l);
			lg( pre_or_post )
			lg( "Trimmed" )
			lg( l[10] )
			if( l[10] != before_trim ){
#				if( and(l[2],16) ) b[10]=reverse_complement(b[10])
#	Reversing mate to match makes bowtie infer discordant pair alignment.
#	Expects opposite orientation so reverse read back instead.
#	Reverse this read rather than its mate to match.
				dir="F"
				if( and(l[2],16) ){
					l[10]=reverse_complement(l[10])
					dir="R"
				}
				print_to_fasta( b, dir )
				print_to_fasta( l, dir )
			} else {
				lg( "No change in trimming so not written to fasta files" )
			}
		} else {
			lg( "Shouldn't be here" )
		}
	} else {
		lg( "Either neither or both reads aligned" )
	}

	delete b; delete l;
	next; #	do not buffer this line
}

#	Non-sequence reference lines, with a mapped reference, not matching previous sequence name
#	Buffer first occurence of sequence name.
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] != $1 ) ){
	lg( "Buffering : "$0 )
	for(i=0;i<=NF;i++)b[i]=$i;
}
#'
#	-> .1.fasta
#	-> .2.fasta
