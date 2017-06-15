#	samtools view -h $base.sam | gawk -v base=$aligned \

#	This script requires that the input bam file be sorted in such a fashion
#	that the reads are grouped together in pairs. (sort by name)
#	As the output fasta files will be used with bowtie2,
#	they too need to be in sync by actual pair.

function reverse_complement(s){
	x=""
	for(i=length(s);i!=0;i--)
		x=x comp[substr(s,i,1)];
	return x;
}
function print_to_fasta(a){
	lane=(and(a[2],64))?"1":"2";
	print ">"a[1]"/"lane >> base"."pre_or_post"_"lane".fasta"
	print a[10]          >> base"."pre_or_post"_"lane".fasta"
}
function trim(r){
	#	Ensure at least 2-digit soft clip and ensure matches near the beginning of the reference.
	#	"near the beginning" means starts at position <= 5
	if( ( r[6] ~ /^[0-9]{2,}S[0-9IDM]*$/ ) && ( r[4] <= 5 ) ){
		split(r[6],a,"S");
		clip=a[1]-r[4]+1;
		r[10]=substr(r[10],1,clip);
		pre_or_post="pre";
	}

	#	Ensure at least 2-digit soft clip and ensure matches near the end of the reference.
	#	"near the end" means starts at position >= 5 more than the reference minus the length of the read
	if( ( r[6] ~ /^[0-9IDM]*[0-9]{2,}S$/ ) && ( r[4] >= ( ref[r[3]] - length(r[10]) + 5 ) ) ){
		clip=ref[r[3]]-r[4]+2;
		r[10]=substr(r[10],clip);
		pre_or_post="post";
	}
#			return r;	#	cannot return arrays. Arrays passed as reference so mods made are actual.
	return pre_or_post;
}
BEGIN {
	comp["A"]="T";
	comp["T"]="A";
	comp["C"]="G";
	comp["G"]="C";
	split("",b);split("",l);
}
#	Store all of the reference lengths
#	Ex. @SQ	SN:chr1	LN:249250621
#	... ref[chr1] = 249250621
( /^@SQ/ ){ ref[substr($2,4)] = substr($3,4); next; }

#	Simply for progress
( ( !/^@/ ) && ( ( NR % 10000000 ) == 0 ) ){ print "Read "NR" records" }

#	Non-sequence reference lines, with a mapped reference, matching previous sequence name
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] == $1 ) ){
	for(i=0;i<=NF;i++)l[i]=$i;

	#	1 and only 1 read aligned. Many ways to check this.
	if( ( and(l[2],4) && !and(l[2],8) ) || ( !and(l[2],4) && and(l[2],8) ) ){
		#	If this read unmapped and mate not unmapped (mate mapped) ...
		if( and(l[2],4) && !and(l[2],8) ){
			before_trim=b[10]
			pre_or_post=trim(b);
			if( b[10] != before_trim ){
				if( and(b[2],16) ) l[10]=reverse_complement(l[10])
				print_to_fasta( b )
				print_to_fasta( l )
			}
		#	If mate read unmapped and this not unmapped (this mapped) ...
		}else if( !and(l[2],4) && and(l[2],8) ){
			before_trim=l[10]
			pre_or_post=trim(l);
			if( l[10] != before_trim ){
				if( and(l[2],16) ) b[10]=reverse_complement(b[10])
				print_to_fasta( b )
				print_to_fasta( l )
			}
		}
	}

	delete b; delete l;
	next; #	do not buffer this line
}

#	Non-sequence reference lines, with a mapped reference, not matching previous sequence name
#	Buffer first occurence of sequence name.
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] != $1 ) ){
	for(i=0;i<=NF;i++)b[i]=$i;
}
#'
#	-> .1.fasta
#	-> .2.fasta
