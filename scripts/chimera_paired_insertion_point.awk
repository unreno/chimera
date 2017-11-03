
#	Expecting direction (F or R), pre_or_post (pre or post)

BEGIN{
	#	input is a sam file so specifying tab may be unnecessary, but clarity is nice.
	FS="\t";
#	log_file="chimera_paired_insertion_point.awk";
}

#  Non-reference lines, with a mapped reference, not matching previous sequence name
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] != $1 ) ){
	#  Buffer first occurence of sequence name.
	for(i=0;i<=NF;i++)b[i]=$i;

#	print "Buffer new" >> log_file
#	print $0 >> log_file
	next;	#	WILL run next block too, unless next
		#	previously, this resulted in BOTH points included
}  
#  Non-reference lines, with a mapped reference, MATCHING previous sequence name
#
#		Does the pair matching the same reference chromosome matter here?
#		If so, would need to check extended columns for concordant pairing (YT:Z:CP, not YT:Z:DP, YT:Z:UP, ...?)
#		PROPER_PAIR is good, but would result is very small numbers. Trying.
#
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] == $1 ) ){
	for(i=0;i<=NF;i++)l[i]=$i;

#	print "Match" >> log_file
#	print $0 >> log_file

#	For the shortest read ...

#	print "shortest" >> log_file
	if( length(l[10]) < length(b[10]) ){
#		print l[10] >> log_file
		for(i=0;i<=NF;i++)s[i]=l[i];
	} else {
#		print b[10] >> log_file
		for(i=0;i<=NF;i++)s[i]=b[i];
	}

#	print direction >> log_file
#	print pre_or_post >> log_file
#	print "insertion point" >> log_file

#	s[2] is the flag field
#	 and(s[2],16) = REVERSE
#	!and(s[2],16) = NOT REVERSE = FORWARD

	if( \
		( !and(s[2],16) && direction == "F" && pre_or_post == "pre" ) ||
		(  and(s[2],16) && direction == "R" && pre_or_post == "post" ) ){
			print s[3]"|"s[4]+length(s[10])
#			print "position plus length" >> log_file
#			print s[4]+length(s[10]) >> log_file
	} else if( \
		(  and(s[2],16) && direction == "R" && pre_or_post == "pre" ) ||
		( !and(s[2],16) && direction == "F" && pre_or_post == "post" ) ){
			print s[3]"|"s[4]
#			print "position" >> log_file
#			print s[4] >> log_file
	}
	next;
}
