
#	Expecting direction (F or R), pre_or_post (pre or post)

BEGIN{
	#	input is a sam file so specifying tab may be unnecessary, but clarity is nice.
	FS="\t"
	log_file="chimera_paired_insertion_point.awk"
	if( logging ) print "\nNew file\n\n" >> log_file
}

#  Non-reference lines, with a mapped reference, not matching previous sequence name
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] != $1 ) ){
	#  Buffer first occurence of sequence name.
	for(i=0;i<=NF;i++)b[i]=$i;

	if( logging ) print "Buffer new" >> log_file
	if( logging ) print $0 >> log_file
	next;	#	WILL run next block too, unless next
		#	previously, this resulted in BOTH points included
}  
#  Non-reference lines, with a mapped reference, MATCHING previous sequence name
#
#		Does the pair matching the same reference chromosome matter here?
#		If so, would need to check extended columns for concordant pairing (YT:Z:CP, not YT:Z:DP, YT:Z:UP, ...?)
#		PROPER_PAIR flag is good (same as CP), but would result is very small numbers. Trying.
#
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] == $1 ) ){
	for(i=0;i<=NF;i++)l[i]=$i;

	if( logging ) print "Match" >> log_file
	if( logging ) print $0 >> log_file

#	For the shortest read ...

	if( logging )	print "shortest" >> log_file
	if( length(l[10]) < length(b[10]) ){
		if( logging ) print l[10] >> log_file
		for(i=0;i<=NF;i++)s[i]=l[i];
	} else {
		if( logging ) print b[10] >> log_file
		for(i=0;i<=NF;i++)s[i]=b[i];
	}

	if( logging ) print direction >> log_file
	if( logging ) print pre_or_post >> log_file
	if( logging ) print "insertion point" >> log_file

#	s[2] is the flag field
#	 and(s[2],16) = REVERSE
#	!and(s[2],16) = NOT REVERSE = FORWARD

	#	PRE means that it was BEFORE the reference and trimmed on the RIGHT
	#	POST means that it was AFTER the reference and trimmed on the LEFT

	if( \
		( !and(s[2],16) && direction == "F" && pre_or_post == "pre" ) ||
		(  and(s[2],16) && direction == "R" && pre_or_post == "post" ) ){
			print s[3]"|"s[4]+length(s[10])
			if( logging ) print "position plus length" >> log_file
			if( logging ) print s[4]+length(s[10]) >> log_file
	} else if( \
		(  and(s[2],16) && direction == "R" && pre_or_post == "pre" ) ||
		( !and(s[2],16) && direction == "F" && pre_or_post == "post" ) ){
			print s[3]"|"s[4]
			if( logging ) print "position" >> log_file
			if( logging ) print s[4] >> log_file
	}
	next;
}
