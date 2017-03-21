#  Non-reference lines, with a mapped reference, not matching previous sequence name
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] != $1 ) ){
	#  Buffer first occurence of sequence name.
	for(i=0;i<=NF;i++)b[i]=$i;
}  
#  Non-reference lines, with a mapped reference, MATCHING previous sequence name
( ( !/^@/ ) && ( $3 != "*" ) && ( b[1] == $1 ) ){
	for(i=0;i<=NF;i++)l[i]=$i;

#	For the shortest read ...

#	IF Forward and pre, print
#			| awk '( $6 != "*" && $6 != "101M" ){print $3"|"$4+length($10)}' \

#	If forward and post, print
#			| awk '( $6 != "*" && $6 != "101M" ){print $3"|"$4}' \

#	If reverse and pre, print
#			| awk '( $6 != "*" && $6 != "101M" ){print $3"|"$4}' \

#	If reverse and post, print
#			| awk '( $6 != "*" && $6 != "101M" ){print $3"|"$4+length($10)}' \

}
