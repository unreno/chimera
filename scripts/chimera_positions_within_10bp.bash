#!/usr/bin/env bash

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0` file1 file2"
	echo
	echo "Expecting file content in the format of ..."
	echo
	echo "..."
	echo "chr8|99063786"
	echo "chr9|105365407"
	echo "chrX|74554211"
	echo "chrY|14564844"
	echo "... OR ..."
	echo "chrX:154433612-155433612|68228"
	echo "chrX:154433612-155433612|790825"
	echo "chrX:93085499-94085499|110644"
	echo "chrX:93085499-94085499|112146"
	echo "..."
	echo
	echo "Example:"
	echo "  `basename $0` file1 file2"
	echo
	exit 1
}

#	Basically, this is TRUE AND DO ...
[ $# -ne 2 ] && usage
if [ ! -f $1 ] ; then 
	echo "$1 is not a file" 
	usage
fi
if [ ! -f $2 ] ; then
	echo "$2 is not a file" 
	usage
fi

#	$1 = queries	
#	$2 hg19 reference

for line in `cat $1` ; do

#	I think that the format changed.
#	chr is first this. pos is last thing.
#	However, I think that it CAN be chr:position:direction

#	echo "line:" $line
	chr=${line%%|*}	#	remove everything after first pipe (including pipe)
#	echo "chr: " $chr
	line=${line#*|}	#	remove everything before first pipe (including pipe)
#	echo "line:" $line
	pos=${line%%|*}	#	remove everything after first pipe (including pipe)
#	echo "pos: " $pos

	#	Expecting file content format like so ...
	#	chrY|6616930|

	awk -F\| -v chr="$chr" -v pos="$pos" '
		( ( $1 == chr ) && ( (pos-10) < $NF ) && ( (pos+10) > $NF ) ){
			print
		}' $2

done
