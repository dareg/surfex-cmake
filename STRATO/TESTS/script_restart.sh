#!/bin/ksh
 
#set -x 

ycas=$1



for ext in lfi txt fa nc
do
	if [ -f PGD_$ycas.$ext ] 
	then
		cp -f PGD_$ycas.$ext PGD.$ext
	fi
	if [ -f PREP_$ycas.$ext ] 
	then
		cp -f PREP_$ycas.$ext PREP.$ext
	fi
done


