for file in *.F90
do
awk  '{if (/^USE LFI_PRECISION$/) {print "USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM"}
	else
         {print $0}}' $file > RES/$file
done
