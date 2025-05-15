#!/bin/sh
set -x
FILEPATH=$1
FILENAME=$(basename "$FILEPATH")

OUTFILE="modi_${FILENAME%.*}.intfb.h"
EXTENSION="${filename##*.}"
if [ ! $EXTENSION = ".F90" ]; then
	OUTFILE="modi_$FILENAME"
fi
FINAL_FILENAME="modi_$FILENAME"

echo "Generate $OUTFILE from $FILENAME"

./make_intfbl_f90.pl "$FILEPATH" 
if [ ! -f "$LOC_INTFBDIR/$OUTFILE" ]; then
	touch "$LOC_INTFBDIR/$OUTFILE"
fi

#echo "TO MOVE $LOC_INTFBDIR/modi_$FILENAME.intfb.h $LOC_INTFBDIR/modi_$FILENAME.F90"
mv "$LOC_INTFBDIR/$OUTFILE" "$LOC_INTFBDIR/$FINAL_FILENAME"
