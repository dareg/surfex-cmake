#!/bin/sh

set -e

MODI_PATH="modi/"
if [ $# -eq 1 ]; then
	MODI_PATH=$1
fi

[ ! -d "$MODI_PATH" ] && mkdir "$MODI_PATH"
python3 gen_modi.py -o "$MODI_PATH" $(find SURFEX/ GELATO/ ASSIM/ OFFLIN/ TRIP/ TOPD/ -name "*.F90")
[ -f "$MODI_PATH/CMakeLists.txt" ] && rm "$MODI_PATH/CMakeLists.txt"
echo "set(modi" >  "$MODI_PATH/CMakeLists.txt"
find "$MODI_PATH"/ | sort >> "$MODI_PATH/CMakeLists.txt"
echo ")" >> "$MODI_PATH/CMakeLists.txt"
