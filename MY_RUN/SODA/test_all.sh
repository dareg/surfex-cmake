#!/bin/bash

tests="OI OI_POINT OI_POINTS EKF EKF_PATCH EKF_POINT EKF_POINTS EKF_PATCH_POINT EKF_PATCH_POINTS"
formats="NC LFI ASCII FA"
nprocs="1 2"

ok=""
failed=""
for t in $tests; do
  for f in $formats; do
    for proc in $nprocs; do
      skip=0
      [ "$f" == "FA" -a "$t" == "EKF_POINT" ] && skip=1
      [ "$f" == "FA" -a "$t" == "EKF_POINTS" ] && skip=1
      [ "$f" == "FA" -a "$t" == "EKF_PATCH_POINT" ] && skip=1
      [ "$f" == "FA" -a "$t" == "EKF_PATCH_POINTS" ] && skip=1

      if [ "$skip" -eq "0" ]; then
        echo "Testing $t with format $f NPROC=$proc"
        OK=1
        ./run_soda.sh $t $f $proc > OUTPUT/${t}_${f}_${proc}.log 2>&1 || OK=0
        [ "$OK" -eq "1" ] && ok="$ok $t $f $proc\n" && echo "Test was OK"  
        [ "$OK" -eq "0" ] && failed="$failed $t $f $proc\n" && echo "Test failed"
  
        if [ -f REF/${t}_${f}_1.log -a $proc -eq 1 ]; then
          diff=${GRAPHICAL_DIFF_PROGRAM:-diff}
          $diff OUTPUT/${t}_${f}_1.log REF/${t}_${f}_1.log
        fi
      else
        echo "Skipping $t with format $f NPROC=$proc because of no implementation"
      fi
    done
  done
done

echo "Testing finished..."
if [ "$ok" != "" ]; then
  echo "==========================="
  echo "The follwing tests were OK:"
  printf "$ok"
  echo "==========================="
fi
if [ "$failed" != "" ]; then
  echo "==========================="
  echo "The follwing tests failed:"
  printf "$failed"
  echo "==========================="
else
  echo
  echo "All tests were OK!"
fi
