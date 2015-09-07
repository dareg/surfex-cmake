#!/bin/bash

tests="OI OI_POINT OI_POINTS EKF EKF_PATCH EKF_POINT EKF_POINTS EKF_PATCH_POINT EKF_PATCH_POINTS"
formats="NC LFI ASCII FA"
nprocs="1 2"

ok=""
failed=""
for t in $tests; do
  for f in $formats; do
    for proc in $nprocs; do
      echo "Testing $t with format $f NPROC=$proc"
      OK=1
      ./run_soda.sh $t $f $proc > OUTPUT/${t}_${f}_${proc}.log 2>&1 || OK=0
      [ "$OK" -eq "1" ] && ok="$ok $t $f\n" && echo "Test was OK"  
      [ "$OK" -eq "0" ] && failed="$failed $t $f\n" && echo "Test failed"

      if [ -f REF/${t}_${f}.log ]; then
        diff=${GRAPHICAL_DIFF_PROGRAM:-diff}
        $diff OUTPUT/${t}_${f}.log REF/${t}_${f}_1.log
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
