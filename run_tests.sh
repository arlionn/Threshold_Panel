#!/bin/bash
set -e

DIR=$(dirname $(realpath "$0")) 	# locate folder where this sh-script is located in

PROJECT="Threshold_Panel"

SCRIPT1="./tests/run_tests.inp"

cd $DIR
echo "Switched to ${DIR}"


gretlcli -b -e -q ${SCRIPT1}
if [ $? -eq 0 ]
then
  echo "Success: Tests for script ${SCRIPT1} passed."
  exit 0
else
  echo "Failure: Tests for script ${SCRIPT1} not passed."
  exit 1
fi

