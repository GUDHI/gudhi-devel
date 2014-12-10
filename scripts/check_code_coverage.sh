#!/bin/bash
#usage bash check_google_style.sh

# UNITARY TEST DIRECTORY CHECK
ROOT_DIR=..
UT_DIR_TO_CHECK="$1"
LCOV_DIR=$UT_DIR_TO_CHECK/lcov
COVERAGE_FILE=$LCOV_DIR/coverage.info
LOG_FILE=$LCOV_DIR/coverage.log
MIN_PERCENT=95

rm -rf $LCOV_DIR
mkdir $LCOV_DIR

lcov --capture --directory $UT_DIR_TO_CHECK --no-external --output-file $COVERAGE_FILE 2>&1 | tee -a $LOG_FILE
lcov --summary $COVERAGE_FILE 2>&1 | tee -a $LOG_FILE
genhtml $COVERAGE_FILE --output-directory $UT_DIR_TO_CHECK/lcov
# CLEAN AFTER USE
lcov --directory $UT_DIR_TO_CHECK --zerocounters 2>&1 | tee -a $LOG_FILE

LINE_PERCENTAGE=`grep "lines......:" $LOG_FILE`
PERC_PER_LINE=${LINE_PERCENTAGE:14:3}

if [ "$PERC_PER_LINE" -lt "$MIN_PERCENT" ]
then
  echo "## Lines not enough covered ($PERC_PER_LINE) in $UT_DIR_TO_CHECK" 2>&1 | tee -a $LOG_FILE
  exit 1
else
  FONC_PERCENTAGE=`grep "functions..:" $LOG_FILE`
  PERC_PER_FUNC=${FONC_PERCENTAGE:14:3}

  if [ "$PERC_PER_FUNC" -lt "$MIN_PERCENT" ]
  then
    echo "## Functions not enough covered ($PERC_PER_FUNC) in $UT_DIR_TO_CHECK" 2>&1 | tee -a $LOG_FILE
    exit 1
  fi
fi

echo "Acceptable coverage values (lines:$PERC_PER_LINE% - functions:$PERC_PER_FUNC%) in $UT_DIR_TO_CHECK" 2>&1 | tee -a $LOG_FILE
exit 0
