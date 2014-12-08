#!/bin/bash
#usage bash check_google_style.sh

# VERSION CHECK
ROOT_DIR=..
FILE_TO_CHECK="$1"
LOG_FILE=$FILE_TO_CHECK.cpplint

if [ -f $LOG_FILE ]; then
    rm -f $LOG_FILE
fi

if [ ! -f $FILE_TO_CHECK ]; then
    echo "File not found! : $FILE_TO_CHECK" | tee -a $LOG_FILE
    exit 1
fi

# CPPLINT FILE
echo "File: $FILE_TO_CHECK" 2>&1 | tee -a $LOG_FILE
python ~/cpplint.py --linelength=120 $FILE_TO_CHECK 2>&1 | tee -a $LOG_FILE

LINE_ERRORS=`grep "Total errors found:" $LOG_FILE`
NB_ERRORS=${LINE_ERRORS:20}

if [ "$NB_ERRORS" -gt 20 ]
then
  echo "## Too many errors ($NB_ERRORS) in $FILE_TO_CHECK" 2>&1 | tee -a $LOG_FILE
  exit 1
fi

echo "Acceptable number of errors ($NB_ERRORS) in $FILE_TO_CHECK" 2>&1 | tee -a $LOG_FILE
exit 0
