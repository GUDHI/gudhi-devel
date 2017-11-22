#!/bin/bash

rm README_FOR_UTILITIES.txt
locate utilities/README | grep `svn info | grep '^URL:' | egrep -o '(tags|branches)/[^/]+|trunk' | egrep -o '[^/]+$'`  | xargs cat -- >> README_FOR_UTILITIES.txt