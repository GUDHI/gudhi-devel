#!/usr/bin/env python

# Convert output from Google's cpplint.py to the cppcheck XML format for
# consumption by the Jenkins cppcheck plugin.

# Reads from stdin and writes to stderr (to mimic cppcheck)


import sys
import re

def cpplint_score_to_cppcheck_severity(score):
    # I'm making this up
    if score == 1:
        return 'style'
    elif score == 2:
        return 'style'
    elif score == 3:
        return 'warning'
    elif score == 4:
        return 'warning'
    elif score == 5:
        return 'error'
  

def parse():
    # TODO: do this properly, using the xml module.
    # Write header
    sys.stderr.write('''<?xml version="1.0" encoding="UTF-8"?>\n''')
    # VR : sys.stderr.write('''<results>\n''')
    # Add from VR  + [
    sys.stderr.write('''<results version=2>\n''') 
    sys.stderr.write('''<cppcheck version="1.63"/>\n''')
    sys.stderr.write('''<errors>\n''')
    # -]

    # Do line-by-line conversion
    r = re.compile('([^:]*):([0-9]*):  ([^\[]*)\[([^\]]*)\] \[([0-9]*)\].*')

    for l in sys.stdin.readlines():
        m = r.match(l.strip())
        if not m:
            continue
        g = m.groups()
        if len(g) != 5:
            continue
        fname, lineno, msg, label, score = g  
        severity = cpplint_score_to_cppcheck_severity(int(score))
        # VR : sys.stderr.write('''<error file="%s" line="%s" id="%s" severity="%s" msg="%s"/>\n'''%(fname, lineno, label, severity, msg))
        # Add from VR  + [
        sys.stderr.write('''<error file="%s" line="%s" id="%s" severity="%s" msg="%s"/>\n</error>\n'''%(fname, lineno, label, severity, msg))
        # -]

    # Write footer
    # Add from VR  + [
    sys.stderr.write('''</errors>\n'''%(fname, lineno, label, severity, msg))
    # -]
    sys.stderr.write('''</results>\n''')


if __name__ == '__main__':
    parse()

