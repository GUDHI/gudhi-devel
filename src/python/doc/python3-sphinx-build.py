#!/usr/bin/env python3

"""
Emulate sphinx-build for python3
"""

from sys import exit, argv
from sphinx import main

if __name__ == '__main__':
    exit(main(argv))
