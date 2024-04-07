#!/bin/sh
'''exec' "\\tamuqfiles\project\Chemical Engineering\Air_Quality_Eng_R_Team\!GITHUB\tamuq-chen-secarelab-receiver\aysha\.CondaPkg\env\python.exe" "$0" "$@" #'''
# -*- coding: utf-8 -*-
import re
import sys

from wheel.cli import main

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(main())
