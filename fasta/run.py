#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A conditional import that picks either the classical `sh` package or the `runps`
replacement, so that this package can work on all platforms including Windows.
"""

import platform
if platform.system() == 'Windows': import runps as sh
else:
    # After sh==1.14.3 the object returned changed #
    import sh
    sh_version = int(sh.__version__.split('.')[0])
    if sh_version > 1: sh = sh.bake(_return_cmd=True)