# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 11:08:58 2023

@author: jackd
"""

import os, sys, subprocess
from subprocess import Popen, PIPE, STDOUT

p = Popen(['test.exe', '-i', '47', '8', '5'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
