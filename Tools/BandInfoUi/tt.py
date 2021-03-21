import sys, os

import ctypes

lib = ctypes.CDLL("gracedyn.so")
lib.GraceOpen(2048);

     
