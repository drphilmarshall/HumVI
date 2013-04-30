"""
Some things that all modules might need.
"""

##========================================================

## Precision cutoff
cutoff = 1e-8

##========================================================

## Current line number
import inspect
def lineno():
	return inspect.currentframe().f_back.f_lineno
