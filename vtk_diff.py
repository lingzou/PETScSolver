#!/usr/bin/env python

import sys

##################################
##################################
def is_number(s):
  try:
    float(s)
    return True
  except ValueError:
    return False
##################################
##################################


# Command line check, it takes two files to compare
if len(sys.argv) != 3 :
  print ("Usage: vtk_diff.py <file> <expected_file>")
  sys.exit(1)

# Read in the two files to be compared
try :
  file = open(sys.argv[1], "r")
except IOError:
  print ("Failed to open file: " + sys.argv[1])
  exit(1)
lines = file.readlines()
file.close()

try :
  exp_file = open(sys.argv[2], "r")
except IOError:
  print ("Failed to open file: " + sys.argv[2])
  exit(1)
exp_lines = exp_file.readlines()
exp_file.close()

# First a quick check to see the two files have the same number of lines
if (len(lines) != len(exp_lines)) :
  print ("The two files do not have the same number of lines.")
  sys.exit(1)

diff = False
for i in range(0, len(lines)) :
  if (is_number(lines[i]) and is_number(exp_lines[i])) :
    val     = float(lines[i])
    exp_val = float(exp_lines[i])
    if (abs(exp_val) < 1.0e-8) :
      # FIXME: This criterion seems still overly rigorous
      if (abs(val - exp_val) > 1.0e-8) :
        diff = True
        print ("Line %d: Numbers are different between the two files:" % (i+1))
        print ("Actual value = %e; Expected value = %e" % (val, exp_val))
    else :
      if (abs((val - exp_val) / exp_val) > 1.0e-8) :
        diff = True
        print ("Line %d: Numbers are different between the two files:" % (i+1))
        print ("Actual value = %f; Expected value = %f" % (val, exp_val))

  else :
    if (lines[i] != exp_lines[i]) :
      diff = True
      print ("Line %d: String is different between the two files." % (i+1))
      print ("  <File>         :", lines[i].rstrip())
      print ("  <Expected File>:", exp_lines[i].rstrip())


if (not diff) :
  print ("The two files are the same [with zero_tol = %e, and rel_tol = %e]" % (1.0e-8, 1.0e-8))
  sys.exit(0)
else :
  print ("The two files are different.")
  sys.exit(1)
