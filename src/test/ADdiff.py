#!/usr/bin/env python
""" usage: ADdiff.py [-h] [-t TOLERANCE] [-n NDATA] [filereference] [filenew]

Some description

positional arguments:
  filereference             name of the file with reference sensitivities
  filenew                   name of the file with new sensitivities

optional arguments:
  -h, --help            show this help message and exit
  -t TOLERANCE, --tolerance TOLERANCE
                        numerical differences below this treshold will be
                        ignored (default: 1e-13)
  -n NDATA, --ndata NDATA
                        number of sensitivity fields to compare
                        (default: 10)
Examples:
"""

from sys import exit
from sys import stdin
import argparse
import re

__author__ = 'Stefano Carli'
__email__ = 'stefano.carli@kuleuven.be'

tolerance = 1e-13 # default tolerance, can be overriden by the -t option
ndata = 10 # default number of sensitivity fields, can be overriden by the -n option

parser = argparse.ArgumentParser(description = "Some description",
            epilog=__doc__[1914:], 
            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("filereference", type=str, default="(stdin)", nargs='?',
           help="name of the file with reference sensitivities")
parser.add_argument("filenew", type=str, default="(stdin)", nargs='?',
           help="name of the file with new sensitivities")
parser.add_argument("-t", "--tolerance", type=float, default=tolerance,
           help="numerical differences below this treshold will be "
                "ignored (default: " + str(tolerance) + ")")
parser.add_argument("-n", "--ndata", type=int, default=ndata,
           help="number of sensitivity fields to compare"
                "(default: " + str(ndata) + ")")

args = parser.parse_args()
if args.tolerance:
    tolerance = args.tolerance
if args.ndata:
    ndata = args.ndata

## Collect reference data
# Open either stdin or the file given as argument
if args.filereference=='(stdin)':
    inputfile = stdin
else:
    try:
        inputfile = open(args.filereference, "rt")
    except IOError as err:
        print("Error opening file", args.filereference, err.strerror)
        exit(2)
   
reference = []
for line in inputfile:
    val = float(line)
    reference.append(val)

## Collect new data
if args.filenew=='(stdin)':
    inputfile = stdin
else:
    try:
        inputfile = open(args.filenew, "rt")
    except IOError as err:
        print("Error opening file", args.filenew, err.strerror)
        exit(2)

new = []
for line in inputfile:
    val = float(line)
    new.append(val)

correct = True

if (not len(reference)==ndata):
    print("Expected ",str(ndata)," sensitivities for reference, got instead ", str(len(reference)))
    exit(2)
if (not len(new)==ndata):
    print("Expected ",str(ndata)," sensitivities for new, got instead ", str(len(new)))
    exit(2)

# now check correctness
for ii in range(ndata):
    err = abs((reference[ii]-new[ii])/reference[ii])
    if (err>tolerance):
      correct = False
      print("WARNING! Sensitivity error for variable ",str(ii),":",str(err))

if correct:
    print("Results are correct")
    exit(0)
else:
    print("Results are NOT correct")
    exit(1)
