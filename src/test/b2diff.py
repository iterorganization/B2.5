#!/usr/bin/env python
""" usage: b2diff.py [-h] [-t TOLERANCE] [-s SPECIFIC_TOLERANCE] [-i IGNORE]
                 [-n NOUT] [-l LIST_ALWAYS] [-m MAXERR] [-a AVGERR] [-v]
                 [filename]

Post process the results from check_b2_output, and list the variables with largest error. Returns with exit status 1 if errors are larger than tolerance, otherwise returns 0. By default the value (max_error/avg_value) is compared to the tolerance.

positional arguments:
  filename              name of the input files

optional arguments:
  -h, --help            show this help message and exit
  -t TOLERANCE, --tolerance TOLERANCE
                        numerical differences below this treshold will be
                        ignored (default: 1e-13)
  -s SPECIFIC_TOLERANCE, --specific-tolerance SPECIFIC_TOLERANCE
                        you can set a specific tolerance for certain
                        variables, like -s 'var1 1e-7 var2 1e-9'
  -i IGNORE, --ignore IGNORE
                        regex pattern to ignore certain variables
  -n NOUT, --nout NOUT  number of variables to list
  -l LIST_ALWAYS, --list-always LIST_ALWAYS
                        List of variable names that should always be listed.
                        The variables 'te ti po ne ni ua na' are included by
                        default
  -m MAXERR, --maxerr MAXERR
                        List of variables for which the max relative error
                        should be checked against tolerance
  -a AVGERR, --avgerr AVGERR
                        List of variables for which the average error should
                        be checked against tolerance
  -v, --verbose         More verbose output. Use -v to get the name of the
                        variables that are not correct. Use -vv to list the
                        ignored variables too. The input parameters are also
                        listed if -vvv is specified.

Examples:
  - Ignore variables del* and res*:
    ./b2diff.py --ignore 'del*|res*' input.txt
    
  - Use different tolerances: 1e-7 for ua, and 1e-10 for fch
    ./b2diff.py --specific-tolerance 'ua 1e-7 fch 1e-10' input.txt
    
  - Check max relative error for ni and te, also use tolerance 1e-12 for them:
    ./b2diff.py --maxerr 'ni te' --specific-tolerance 'ni 1e-12 te 1e-12'
"""

from sys import exit
from sys import stdin
import argparse
import re
__author__ = 'Tamas Feher'
__email__ = 'tamas.bela.feher@ipp.mpg.de'

tolerance = 1e-13 # default tolerance, can be overriden by the -t option

parser = argparse.ArgumentParser(description = "Post process the results from "
            "check_b2_output, and list the variables with largest error. "
            "Returns with exit status 1 if errors are larger than tolerance, "
            "otherwise returns 0. By default the value (max_error/avg_value) "
            "is compared to the tolerance.",
            epilog=__doc__[1914:], 
            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("filename", type=str, default="(stdin)", nargs='?',
           help="name of the input files")
parser.add_argument("-t", "--tolerance", type=float, default=tolerance,
           help="numerical differences below this threshold will be "
                "ignored (default: " + str(tolerance) + ")")
parser.add_argument("-s", "--specific-tolerance", type=str, help="you can set "
           "a specific tolerance for certain variables, like -s 'var1 1e-7 "
           "var2 1e-9'")
parser.add_argument("-i", "--ignore", type=str, help="regex pattern to ignore " 
                    "certain variables ")
parser.add_argument("-n", "--nout", type=int, default=10, 
           help="number of variables to list")
parser.add_argument("-l", "--list-always", type=str, help="List of variable "
           "names that should always be listed. The variables 'te ti po ne ni "
           "ua na' are included by default")                    
parser.add_argument("-m", "--maxerr", type=str, help="List of variables for "
              "which the max relative error should be checked against "
              "tolerance")
parser.add_argument("-a", "--avgerr", type=str, help="List of variables for "
                    "which the average error should be checked against tolerance")
parser.add_argument("-v", "--verbose", action='count', help="More verbose output. "
                    "Use -v to get the name of the variables that are not correct. "
                    "Use -vv to list the ignored variables too. The input "
                    "parameters are also listed if -vvv is specified.")


args = parser.parse_args()
if args.tolerance:
    tolerance = args.tolerance

# The error in these arrays will be always printed
needit = {'te', 'ti', 'po', 'ne', 'ni', 'ua', 'na'}
if args.list_always:
    needit |= set(args.list_always.split())
if args.verbose > 2: print('Variables always printed (except if error is 0):', " ".join(needit))
    
spec_tol = dict()
if args.specific_tolerance:
    tlist = args.specific_tolerance.split()
    for i in range(0,len(tlist),2):
        spec_tol[tlist[i]] = float(tlist[i+1])
    if args.verbose > 2: print('Specific tolerance', spec_tol)

if args.maxerr:
    use_maxerr = args.maxerr.split()
    if args.verbose > 2: print('using maxerr for', use_maxerr)
else:
    use_maxerr = list()

if args.avgerr:
    use_avgerr = args.avgerr.split()
    if args.verbose > 2: print('using average error for', use_avgerr)
else:
    use_avgerr = list()
        
if args.ignore:
    ignore = re.compile(args.ignore)
else:
    ignore = None
    
# Open either stdin or the file given as argument
if args.filename=='(stdin)':
    inputfile = stdin
else:
    try:
        inputfile = open(args.filename, "rt")
    except IOError as err:
        print("Error opening file", args.filename, err.strerror)
        exit(2)
    
# We collect the results form check_b2_outptut into a dictionary, where the key
# is the variable name, and the value is a list of
# [Max_relative_err, max_absolute_err, avg_relative_err, avg_array_value, maxerr/avgval]
properties = ['Max absolute error in ', 'Average relative error in',
              'Average array value (abs) in', 'Max error / avg array value']
# note that property 'Max error in' is handled separately 

errors = dict()

for line in inputfile:
    if line.isspace(): continue
    if line.find('Max error in') == 0:
        name = line.split()[3]
        val = float(line.split()[4])
        if (not name in errors) or (val > errors[name][0]):
            errors[name] = [val]
    # We will append the other values to the list errors[name]
    for p in properties:
        if line.find(p) == 0:
            n = len(p.split())
            name = line.split()[n]
            val = float(line.split()[n+1])
            if len(errors[name]) < len(properties)+1:
                # only append if the list is not yet full
                errors[name].append(val)

print_these = list()
# Add the arrays with the largest max errors to the print_these dictionary
# also check the correctness
N = args.nout
i = 0
correct = True
for x in sorted(list(errors.items()), key=lambda elem: elem[1][0], reverse=True):
    # x[0] is the key of the dictionary == variable name
    # x[1] is the list of properties that we collected [max err, max_abs_err, etc]
    if ignore and ignore.match(x[0]):
        if args.verbose>1:
            print('Ignoring variable', x[0])
    else:
        i += 1
        if i <= N or x[0] in needit:
            # every element of print_these will be a flat list of [name, max_err, max_abs_err, ...]
            print_these.append([x[0]] + x[1])
        if x[0] in use_maxerr:
            val = x[1][0] # max relative error
        elif x[0] in use_avgerr:
            val = x[1][2] # average relative error
        else:
            val = x[1][4] # max error / avg array value
        if x[0] in spec_tol:
            tol = spec_tol[x[0]]
        else:
            tol = tolerance
        if val > tol:
            correct = False
            if args.verbose > 0:
                print('Error in', x[0], '(' + str(val) + ')', 
                      'is larger than tolerance (' + str(tol) + ')')

print("# name     max_rel_err   maxerr  avg.rel.err  avg.val  maxerr/avg.val")
for y in print_these:
    print("{0:12s}  {1: 7.1e} {2: 7.1e}   {3: 7.1e}   {4: 7.1e} {5: 7.1e}".format(y[0], y[1], y[2], y[3], y[4], y[5]))

if correct:
    if args.verbose >0:
        print("Results are correct")
    exit(0)
else:
    if args.verbose >0:
        print("Results are NOT correct")
    exit(1)
