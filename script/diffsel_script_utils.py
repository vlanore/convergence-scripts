# Copyright or Copr. Centre National de la Recherche Scientifique (CNRS) and École Normale Supérieure de Lyon (2018)
# Contributors:
# - Vincent Lanore <vincent.lanore@gmail.com>
# - Carine Rey <carine.rey@ens-lyon.org>

# This software is a computer program whose purpose is to provide a set of scripts for pre and post processing of data for
# convergence detection programs.

# This software is governed by the CeCILL-C license under French law and abiding by the rules of distribution of free software.
# You can use, modify and/ or redistribute the software under the terms of the CeCILL-C license as circulated by CEA, CNRS and
# INRIA at the following URL "http://www.cecill.info".

# As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users
# are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive
# licensors have only limited liability.

# In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or
# reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated
# to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth
# computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements
# in conditions enabling the security of their systems and/or data to be ensured and, more generally, to use and operate it in
# the same conditions as regards security.

# The fact that you are presently reading this means that you have had knowledge of the CeCILL-C license and that you accept
# its terms.

import sys
import random

## Trick to use input in python 2 et 3.*
try: input = raw_input
except NameError: pass

# Color-related functions
if sys.stdout.isatty():
    class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[34m'
        OKGREEN = '\033[32m'
        YELLOW = '\033[33m'
        WARNING = '\033[93m'
        CYAN = '\033[33m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'
else:
    class bcolors:
        HEADER = ''
        OKBLUE = ''
        OKGREEN = ''
        WARNING = ''
        CYAN = ''
        FAIL = ''
        ENDC = ''
        BOLD = ''
        UNDERLINE = ''

def boldred(string):
    return bcolors.FAIL+bcolors.BOLD+string+bcolors.ENDC

def red(string):
    return bcolors.FAIL+string+bcolors.ENDC

def yellow(string):
    return bcolors.YELLOW+string+bcolors.ENDC

def boldcyan(string):
    return bcolors.CYAN+bcolors.BOLD+string+bcolors.ENDC

def param(myparam):
    return bcolors.OKBLUE+str(myparam)+bcolors.ENDC

def data(myparam):
    return bcolors.CYAN+str(myparam)+bcolors.ENDC

def step(string):
    return bcolors.BOLD+bcolors.HEADER+string+bcolors.ENDC

def boldgreen(string):
    return bcolors.BOLD+bcolors.OKGREEN+string+bcolors.ENDC

def green(string):
    return bcolors.OKGREEN+string+bcolors.ENDC

def good(string):
    return "-- ("+bcolors.OKGREEN+"Good"+bcolors.ENDC+") "+string

def bad(string):
    return "-- ("+bcolors.FAIL+"Bad"+bcolors.ENDC+") "+string

def success_str(string):
    return "-- ["+boldgreen("SUCCESS")+"] "+string

def success():
    return "["+boldgreen("SUCCESS")+"] "

def failure_str(string):
    return "-- ["+boldred("FAILURE")+"] "+string

def failure():
    return "["+boldred("FAILURE")+"] "

def warning_str(string):
    return "-- ["+boldred("WARNING")+"] "+string

def warning():
    return "["+boldred("WARNING")+"] "

def ask_input(string):
    return "-- ["+boldcyan("INPUT")+"] "+str(string)

def STEP(string):
    print(step(string))

def MESSAGE(string):
    print("-- "+string)

def SUBMESSAGE(string):
    print("  * "+string)

def WARNING(string):
    print(warning_str(string))
