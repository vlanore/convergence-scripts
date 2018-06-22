# Copyright or Copr. Centre National de la Recherche Scientifique (CNRS) (2018)
# Contributors:
# - Vincent Lanore <vincent.lanore@gmail.com>

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

"""Temporary script to run readdiffsel and extract a list of max convergence probability per site.

Usage:
  diffsel_analyze_result.py [options...] <chainname>

Positional arguments:
  chainname             the name of the diffsel chain to analyze (without extension)

Options:
  -h, --help            show this help message and exit
  -r, --readdiffsel-location <filename>
                        location of readdiffsel executable [default: _build/readdiffsel]
  -o, --output-file <filename>
                        output file [default: <chainname>.out]"""

from diffsel_script_utils import *

#===================================================================================================
print(step("Parsing command line arguments"))
from docopt import docopt

args = docopt(__doc__)
chainname = args["<chainname>"]
MESSAGE("Chain name is " + param(chainname))
out = args["--output-file"]
out = out[0] if out!=["<chainname>.out"] else chainname+".out"
MESSAGE("Output file is " + param(out))
rdexec = args["--readdiffsel-location"][0]
MESSAGE("Readdiffsel is here: "+param(rdexec))

#===================================================================================================
print(step("Running readdiffsel"))
from subprocess import Popen, PIPE

iterations = int(Popen("wc -l < " + chainname + ".trace", stdout=PIPE, shell=True).stdout.read())-1 # fragile
MESSAGE("Found trace with " + data(iterations) + " iterations")
burnin = int(iterations * 0.2)
MESSAGE("Burn-in set to " + data(burnin))

rdcommand = "%s -x %d 1 %d %s" % (rdexec, burnin, iterations, chainname)
MESSAGE("Command to run is " + data(rdcommand))
rdprocess = Popen(rdcommand, stdout=PIPE, stderr=PIPE, shell=True)
rdcode = rdprocess.wait()
if rdcode:
    MESSAGE(failure() + "Command produced an error")
    MESSAGE("stderr is:")
    print(rdprocess.stderr.read().decode('ascii'))
    exit(1)
else:
    MESSAGE(success() + "Command was a success")

#===================================================================================================
print(step("Reading readdiffsel results"))
import pandas as pd

meanfile = chainname + "_1.meandiffsel"
MESSAGE("Reading from " + data(meanfile))
rddata = pd.read_csv(meanfile, sep="\t", header=None, usecols=range(1,21))
MESSAGE("Data read contains %s rows and %s columns" % (data(rddata.shape[0]), data(rddata.shape[1])))
if rddata.shape[1] != 20:
    MESSAGE(failure() + "Expected 20 columns!")
    exit(1)
MESSAGE("Computing maximum probability per site")
maxes = rddata.max(axis=1)

#===================================================================================================
print(step("Writing result to file"))

MESSAGE("Output file is " + param(out))
maxes.index.name="Sites"
maxes.to_csv(out, sep='\t', header=["Diffsel"])
MESSAGE("Done :)")
