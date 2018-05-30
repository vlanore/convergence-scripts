# Copyright or Copr. Centre National de la Recherche Scientifique (CNRS) (2018)
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


from ete3 import Tree, NodeStyle, TreeStyle, TextFace
from diffsel_script_utils import *


import pandas as pd
import sys

#===================================================================================================
STEP("Parsing command line arguments")

from argparse import ArgumentParser, FileType
parser = ArgumentParser(description='Tool to plot convergent sites from a result table, a tree and a MSA.')
parser.add_argument('-tsv', metavar="table.tsv", type=FileType('r'),
    help='the result table (tabular file, with a column named "Sites" and the others named by the name of the convergent detection tool)', required=True)
parser.add_argument('-msa', metavar="msa.fa", type=FileType('r'),
    help='the msa file (fasta format)', required=True)
parser.add_argument('-tree', metavar="tree.nhx", type=FileType('r'),
    help='the tree file (NHX format with the "Condition" tag)', required=True)
parser.add_argument('-out', metavar="output.svg", type=str,
    help='the output file (format svg)', required=True)

parser.add_argument('-meth', dest="methods_to_be_plot", type=str,
    metavar="\"meth1,meth3\"", help="columns of the table file to use (default:all)",
    default=None)
parser.add_argument('-t', dest="threshold_by_method", type=str,
    metavar="\"meth1:0.85,meth3:70\"", help="Threshold to filter site by method in the table file (default:0.99 or 99)",
    default=None)

args = parser.parse_args()

tsv_file = args.tsv
MESSAGE("Tabular file is "+param(tsv_file.name))
msa_file = args.msa
MESSAGE("MSA file is     "+param(msa_file.name))
tree_file = args.tree
MESSAGE("Tree file is    "+param(tree_file.name))
out_file = args.out
MESSAGE("Output file is  "+param(out_file))

methods_to_be_plot = args.methods_to_be_plot
if methods_to_be_plot:
    methods_to_be_plot = methods_to_be_plot.split(",")
    MESSAGE("methods_to_be_plot: "+param(methods_to_be_plot))
threshold_by_method = args.threshold_by_method
if threshold_by_method:
    dic_threshold_by_method = {mt.split(":")[0]:float(mt.split(":")[1]) for mt in threshold_by_method.split(",") }
    MESSAGE("threshold_by_method: "+param(dic_threshold_by_method))


#===================================================================================================
STEP("Read input files:")
MESSAGE("read the table")

df = pd.read_table(tsv_file)
colnames = list(df)

if "Sites" in colnames:
    SUBMESSAGE(success() + "\"Sites\" column is present")
    colnames.remove("Sites")
else:
    SUBMESSAGE(failure() + "\"Sites\" column is not present in the input table")
    sys.exit(1)

SUBMESSAGE("Detected methods: " + ", ".join([param(m) for m in colnames]))

if methods_to_be_plot:
    methods_to_be_plot = list(set(colnames) & set(methods_to_be_plot))
    SUBMESSAGE("Methods which will be used after filtering: " + ", ".join([param(m) for m in methods_to_be_plot]))
else:
    methods_to_be_plot = colnames
    SUBMESSAGE("All methods will be used")

SUBMESSAGE("Threshold to filter sites:")
for meth in methods_to_be_plot:
    max_meth = max(df[meth])
    if meth in dic_threshold_by_method:
        threshold = dic_threshold_by_method[meth]
    elif max_meth > 1:
        threshold = 99
    else:
        threshold = 0.99
    nb_sites = sum(df[meth]>threshold)
    print("    - %s: > %s (%s sites)" %(param(meth), data(threshold), nb_sites))



MESSAGE("read the msa")


MESSAGE("read the tree")


#===================================================================================================
STEP("Draw plot")

#===================================================================================================
STEP("Writing result to file: ")
MESSAGE("-- Output file is " + data(out_file))
