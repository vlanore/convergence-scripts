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
from plot_data import *

import pandas as pd
from Bio import AlignIO, SeqIO, Seq, SeqRecord

import sys

#===================================================================================================
STEP("Parsing command line arguments")

from argparse import ArgumentParser, FileType
parser = ArgumentParser(description='''Tool to plot convergent sites from a result table, a tree and a MSA.
ex:\n\n
python script/plot_convergent_sites.py  -tsv example/tree4plot.tsv -msa example/tree4plot.fa -tree example/tree4plot.nw.annotated -out example/tree4plot.svg -meth Meth3,Meth1 -t Meth3:0.8,Meth1:70
''')
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
dic_threshold_by_method = {}
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

if not methods_to_be_plot:
    SUBMESSAGE(failure() + "No method to be plot. Check your table file and your option.")
    sys.exit(1)


SUBMESSAGE("Threshold to filter sites:")
methode_scales = {}
for meth in methods_to_be_plot:
    max_meth = max(df[meth])
    if max_meth > 1:
        methode_scales[meth] = 100
    else:
        methode_scales[meth] = 1
    if meth in dic_threshold_by_method:
        threshold = dic_threshold_by_method[meth]
    elif max_meth > 1:
        threshold = 99
        dic_threshold_by_method[meth] = threshold
    else:
        threshold = 0.99
        dic_threshold_by_method[meth] = threshold
    nb_sites = sum(df[meth]>threshold)
    print("    - %s: > %s (%s sites)" %(param(meth), data(threshold), nb_sites))

MESSAGE("read the msa")

try:
    alignment = AlignIO.read(msa_file, "fasta")
except Exception as exc:
    SUBMESSAGE(failure() + str(exc))
    sys.exit(1)

nb_seq=len(alignment)
nb_sites=len(alignment[0].seq)
SUBMESSAGE(success() + "There are %i sequences of %i sites" %(nb_seq, nb_sites))

seq_names = [ s.name for s in alignment]

MESSAGE("read the tree")

try:
    t=Tree(tree_file.name)
except Exception as exc:
    SUBMESSAGE(failure() + str(exc))
    sys.exit(1)

leaves_names = [ l.name for l in t.get_leaves()]
nb_leaves = len(leaves_names)

SUBMESSAGE(success() + "There are %i leaves" %(nb_leaves))

MESSAGE("Check link between sequence names and leaf names")

seq_names_not_in_leaves = list(set(seq_names) - set(leaves_names))
leaf_names_not_in_seq = list(set(leaves_names) - set(seq_names))
if seq_names_not_in_leaves:
    SUBMESSAGE(warning() + "These sequences are not in the tree:\n\t-%s" %("\n\t-".join(seq_names_not_in_leaves)))
if leaf_names_not_in_seq:
    SUBMESSAGE(warning() + "These leaves are not in the alignment:\n\t-%s" %("\n\t-".join(leaf_names_not_in_seq)))

if set(leaves_names) != set(seq_names) or len(leaves_names) != len(seq_names):
    if len(leaves_names) != len(seq_names):
        SUBMESSAGE(failure() + "There are not the same number of leaves and sequences")
    if set(leaves_names) != set(seq_names):
        SUBMESSAGE(failure() + "Sequence names and leaf names are not identical.")
    sys.exit(1)
else:
    SUBMESSAGE(success())


MESSAGE("Detect existing tags in the tree")
features = []
for n in t.traverse("postorder"): # get all features:
    features.extend(list(set(dir(n)) - set(dir(Tree()))))
features = list(set(features)) # list(set(*)) = remove duplicates

if "Condition" in features:
    SUBMESSAGE(success() + "Tag \"Condition\" detected")
else:
    SUBMESSAGE(warning() + "Tag \"Condition\" not detected")

if "Transition" in features:
    SUBMESSAGE(success() + "Tag \"Transition\" detected")



#===================================================================================================
STEP("Prepare plot")

def unlist(list2d):
    return [item for sublist in list2d for item in sublist]

def filter_l(l, pos):
    new_l = [None]*len(pos)
    for i in range(len(pos)):
        p = pos[i]
        new_l[i] = l[p-1]
    return new_l

SUBMESSAGE("Filter alignment")
# filter position:
bilan_f = {}
all_pos = range(1, nb_sites +1)
dict_pos_filtered = {}
for meth in methods_to_be_plot:
    dict_pos_filtered[meth] = df["Sites"][df[meth]>dic_threshold_by_method[meth]].tolist()

all_filtered_position = list(set(unlist(dict_pos_filtered.values())))
all_filtered_position.sort()
dict_pos_filtered["union"] = all_filtered_position

# filtered ali:
for meth in methods_to_be_plot + ["union"]:
    filtered_ali = []
    for seq in alignment:
        new_seq = SeqRecord.SeqRecord(Seq.Seq("".join(filter_l(list(seq.seq),dict_pos_filtered[meth]))), seq.id, "", "")
        filtered_ali.append(new_seq)
    SeqIO.write(filtered_ali, "filtered_ali."+meth+".faa", "fasta")
    if meth == "union":
        methstr = "union"
    else:
        methstr = meth

SUBMESSAGE("Conserved sited: %s" %(",".join(map(str,dict_pos_filtered ["union"]))))
#===================================================================================================
STEP("Draw plot")

positions_to_highlight = None

# all model
if dict_pos_filtered["union"]:
    dict_values_pcoc_filtered_model = {}
    for meth in methods_to_be_plot:
        dict_values_pcoc_filtered_model[meth] = df[meth][df["Sites"].isin(dict_pos_filtered["union"])].tolist()
    meth = "union"
    make_tree_ali_detect_combi(tree_file.name, "filtered_ali."+meth+".faa", out_file,
                               dict_benchmark = dict_values_pcoc_filtered_model,
                               x_values= dict_pos_filtered[meth], hp=positions_to_highlight,
                               methode_scales = methode_scales, methode_thresholds=dic_threshold_by_method)


#===================================================================================================
STEP("Writing result to file: ")
MESSAGE("-- Output file is " + data(out_file))
