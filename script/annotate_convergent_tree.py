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

#===================================================================================================
STEP("Parsing command line arguments")

from argparse import ArgumentParser, FileType
parser = ArgumentParser(description='Annotates a phylogenetic trees with condition numbers. Output file uses the NHX format with the "Condition" tag.')
parser.add_argument('inputFile', metavar="input", type=FileType('r'), nargs=1, help='the tree file (newick format)')
parser.add_argument('-s', '--sister-branch-cond', dest="sister", action='store_true', help="toggle the use of a different condition for the sister branches of convergent branches")
parser.add_argument('-t', '--transition', dest="add_transition", action='store_true', help="add the tag Transition where you put a convergent event.")

args = parser.parse_args()

tree_file = args.inputFile[0]
MESSAGE("Sequence file is "+param(tree_file.name))
sister = args.sister
MESSAGE("Sister branch condition: "+param(sister))
add_transition = args.add_transition
MESSAGE("Add transtion condition: "+param(add_transition))
out_file = tree_file.name+".annotated"
MESSAGE("Output file is "+data(out_file))

#===================================================================================================
STEP("Setting tree and node styles")

condi_color_dic = {"0":"#E6E6FA", "1":"#ADD8E6", "2":"#90EE90"}

MESSAGE("Setting node style")
nstyle = NodeStyle()
nstyle["fgcolor"] = "black"
nstyle["size"] = 1

MESSAGE("Setting tree style")
tree_style = TreeStyle()
tree_style.show_leaf_name = False
tree_style.show_branch_length = False
tree_style.draw_guiding_lines = True
tree_style.complete_branch_lines_when_necessary = True
tree_style.legend_position = 1

MESSAGE("Setting legend with condition numbers and colors")
for condi_i in sorted(condi_color_dic.keys()):
    tf = TextFace("Condition      " + condi_i)
    tf.background.color = condi_color_dic[condi_i]
    tf.margin_right = 2
    tf.margin_top = 1
    tf.margin_left = 2
    tf.margin_bottom = 1
    tf.border.width = 1
    tree_style.legend.add_face(tf, column=1)

if add_transition:
    MESSAGE("Setting transition style")
    tf = TextFace("Transition -> x")
    tf.background.color = "white"
    tf.margin_right = 2
    tf.margin_top = 1
    tf.margin_left = 2
    tf.margin_bottom = 1
    tf.border.color = "red"
    tf.border.width = 2
    tree_style.legend.add_face(tf, column=1)

#===================================================================================================
STEP("Tree retrieval and preparation")

MESSAGE("Reading tree from file")
t = Tree(tree_file.name)

print("-- Detect existing tags")
# get all features:
features = []
for n in t.traverse("postorder"):
    features.extend(list(set(dir(n)) - set(dir(Tree()))))
features = list(set(features))

if not features:
    print("  * No detected tag")
    print("-- Setting all nodes to Condition = "+data(0))
else:
    print("  * detected tags: "+",".join(features))


if not "Condition" in features:
    print("-- Setting all nodes to Condition = "+data(0))
else:
    print("-- Setting all nodes without tag Condition to "+data(0))

def set_if_no_tag(node, tag, value):
    if not hasattr(node, tag):
        node.add_feature(tag, value)

for n in t.traverse("postorder"):
    set_if_no_tag(n,"Condition",0)


print("-- Numberings nodes")

i=0
for n in t.traverse("postorder"):
    n.add_feature("i",str(i))
    i+=1

#===================================================================================================
STEP("Convergent branch selection")

def draw_tree(tree):
    tree_copy = tree.copy("newick-extended")
    tree_copy.add_feature("i", tree.i)
    tree_copy.add_feature("Condition", tree.Condition)

    for n in tree_copy.traverse():
        if n.is_leaf():
                n.set_style(nstyle)
                n.add_face(TextFace(str(n.name)), column=0, position="aligned")
        else:
            n.set_style(nstyle)
        nd = TextFace(str(n.i))
        nd.background.color = condi_color_dic[str(n.Condition)]
        nd.margin_right = 2
        nd.margin_top = 1
        nd.margin_left = 2
        nd.margin_bottom = 1
        nd.border.width = 1
        if add_transition:
            if hasattr(n, "Transition"):
                nd.border.color = "red"
                nd.border.width = 2
        n.add_face(nd, column=0, position="float")
        n.add_face(TextFace("       "), column=0, position="branch-bottom")
    return tree_copy

def set_tag(node, tag, value):
    if hasattr(node, tag):
        setattr(node, tag, value)
    else:
        node.add_feature(tag, value)

def mark_subtree(node, condition):
    SUBMESSAGE("adding tag Condition = " + data(condition) + " to the subtree rooted at node " + data(node.i))
    node.Condition = condition
    for child in node.get_descendants():
        child.Condition = condition

    if add_transition:
        SUBMESSAGE("adding tag Transition = " + data(condition) + " at node " + data(node.i))
        set_tag(node, "Transition", condition)

draw_tree(t)

MESSAGE("Starting subtree selection")
pdf_file = tree_file.name+".pdf"
while True:
    # asking user input
    MESSAGE("Please look at "+data(pdf_file)+" to see node numbers")
    draw_tree(t).render(pdf_file, tree_style=tree_style)

    user_input = input(ask_input("Please enter start of convergent subtree (type "+green("s")+" to save and quit):"))

    #processing input
    if user_input.isdigit(): #if input is an int
        n_i = t.search_nodes(i=user_input) # locating subtree whose root is at nb
        if n_i: # if found
            n_i = n_i[0] # we expect only one result as indices are supposed to be unique

            mark_subtree(n_i, 1)

            if sister: # if sister is specified, handle sister trees
                n_s = n_i.get_sisters()
                for n_s_i in n_s:
                    mark_subtree(n_s_i, 2)

        draw_tree(t)

    elif user_input == "s":
        break
    else:
        print(failure("Input was not an integer; try again"))

#===================================================================================================
STEP("Writing result to file: ")

print("-- Output file is " + data(out_file))

if add_transition:
    features.append("Transition")
t.write(format=1, features=features, outfile = out_file)
