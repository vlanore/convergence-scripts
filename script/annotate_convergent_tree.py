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
print(step("Parsing command line arguments"))

from argparse import ArgumentParser, FileType
parser = ArgumentParser(description='Annotates a phylogenetic trees with condition numbers compatible with diffsel.')
parser.add_argument('inputFile', metavar="input", type=FileType('r'), nargs=1, help='the tree file (newick format)')
parser.add_argument('-s', '--sister-branch-cond', dest="sister", action='store_true', help="toggle the use of a different condition for the sister branches of convergent branches")
parser.add_argument('--pdf', dest="pdf_window", action='store_true', help="use a pdf file instead of a interactive window.")
parser.add_argument('--transition', dest="add_transition", action='store_true', help="add the tag Transition where you put a convergent event.")

args = parser.parse_args()

tree_file = args.inputFile[0]
print("-- Sequence file is "+param(tree_file.name))
sister = args.sister
print("-- Sister branch condition: "+param(sister))
add_transition = args.add_transition
print("-- Add transtion condition: "+param(add_transition))
pdf_window = args.pdf_window
if not pdf_window:
    print("-- Use a pop-up widow instead of a pdf file")
else:
    pdf_file = tree_file.name+".pdf"
    print("-- Use a pdf file instead of a pop-up widow  ("+data(pdf_file)+")")
out_file = tree_file.name+".annotated"
print("-- Output file is "+data(out_file))

#===================================================================================================
print(step("Setting tree and node styles"))

condi_color_dic = {"0":"#E6E6FA", "1":"#ADD8E6", "2":"#90EE90"}

print("-- Setting node styles")
nstyle = NodeStyle()
nstyle["fgcolor"] = "black"
nstyle["size"] = 1

nstyle_L = NodeStyle() # VL: isn't it identical to nstyle?
nstyle_L["fgcolor"] = "black"
nstyle_L["size"] = 1

print("-- Setting tree style")
tree_style = TreeStyle()
tree_style.show_leaf_name = False
tree_style.show_branch_length = False
tree_style.draw_guiding_lines = True
tree_style.complete_branch_lines_when_necessary = True
tree_style.legend_position = 1

print("-- Setting legend with condition numbers and colors")
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
    print("-- Setting transition style")
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
print(step("Tree retrieval and preparation"))

print("-- Reading tree from file")
t = Tree(tree_file.name)

print("-- Setting all node  to Condition = "+data(0))
i=0
for n in t.traverse("postorder"):
    n.add_feature("Condition",0)
    n.add_feature("i",str(i))
    i+=1

#===================================================================================================
print(step("Convergent branch selection"))

def draw_tree(tree):
    tree_copy = tree.copy("newick-extended")
    tree_copy.add_feature("i", tree.i)
    tree_copy.add_feature("Condition", tree.Condition)

    for n in tree_copy.traverse():
        if n.is_leaf():
                n.set_style(nstyle_L)
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
    if hasattr(node, "Transition"):
        setattr(node, tag, value)
    else:
        node.add_feature(tag, value)

def mark_subtree(node, condition):
    print("  * adding tag " + data("Condition="+str(condition)) + " to the subtree rooted at node " + data(node.i))
    node.Condition = condition
    for child in node.get_descendants():
        child.Condition = condition

    if add_transition:
        print("  * adding tag " + data("Transition="+str(condition)) + " at node " + data(node.i))
        set_tag(node, "Transition", condition)

draw_tree(t)

print("-- Starting subtree selection")
continue_flag = True
new_tree = t.copy()
while continue_flag:
    # asking user input
    if not pdf_window:
        print("-- Choose your node and close the window")
        draw_tree(new_tree).show(tree_style=tree_style)
    else:
        print("-- Please look at "+data(pdf_file)+" to see node numbers")
        draw_tree(new_tree).render(pdf_file, tree_style=tree_style)

    user_input = input(ask_input("Please enter start of convergent subtree (type "+green("s")+" to save and quit):"))

    #processing input
    if user_input.isdigit(): #if input is an int
        n_i = new_tree.search_nodes(i=user_input) # locating subtree whose root is at nb
        if n_i: # if found
            n_i = n_i[0] # we expect only one result as indices are supposed to be unique

            mark_subtree(n_i, 1)

            if sister: # if sister is specified, handle sister trees
                n_s = n_i.get_sisters()
                for n_s_i in n_s:
                    mark_subtree(n_s_i, 2)

        draw_tree(new_tree)

    elif user_input == "s":
        continue_flag = False
    else:
        print(failure("Input was not an integer; try again"))

#===================================================================================================
print(step("Writing result to file: "))

print("-- Output file is " + data(out_file))
features = ["Condition"]
if add_transition:
    features.append("Transition")
new_tree.write(format=1, features=features, outfile = out_file)
