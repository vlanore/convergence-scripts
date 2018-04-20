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

## Trick to use input in python 2 et 3.*
try: input = raw_input
except NameError: pass

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
print("-- add transtion condition: "+param(add_transition))
pdf_window = args.pdf_window
if not pdf_window:
    print("-- Use a pop-up widow instead of a pdf file")
else:
    pdf_file = tree_file.name+".pdf"
    print("-- Use a pdf file instead of a pop-up widow  ("+param(pdf_file)+")")
out_file = tree_file.name+".annotated"
print("-- Output file is "+param(out_file))

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

def draw_tree(t):
    for n in t.traverse():
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
            if hasattr(n,"Transition"):
                nd.border.color = "red"
                nd.border.width = 2
        n.add_face(nd, column=0, position="float")
        n.add_face(TextFace("       "), column=0, position="branch-bottom")

draw_tree(t)

print("-- Starting subtree selection")
continue_flag = True
t_new = t.copy()
while continue_flag:
    if not pdf_window:
        print("-- Choose your node and close the window")
        t_new.show(tree_style=tree_style)
    else:
        print("-- Choose your node in "+data(pdf_file))
        t_new.render(pdf_file, tree_style=tree_style)
    print("(to save and quit type "+param("s")+")")
    testVar = input(ask_input("Please enter start of convergent subtree: "))
    if testVar.isdigit():
        t_new = t_new.copy("newick-extended")
        t_new.add_feature("i", t.i)
        t_new.add_feature("Condition", t.Condition)
        nb = int(testVar)
        print("-- Selected subtree rooted at node "+data(nb))
        print("  -- add tag Condition="+data(1)+" to the subtree of node " +data(nb))
        n_i = t_new.search_nodes(i=str(nb))
        if n_i:
            n_i = n_i[0]
            n_i.Condition = 1
            if add_transition:
                if hasattr(n_i,"Transition"):
                    n_i.Transition = 1
                else:
                    n_i.add_feature("Transition",1)
                print("  -- add tag Transition="+data(1)+" at node "+data(nb))
            for n_d in n_i.get_descendants():
                n_d.Condition = 1
            if sister:
                n_s = n_i.get_sisters()
                for n_s_i in n_s:
                    n_s_i.Condition = 2
                    if add_transition:
                        if hasattr(n_i,"Transition"):
                            n_i.Transition = 2
                        else:
                            n_i.add_feature("Transition",2)
                        print("  -- add tag Transition="+data(1)+" at node "+data(nb))
                    if add_transition:
                        n_s_i.Transition = 2
                    for n in n_s_i.get_descendants():
                        n.Condition = 2
        else:
            print(failure("Input node was not in the tree; try again"))
        draw_tree(t_new)
    elif testVar == "s":
        continue_flag = False
    else:
        print(failure("Input was not an integer; try again"))

#===================================================================================================
print(step("Writing result to file: "))
print("-- in: " + data(out_file))
features = ["Condition"]
if add_transition:
    features.append("Transition")
t_new.write(format=1,features=features, outfile = out_file)
