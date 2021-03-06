# Convergence scripts

[![licence CeCILL](https://img.shields.io/badge/license-CeCILL--C-blue.svg)](http://www.cecill.info/licences.en.html)

This is a collection of scripts to perform pre and post processing of data for convergence detection tools.

Notable scripts include:

 * [annotate_diffsel.py](script/annotate_diffsel.py), a tool to annotate phylogenetic trees with conditions (typically, convergent or non-convergent). It uses the graphical interface of ete3 to help annotate the trees. Currently, the output format is that of [diffsel](https://github.com/vlanore/diffsel) (i.e., it encodes conditions as branch lengths).

 * [annotate_convergent_tree.py](script/annotate_convergent_tree.py), a tool to annotate phylogenetic trees with conditions (typically, convergent or non-convergent). It uses the graphical interface of ete3 to help annotate the trees or a pdf file. Currently, the output format is a Nhx tree with the tag "Condition" ( [&&NHX:Condition=0] means non-convergent and [&&NHX:Condition=1] means convergent).
