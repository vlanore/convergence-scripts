.PHONY: annot plot

# set env variable PYTHON=python3 to use python3
PY=$(shell echo $${PYTHON-python})

annot:
	$(PY) script/annotate_convergent_tree.py example/newick.tree

plot:
	$(PY) script/plot_convergent_sites.py -tsv example/tree4plot.tsv -msa example/tree4plot.fa -tree example/tree4plot_annotated.nw -out example/tree4plot.svg -meth Meth3,Meth1 -t Meth3:0.8,Meth1:70
