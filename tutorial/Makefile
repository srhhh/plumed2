tutorial.pdf : tutorial.ind tutorial.bbl
	cd figures ; make
	pdflatex tutorial
tutorial.ind : tutorial.idx tutorial.toc
	latex tutorial
	makeindex tutorial
tutorial.idx : *.tex
	latex tutorial
tutorial.bbl : *.tex
	bibtex tutorial
	latex tutorial
tutorial.toc : *.tex
	latex tutorial
clean :
	rm  tutorial.dvi tutorial.aux  tutorial.idx  tutorial.ilg  tutorial.ind  tutorial.log  tutorial.out  tutorial.pdf tutorial.toc  tutorial.bbl  tutorial.blg
