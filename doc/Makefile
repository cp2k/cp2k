manual.pdf: manual.tex manual.bib
	pdflatex manual.tex
	bibtex manual
	makeindex manual
	pdflatex manual.tex
	pdflatex manual.tex
clean:
	rm manual.ind manual.ilg manual.blg manual.bbl \
           manual.toc manual.out manual.log manual.idx manual.aux
