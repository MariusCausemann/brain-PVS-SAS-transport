FILE=main

default:
	pdflatex $(FILE).tex

final:
	pdflatex $(FILE).tex
	bibtex $(FILE)
	pdflatex $(FILE).tex
	pdflatex $(FILE).tex

	pdflatex supplementary.tex
	bibtex supplementary
	pdflatex supplementary.tex
	pdflatex supplementary.tex

diff:
	latexdiff old.tex paper.tex > diff.tex
	pdflatex diff.tex
	bibtex diff
	pdflatex diff.tex
	pdflatex diff.tex

clean:
	rm *~ *.aux *.out *.toc *.blg *.bbl *.log *.tdo *.dvi
