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
	latexdiff main_r0.tex main.tex > diff_main_r0.tex
	pdflatex diff_main_r0.tex
	bibtex diff_main_r0
	pdflatex diff_main_r0.tex
	pdflatex diff_main_r0.tex

	latexdiff supplementary_r0.tex main.tex > diff_supplementary_r0.tex
	pdflatex diff_supplementary_r0.tex
	bibtex diff_supplementary_r0
	pdflatex diff_supplementary_r0.tex
	pdflatex diff_supplementary_r0.tex

clean:
	rm *~ *.aux *.out *.toc *.blg *.bbl *.log *.tdo *.dvi
