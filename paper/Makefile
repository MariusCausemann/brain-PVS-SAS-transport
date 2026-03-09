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

diff: diff_main diff_supp

diff_main: 
	latexdiff --subtype=INVISIBLE --graphics-markup=none main_r1.tex main.tex > diff_main_r1.tex
	pdflatex -interaction=nonstopmode diff_main_r1.tex
	bibtex diff_main_r1
	pdflatex -interaction=nonstopmode diff_main_r1.tex
	pdflatex -interaction=nonstopmode diff_main_r1.tex

diff_supp:
	latexdiff --graphics-markup=none supplementary_r1.tex supplementary.tex > diff_supplementary_r1.tex
	pdflatex -interaction=nonstopmode diff_supplementary_r1.tex
	bibtex diff_supplementary_r1
	pdflatex -interaction=nonstopmode diff_supplementary_r1.tex
	pdflatex -interaction=nonstopmode diff_supplementary_r1.tex

clean:
	rm *~ *.aux *.out *.toc *.blg *.bbl *.log *.tdo *.dvi
