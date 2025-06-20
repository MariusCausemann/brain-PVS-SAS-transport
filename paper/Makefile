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
	latexdiff --graphics-markup=none main_r0.tex main.tex > diff_main_r0.tex
	pdflatex -interaction=nonstopmode diff_main_r0.tex
	bibtex diff_main_r0
	pdflatex -interaction=nonstopmode diff_main_r0.tex
	pdflatex -interaction=nonstopmode diff_main_r0.tex

diff_supp:
	latexdiff --graphics-markup=none supplementary_r0.tex supplementary.tex > diff_supplementary_r0.tex
	pdflatex -interaction=nonstopmode diff_supplementary_r0.tex
	bibtex diff_supplementary_r0
	pdflatex -interaction=nonstopmode diff_supplementary_r0.tex
	pdflatex -interaction=nonstopmode diff_supplementary_r0.tex

clean:
	rm *~ *.aux *.out *.toc *.blg *.bbl *.log *.tdo *.dvi
