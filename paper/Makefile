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

DIFF_MAIN_FILE=causemann2026insilico-r2-diff
DIFF_SUPP_FILE=causemann2026insilico-r2-supplementary-diff
diff_main: 
	latexdiff --no-del --graphics-markup=none main_r1.tex main.tex > $(DIFF_MAIN_FILE).tex
	pdflatex -interaction=nonstopmode $(DIFF_MAIN_FILE).tex
	bibtex $(DIFF_MAIN_FILE)
	pdflatex -interaction=nonstopmode $(DIFF_MAIN_FILE).tex
	pdflatex -interaction=nonstopmode $(DIFF_MAIN_FILE).tex

diff_supp:
	latexdiff --no-del --graphics-markup=none supplementary_r1.tex supplementary.tex > $(DIFF_SUPP_FILE).tex
	pdflatex -interaction=nonstopmode $(DIFF_SUPP_FILE).tex
	bibtex $(DIFF_SUPP_FILE)
	pdflatex -interaction=nonstopmode $(DIFF_SUPP_FILE).tex
	pdflatex -interaction=nonstopmode $(DIFF_SUPP_FILE).tex

clean:
	rm *~ *.aux *.out *.toc *.blg *.bbl *.log *.tdo *.dvi 
