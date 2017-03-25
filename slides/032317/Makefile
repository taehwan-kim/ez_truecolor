all: lmkv

PDF       = open -a Preview
PDFLATEX  = pdflatex
LATEXARGS = -interaction batchmode
BIBTEX    = bibtex
BIBARGS   =
IMAGEDIR  = ../images/
LATEXMK   = latexmk

#BIB       = phase1
PROJ      = qual_taehwan

#%.aux: %.tex $(BIB).bib
	#-$(PDFLATEX) $(LATEXARGS) $^
	#-$(BIBTEX) $(BIBARGS) $*.aux

%.pdf: %.tex
	$(LATEXMK) -pdf -bibtex-cond --file-line-error -halt-on-error $^

	#-$(PDFLATEX) $(LATEXARGS) $^
	#-$(PDFLATEX) $(LATEXARGS) $^

#$(IMAGEDIR)/.timestamp:
	#make -C ../images

view: $(PROJ).pdf
	$(PDF) $< &

lmk: $(PROJ).pdf

lmkv:
	make lmk
	make view

#lmk: 

.PHONY: lmk view

clean:
	rm -f *.pdf *.aux *.bbl *.log *.blg *.xml *-blx.bib *.toc *.vrb *.snm *.fls *.nav *.out *.fdb_latexmk
