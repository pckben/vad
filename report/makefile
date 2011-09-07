MAIN = main
BIBFILE = bibliography.bib
OUTDIR  = ./output
LATEX_ARGS = --interaction=nonstopmode -file-line-error -halt-on-error -output-directory=$(OUTDIR) 

pdf : dvi
	dvips $(OUTDIR)/$(MAIN).dvi -o $(OUTDIR)/$(MAIN).ps
	ps2pdf $(OUTDIR)/$(MAIN).ps
	#pdflatex $(LATEX_ARGS) $(MAIN).tex 

dvi : *.tex *.bib
	! test $(OUTDIR) && mkdir $(OUTDIR) || echo
	latex $(LATEX_ARGS) $(MAIN).tex
	cp $(BIBFILE) $(OUTDIR)
	cd $(OUTDIR) ; bibtex $(MAIN)
	latex $(LATEX_ARGS) $(MAIN).tex

pdf_display: pdf
	open -a /Applications/Skim.app output/$(MAIN).pdf

