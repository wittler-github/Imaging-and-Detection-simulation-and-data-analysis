

FILENAME=nadia

#my_paper

all:	simple
	acroread $(FILENAME).pdf

simple: $(FILENAME).dvi
	dvipdf $(FILENAME).dvi
#	dvips $(FILENAME).dvi


$(FILENAME).pdf: $(FILENAME).dvi
	dvipdf $(FILENAME).dvi

$(FILENAME).dvi: $(FILENAME).tex #$(FILENAME).bib
	latex $(FILENAME).tex
	latex $(FILENAME).tex
#	latex $(FILENAME).tex
#	bibtex $(FILENAME)



clean:
	rm -f $(FILENAME).dvi
	rm -f *~ *.aux *.log *.toc *.bak *.out *.blg *.bbl

clobber: clean
	rm -f $(FILENAME).ps $(FILENAME).pdf
