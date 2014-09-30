default: asydo-adass.pdf 

clean:
	rm -f *~ *.aux *.bbl *.blg *.lof *.log *.lot *.toc *.nav *.snm *.idx *.out 

distclean: clean
	rm -f *.pdf *.dvi

.PHONY: all clean distclean

%.pdf: %.tex
	latex -halt-on-error $<
	bibtex $(basename $<)
	latex -halt-on-error $<
	latex -halt-on-error $<
