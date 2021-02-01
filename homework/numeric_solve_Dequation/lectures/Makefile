default: notesAlgebraicTopology.pdf handoutsAlgebraicTopology.pdf

#	$(MAKE) --directory=pst
#	bibtex $*
#	pdflatex $<
#	pdflatex $<
%.pdf : %.tex sec
	pdflatex $<
	pdflatex $<

clean:
	$(RM) *.aux *.bbl *.blg *.log *.dvi *~ 

realclean:
	$(MAKE) clean
	$(RM) eps/*.pdf *.pdf */*~ *~
