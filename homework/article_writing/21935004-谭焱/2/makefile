
run:
	xelatex main.tex
	bibtex main.aux
	xelatex main.tex
	xelatex main.tex
	make clear

clear:
	rm -f main.abs main.aux main.bbl main.blg main.log main.out main.synctex.gz

clean:
	make clear
	rm -f main.pdf