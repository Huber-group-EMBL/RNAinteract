all: RNAinteract clean

RNAinteract: RNAinteract.tex
	texi2dvi --pdf RNAinteract.tex
	R CMD Stangle RNAinteract.Rnw

clean:
	-$(RM) -rf RNAinteractHTML
	-$(RM) -f RNAinteract-*
	-$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	-$(RM) -f Rplots.pdf

