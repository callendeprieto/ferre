
all: ferre.pdf

clean:
	rm ferre.aux ferre.dvi ferre.log ferre.toc ferre.pdf

ferre.pdf: ferre.tex f1.ps
	latex ferre.tex
	latex ferre.tex
	dvipdf ferre
	rm ferre.aux ferre.dvi ferre.log ferre.toc

