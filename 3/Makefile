
R_OPTS=--vanilla

main.pdf: main.tex main.aux 
	pdflatex $<

main.tex: main.Rnw
	R ${R_OPTS} -e "library(knitr);knit('$<')"

view: main.pdf
	xdg-open $< > /dev/null 2>&1 &

%.aux: %.tex
	pdflatex $<

clean:
	rm -f *.pdf
	rm -f *.aux
	rm -f *.log
	rm -f *.toc
	rm -f *.out
