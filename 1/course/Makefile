
R_OPTS=--vanilla

.PHONY: main.tex

main.pdf: main.tex Makefile
	latexmk -pdf $<

main.tex: main.Rnw
	R ${R_OPTS} \
		-e "library(knitr)" \
		-e "knit('$<')"

view: main.pdf
	xdg-open $< > /dev/null 2>&1 &

view-guide: doc/ktad_kp1_2013.pdf
	xdg-open $< > /dev/null 2>&1 &

clean:
	latexmk -f -C main.pdf
	rm -f main.tex
	rm cache/*.RData
	rm cache/*.rdb
	rm cache/*.rdx
	rm figure/*
