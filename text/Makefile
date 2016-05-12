name=DP_Myslivec_Vojtech_2016
target=$(name).pdf
source=$(name).tex
tmp=tmp.pdf

.PHONY: all
all: $(target)

.PHONY: ghost
ghost: $(target)
	gs -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -sOutputFile=$(tmp) $^
	mv $(tmp) $^

$(target): $(source)
	vlna $^
	pdflatex $^
	pdflatex $^
	pdflatex $^


.PHONY: clean
clean:
	rm -f $(name).aux $(name).loa $(name).lof $(name).log $(name).lot $(name).out $(name).te~ $(name).toc

