name=prezentace_Myslivec_Vojtech_2016
target=$(name).pdf
source=$(name).tex

.PHONY: all
all: $(target)

$(target): $(source)
	vlna $^
	pdflatex $^
	pdflatex $^

.PHONY: clean
clean:
	rm -f $(name).aux $(name).loa $(name).lof $(name).log $(name).lot $(name).nav $(name).out $(name).snm $(name).te~ $(name).toc

