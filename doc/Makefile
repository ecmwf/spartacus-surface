# If you have pdflatex installed, type "make" to create the
# documentation, "make clean" to delete it

BASE = spartacus_surface_documentation

documentation: $(BASE).pdf

$(BASE).pdf: $(BASE).tex surface_type_schematic.pdf isolation_schematic.pdf
	pdflatex $(BASE).tex
	pdflatex $(BASE).tex
	pdflatex $(BASE).tex

clean:
	rm -f $(BASE).pdf $(BASE).aux $(BASE).log $(BASE).out

.PHONY: documentation clean
