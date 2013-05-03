TEXS=$(wildcard *.tex)
TARGETS=$(addsuffix .pdf, $(basename $(TEXS)))

all: $(TARGETS)

%.pdf: %.tex
	pdflatex $<
	bibtex doc.aux
	pdflatex $<
	pdflatex $<
