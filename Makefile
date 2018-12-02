##################################################################
# Makefile for LaTeX
##################################################################
TEX:=$(shell ls *.tex)
FILES = *.tex Makefile *.bst *.pdf *.bib *.sh
FOLDER =./home/saul/sauld@cimat.mx/UNISON/Ponencias/JXVInvitation/beamer_source_code

ARCHIVELIST = archive.txt
BIBTEXSRCS := $(shell sh biblio.sh)
BIBTEXSRCSLIST += $(BIBTEXSRCS)
OTHER = *~ *.aux *.dvi *.toc *.bbl *.blg *.out *.thm *.ps *.idx *.ilg *.ind \
*.tdo *.pdf *.tar.gz *.log *.spl *.synctex.gz
NAMETAR1:= $(shell date '+%Y%b%d_%H_%M')
NAMETAR = "$(NAMETAR1)-JXVinvitation-talk.tar.gz"
NAMETARTEX = "$(NAMETAR1)-JXVinvitation-talk.tar.gz"

lualatex: main.tex biblio.sh
	lualatex --synctex=1 -interaction=nonstopmode main.tex
	lualatex --synctex=1 -interaction=nonstopmode main.tex
	lualatex --synctex=1 -interaction=nonstopmode main.tex	
	sh biblio.sh
	lualatex --synctex=1 -interaction=nonstopmode main.tex
	lualatex --synctex=1 -interaction=nonstopmode main.tex

bib: biblio.sh
	sh biblio.sh

clean:
	rm -f $(OTHER)

tar:
	tar -cvf $(NAMETAR) -T $(ARCHIVELIST)

zip:
	tar -cvf $(NAMETAR) -T $(ARCHIVELIST)
