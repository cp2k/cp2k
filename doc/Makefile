SHELL = /bin/sh
.SUFFIXES:

ROBODOC=/usr/local/bin/robodoc
ROBOOPTS=SORT 

#
#
CP2KHOME=$(HOME)/CP2K/
VPATH=$(CP2KHOME)/src

# Your source files.
#
include $(CP2KHOME)/src/OBJECTDEFS
SOURCES=$(OBJECTS_GENERIC:.o=.F)

# The name of your Project
#
PROJECT=cp2k

# The various documentation files, derived from the source files.
# HTML
#
HTMLDOCS=$(SOURCES:=.html)
HTMLXREFS=$(HTMLDOCS:.html=.html.xref)
HTMLXREFSFILE=$(PROJECT)_html.xrefs
# LATEX
#
LATEXDOCS=$(SOURCES:=.tex)
LATEXXREFS=$(LATEXDOCS:.tex=.tex.xref)
LATEXXREFSFILE=$(PROJECT)_tex.xrefs
# ASCII
#
ASCIIDOCS=$(SOURCES:=.txt)
# RTF
#
RTFDOCS=$(SOURCES:=.rtf)
RTFXREFS=$(RTFDOCS:.rtf=.rtf.xref)
RTFXREFSFILE=$(PROJECT)_rtf.xrefs

# Some common targets
xrefall: xhtml xtex xrtf
docall: html tex ascii rtf

# Create the xref files for the various formats.
xhtml: $(HTMLXREFSFILE) 
xtex: $(LATEXXREFSFILE) 
xrtf: $(RTFXREFSFILE)

# Create the documentation files for the various formats.
html: $(HTMLDOCS) $(PROJECT)_mi.html 
tex: $(LATEXDOCS) $(PROJECT)_mi.tex
rtf: $(RTFDOCS)
ascii: $(ASCIIDOCS)

# master index file, currently works only for html and latex documentation.
# Note that you can define the title of the document.
$(PROJECT)_mi.html: $(HTMLXREFSFILE) 
	$(ROBODOC) $< $@ INDEX HTML TITLE "$(PROJECT) Master Index"

$(PROJECT)_mi.tex: $(LATEXXREFSFILE)
	$(ROBODOC) $< $@ INDEX LATEX TITLE "$(PROJECT) API Reference"

# create xrefs file (file with the names of all .xref files).
$(HTMLXREFSFILE) : $(HTMLXREFS)
	/bin/ls $(HTMLXREFS) > $@
$(LATEXXREFSFILE) : $(LATEXXREFS)
	/bin/ls  $(LATEXXREFS) > $@
$(RTFXREFSFILE) : $(RTFXREFS)
	/bin/ls  $(RTFXREFS) > $@

# Rule to create an .xref file from a source file for the various formats.
%.html.xref : %
	$(ROBODOC) $< $(@:.xref=) $(ROBOOPTS) INTERNAL GENXREF $@
%.tex.xref : %
	$(ROBODOC) $< $(@:.xref=) $(ROBOOPTS) INTERNAL GENXREF $@
%.rtf.xref : %
	$(ROBODOC) $< $(@:.xref=) $(ROBOOPTS) INTERNAL GENXREF $@

# Rule to create html documentation from a source file.
%.html : %
	$(ROBODOC) $< $@ HTML $(ROBOOPTS) XREF $(HTMLXREFSFILE)

# Rule to create latex documentation from a source file.
# We do not include source items, and generate laxtex documents
# than can be included in a master document.
%.tex : %
	$(ROBODOC) $< $@ LATEX $(ROBOOPTS) NOSOURCE SINGLEDOC XREF $(LATEXXREFSFILE)

# Rule to create ascii documentation from a source file.
%.txt : %
	$(ROBODOC) $< $@ ASCII 

# Rule to create rtf documentation from a source file.
%.rtf : %
	$(ROBODOC) $< $@ RTF $(ROBOOPTS) XREF $(RTFXREFSFILE)

# Use netscape to view the master index file for our project.
htmlview: html
	netscape $(PROJECT)_mi.html

# Use the latex programs to generate a .dvi from the master index file
# for our prokect. View this .dvi file with xdvi
texview:  tex
	latex $(PROJECT)_mi
	makeindex $(PROJECT)_mi
	latex $(PROJECT)_mi
	latex $(PROJECT)_mi
	xdvi  $(PROJECT)_mi.dvi

# Clean-up the mess we made
#
clean:
	rm -f $(HTMLXREFS) 
	rm -f $(HTMLDOCS) 
	rm -f $(LATEXXREFS)
	rm -f $(LATEXDOCS) 
	rm -f $(PROJECT)_mi.* *.aux
	rm -f $(RTFXREFS)
	rm -f $(RTFDOCS)
	rm -f $(ASCIIDOCS)
	rm -f $(HTMLXREFSFILE) 
	rm -f $(LATEXXREFSFILE) 
	rm -f $(RTFXREFSFILE)
