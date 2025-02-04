LATEX   = pdflatex
LATEX2HTML = latex2html

PDFS = user_guide.pdf tutorial.pdf developer_guide.pdf point_groups.pdf thermo.pdf units.pdf equilibrium.pdf symmetry.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)


doc:  pdf
all:  pdf html
pdf: $(PDFS)
html: user_guide developer_guide tutorial

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) *~
	- rm -rf user_guide/ developer_guide/ thermo/ point_groups/ tutorial/
	- rm -rf figure_hex.ps
	- rm -rf input_xx.xsl


user_guide: user_guide.pdf
	rm -rf user_guide/
	latex2html \
                -t "User's Guide for the Thermo_pw package" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                user_guide.tex
	cd user_guide; \
	for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	cp user_guide.css user_guide
	@echo ""
	@echo "***"
	@echo "*** User's Guide created in user_guide/user_guide.html"
	@echo "***"
	@echo ""

developer_guide: developer_guide.pdf
	rm -rf developer_guide/
	latex2html \
                -t "Developer's Guide for the Thermo_pw package" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                developer_guide.tex
	cd developer_guide; \
	for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	cp user_guide.css developer_guide/developer_guide.css
	@echo ""
	@echo "***"
	@echo "*** Developer's Guide created in developer_guide/developer_guide.html"
	@echo "***"
	@echo ""

thermo: thermo.pdf
	rm -rf thermo/
	latex2html \
                -t "Thermo guide for the Thermo_pw package" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                thermo.tex
	cd thermo; \
	for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	cp user_guide.css thermo/thermo.css
	@echo ""
	@echo "***"
	@echo "*** Thermo Guide created in thermo/thermo.html"
	@echo "***"
	@echo ""

point_groups: point_groups.pdf
	rm -rf point_groups/
	latex2html \
                -t "Crystallographic point groups guide for the Thermo_pw package" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                point_groups.tex
	cd point_groups; \
	for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	cp user_guide.css point_groups/point_groups.css
	@echo ""
	@echo "***"
	@echo "*** Point Groups Guide created in point_groups/point_groups.html"
	@echo "***"
	@echo ""

tutorial: tutorial.pdf
	rm -rf tutorial/
	latex2html \
                -t "Tutorial for the Thermo_pw package" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                tutorial.tex
	cd tutorial; \
	for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	cp user_guide.css tutorial/tutorial.css
	@echo ""
	@echo "***"
	@echo "*** Tutorial created in tutorial/tutorial.html"
	@echo "***"
	@echo ""

