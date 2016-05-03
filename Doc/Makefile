HELPDOC=../../dev-tools/helpdoc

LATEX   = pdflatex
LATEX2HTML = latex2html

PDFS = user_guide.pdf tutorial.pdf developer_guide.pdf point_groups.pdf thermo.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)


doc:  all
all:  pdf html
pdf: $(PDFS)
html: user_guide developer_guide 

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) *~
	- rm -rf user_guide/ developer_guide/ 
	- rm -rf input_xx.xsl
	- rm -f INPUT_*.html INPUT_*.txt INPUT_*.xml
	- rm -rf ../../Doc/INPUT_PW.*


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
	@echo ""
	@echo "***"
	@echo "*** Developer's Guide created in developer_guide/developer_guide.html"
	@echo "***"
	@echo ""

