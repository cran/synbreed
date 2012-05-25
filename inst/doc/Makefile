nothing:
	echo "Doing nothing, to compile the vignette use make vignette" 

vignette: sweave
	pdflatex IntroSyn.tex
	bibtex IntroSyn
	pdflatex IntroSyn.tex
	pdflatex IntroSyn.tex
	rm IntroSyn.blg IntroSyn.log IntroSyn.out IntroSyn.aux IntroSyn.bbl 

sweave:
	echo "Sweave(\"IntroSyn.Rnw\")" | R --no-save --no-restore