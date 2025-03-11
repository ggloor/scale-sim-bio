#!/bin/zsh

awk '{sub(/true \\and true \\and true/,"Gregory B. Gloor$^1$ \\and Michelle Pistner Nixon$^2$ \\and Justin D. Silverman$^{2,3,4,5}$")};1' go3.tex > jnk.tex

awk '{sub(/\\date{}/,"\\date{\\today}")};1' jnk.tex > jnk2.tex

awk '{sub(/\\maketitle/,"\\maketitle\n\n\$^1\$Department of Biochemistry, University of Western Ontario\n\$^2\$Department of Population Health Sciences, Geisinger, Danville, PA\n\$^3\$Department of Medicine, Pennsylvania State University\n\$^4\$Institute for Computational and Data Science, Pennsylvania State University\n\$^5\$Department of Statistics, Pennsylvania State University\n")};1' jnk2.tex > go3-author.tex

pdflatex go3-author.tex -o go3-author.pdf
bibtex go3-author
pdflatex go3-author.tex -o go3-author.pdf

rm jnk.tex
rm jnk2.tex
rm go3.pdf
rm go3.aux
rm go3.out
rm go3.toc
rm go3.log
rm go3-author.aux
rm go3-author.log
rm go3-author.blg

