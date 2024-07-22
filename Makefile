# Makefile

all : go3-author.pdf go3.tex supplement_figs.pdf information-supplement.pdf

go3-author.pdf : go3.tex
	code/fix_authors.sh

# Rule to generate the metatranscriptome output files
go3.tex : go3.Rmd analysis/xt.Rda analysis/xg.Rda\
    analysis/xt.m.Rda analysis/res.Rda \
    analysis/x.all.Rda analysis/x.s.all.Rda\
    analysis/normalizations.Rda\
    analysis/normalizations.nl.Rda\
    analysis/yst.sens.Rda
	Rscript code/markdown.R

analysis/xt.Rda : code/meta-aldex.R
	Rscript code/meta-aldex.R

analysis/xg.Rda : code/meta-aldex.R
	Rscript code/meta-aldex.R

analysis/xt.m.Rda : code/meta-aldex.R
	Rscript code/meta-aldex.R

# Rules to generate the yeast output files

analysis/yst.sens.Rda : code/yst-sens.R
	Rscript code/yst-sens.R

analysis/res.Rda : code/yst-DESeq.R
	Rscript code/yst-DESeq.R

analysis/x.all.Rda : code/yst-aldex.R
	Rscript code/yst-aldex.R
	
analysis/x.s.all.Rda : code/yst-aldex.R
	Rscript code/yst-aldex.R

analysis/x.out.raw.Rda : code/scale-repro-aldex.R
	Rscript code/scale-repro-aldex.R

analysis/x.out.scale.Rda : code/scale-repro-aldex.R
	Rscript code/scale-repro-aldex.R
	
analysis/normalizations.nl.Rda : code/scale_norms.R
	Rscript code/scale_norms.R
	
analysis/normalizations.Rda : code/scale_norms.R
	Rscript code/scale_norms.R


# Rule for supplement
supplement_figs.pdf : supplement_figs.Rmd analysis/x.s.mu.all.Rda\
	analysis/meta.DES.res.Rda\
	analysis/data.out.Rda\
	analysis/out.mat.raw.Rda\
    analysis/out.mat.scale.Rda\
    analysis/sel.2.all.Rda\
    analysis/FP.df.Rda

	Rscript -e "rmarkdown::render('supplement_figs.Rmd')"

analysis/FP.df.Rda : code/yst-FP.R
	Rscript code/yst-FP.R
	
analysis/sel.2.all.Rda : code/selex.aldex.mu.R
	Rscript code/selex.aldex.mu.R
	
analysis/x.s.mu.all.Rda : code/yst.aldex.mu.r
	Rscript code/yst.aldex.mu.r
	
analysis/meta.DES.res.Rda : code/meta.DESeq.R
	Rscript code/meta.DESeq.R
	
analysis/data.out.Rda : code/scale_model.R
	Rscript code/scale_model.R
	
analysis/out.mat.raw.Rda : code/yst-reproducibility.R
	Rscript code/yst-reproducibility.R

analysis/out.mat.scale.Rda : code/yst-reproducibility.R
	Rscript code/yst-reproducibility.R
	
information-supplement.pdf : information-supplement.Rmd\

	Rscript -e "rmarkdown::render('information-supplement.Rmd')"

# Clean rule to remove the output files
clean:
	rm -f analysis/xt.all.Rda 
	rm -f analysis/xg.all.Rda 
	rm -f analysis/xt.m.all.Rda 
	rm -f analysis/res.Rda
	rm -f analysis/x.all.Rda
	rm -f analysis/x.s.all.Rda
	rm -f analysis/x.s.mu.all.Rda
	rm -f analysis/data.out.Rda