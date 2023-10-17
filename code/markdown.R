library(rmarkdown)
rmarkdown::render('go3.Rmd', output_format=pdf_document(keep_tex = T))
