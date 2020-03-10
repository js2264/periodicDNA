To check how the final proof would look like in JOSS, the paper can be compiled into a pdf using the following command (adapted from `rticles` R package):

```sh
## JOSS format
~/.local/bin/pandoc +RTS -K512m -RTS paper.Rmd --to latex --from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash --output paper_JOSS.pdf --template ~/bin/R/x86_64-pc-linux-gnu-library/3.5/rticles/rmarkdown/templates/joss_article/resources/template.tex --highlight-style tango --pdf-engine pdflatex -V logo_path=~/bin/R/x86_64-pc-linux-gnu-library/3.5/rticles/rmarkdown/templates/joss_article/resources/JOSS-logo.png -V 'journal_name=Journal of Open Source Software' -V graphics=true --filter ~/.local/bin/pandoc-citeproc
## Regular PDF file
~/.local/bin/pandoc +RTS -K512m -RTS paper.Rmd --to latex --from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash --output paper.pdf --highlight-style tango --pdf-engine pdflatex --filter ~/.local/bin/pandoc-citeproc
```