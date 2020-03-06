Paper can be compiled into a pdf with the JOSS format using the following command (adapted from `rticles` R package):

```sh
~/.local/bin/pandoc +RTS -K512m -RTS paper.md --to latex --from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash --output paper.pdf --template ~/bin/R/x86_64-pc-linux-gnu-library/3.5/rticles/rmarkdown/templates/joss_article/resources/template.tex --highlight-style tango --pdf-engine pdflatex -V logo_path=~/bin/R/x86_64-pc-linux-gnu-library/3.5/rticles/rmarkdown/templates/joss_article/resources/JOSS-logo.png -V 'journal_name=Journal of Open Source Software' -V graphics=true --filter ~/.local/bin/pandoc-citeproc
```