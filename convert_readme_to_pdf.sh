#!/usr/bin/env bash
# convert README in Markdown format to pdf
# requires Pandoc, https://pandoc.org
pandoc README.md --pdf-engine=xelatex -o README.pdf

