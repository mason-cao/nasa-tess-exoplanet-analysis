# TOI-3505.01 manuscript build record

`manuscript_values.json` freezes the primary numbers, status boundary, author
line, novelty wording, and exact source-file inventory used to build the paper.
The editable prose is in `paper/TOI-3505.01_manuscript.md`, and the build script
is `src/build_toi3505_paper.py`.

The standard build writes a self-contained HTML manuscript. Passing
`--pdf-output default` also prints the final PDF with headless Chrome. The PDF
must be rendered to page images and visually inspected after any change to
prose, figures, tables, or print CSS.
