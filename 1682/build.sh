python -c 'import gpkit; exec gpkit.mdmake("gas_hale.md")'
pandoc --template default.latex gas_hale.md.tex.md -o gas_hale.pdf


