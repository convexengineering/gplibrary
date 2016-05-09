model.pdf: model.md
	python -c "from gpkit.tools import mdmake; exec mdmake('model.md')"
	pandoc -o model.pdf model.md.tex.md --template default.latex
