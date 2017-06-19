ifeq (,$(shell sh -c 'cygpath --version 2> /dev/null'))
  # Unix
  pwd := $$(pwd)
  translate = $1
else
  # Windows mit MSys2/Cygwin
  pwd := $$(cygpath -m "$$(pwd)")
  translate = $(shell echo '$1' | sed 's/:/;/g')
endif

all: build/eigenwerte.pdf

# hier Python-Skripte:
build/eigenwerte.pdf: plot_eigenwerte_vergleich.py | build
	python plot_eigenwerte_vergleich.py
#	TEXINPUTS="$(call translate,$(pwd):)" python plot.py

plot_eigenwerte.py: test.m | build
	octave test.m

test.m: Matlab_funktionen/H_F.m  Matlab_funktionen/H_m_n.m Matlab_funktionen/Hamilton_0.m

#hier weitere Abhängigkeiten für build/main.pdf deklarieren:

# build/main.pdf: FORCE | build
# 	  TEXINPUTS="$(call translate,build:)" \
# 	  BIBINPUTS=build: \
# 	  max_print_line=1048576 \
# 	latexmk \
# 	  --lualatex \
# 	  --output-directory=build \
# 	  --interaction=nonstopmode \
# 	  --halt-on-error \
# 	main.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
