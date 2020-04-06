
drawstars: drawstars.cpp file.h vmath.h image.h image_jpeg.h
	g++ drawstars.cpp -o drawstars -O3 -fdiagnostics-color=never -std=c++17 -ljpeg -lpng

.PHONY: examples
examples: drawstars
	cd examples/celestial-sphere && $(MAKE)
	cd examples/compare-to-real && $(MAKE)

.PHONY: catalogs
catalogs:
	cd gaia && $(MAKE)
	cd hip2 && $(MAKE)

.PHONY: all
all: catalogs examples
