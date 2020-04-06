DRAWSTARS_DIR=../..
DRAWSTARS=$(DRAWSTARS_DIR)/drawstars

equirectangular.jpg:
	(time $(DRAWSTARS) --catalog-dir=$(DRAWSTARS_DIR) --proj=EQUIRECTANGULAR --flux-gamma=1.2 --flux-multiplier=M-4.0 --flux-inc-radius=M3.0 --radius-default=1.20 --radius-max=2.5 --post-gamma=1.5 --jpeg-q=75 --output=equirectangular.jpg --floats=equirectangular.floats)

.PHONY: post
post:
	(time $(DRAWSTARS) --proj=LOAD --post-gamma=1.5 --jpeg-q=75 --output=equirectangular.jpg --floats=equirectangular.floats)

.PHONY: clean
clean:
	rm -f equirectangular.jpg
	rm -f equirectangular.floats

.PHONY: redraw
redraw:
	$(MAKE) clean
	$(MAKE)

.PHONY: serve
serve: equirectangular.jpg
	http-server -o /skyview.html
