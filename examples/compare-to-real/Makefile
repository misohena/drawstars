DRAWSTARS_DIR=../..
DRAWSTARS=$(DRAWSTARS_DIR)/drawstars --catalog-dir=$(DRAWSTARS_DIR)

.PHONY: all
all: 20190729.jpg 20190729-fovy53.jpg


20190729.jpg:
	(time $(DRAWSTARS) --proj=PERSPECTIVE --width=2400 --height=1600 --fovy=44.9 --view-z=0.18 --lat=36.699255 --lng=138.486832 --time="2019-07-29 20:58:56" --az=-5.58 --el=34.35 --roll=4.90 --flux-multiplier=M-7.5 --flux-inc-radius=M7.5 --post-gamma=2.2  --output=20190729.jpg --floats=20190729.floats)

20190729-fovy53.jpg:
	(time $(DRAWSTARS) --proj=PERSPECTIVE --width=2400 --height=1600 --fovy=53                 --lat=36.699255 --lng=138.486832 --time="2019-07-29 20:58:56" --az=-5.58 --el=34.35 --roll=4.90 --flux-multiplier=M-7.5 --flux-inc-radius=M7.5 --post-gamma=2.2  --output=20190729-fovy53.jpg --floats=20190729-fovy53.floats)

.PHONY: gamma
gamma:
	(time $(DRAWSTARS) --proj=LOAD --width=2400 --height=1600 --post-gamma=2.2 --output=20190729.jpg --floats=20190729.floats)

.PHONY: clean
clean:
	rm -f 20190729.jpg
	rm -f 20190729.floats
	rm -f 20190729-fovy53.jpg
	rm -f 20190729-fovy53.floats

.PHONY: redraw
redraw:
	$(MAKE) clean
	$(MAKE)
