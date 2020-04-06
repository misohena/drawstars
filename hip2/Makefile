.PHONY: all
all: hip2_ra_dec_mag_bv.dat

# convert

hip2_ra_dec_mag_bv.dat: hip2.dat
	python hip2_extract.py

# download & gunzip

hip2.dat: hip2.dat.gz
	gunzip -c hip2.dat.gz > hip2.dat

hip2.dat.gz:
	wget http://cdsarc.u-strasbg.fr/ftp/I/311/hip2.dat.gz

# clean

.PHONY: clean
clean:
	rm hip2_ra_dec_mag_bv.dat

.PHONY: cleanall
cleanall: clean
	rm hip2.dat.gz hip2.dat
