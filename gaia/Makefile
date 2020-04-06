
ARCHIVE_DIR := $(shell cat ARCHIVE_DIR)
GAIA_CSV_LAST_FILE = $(ARCHIVE_DIR)/cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source/csv/MD5SUM.txt

gaia_ra_dec_g_bp_rp_teff.dat: $(GAIA_CSV_LAST_FILE) gaia_extract
#	./gaia_extract.py
	./gaia_extract

$(GAIA_CSV_LAST_FILE):
	./gaia_download.sh

gaia_extract: gaia_extract.cpp
	g++ gaia_extract.cpp -o gaia_extract -O3 -fdiagnostics-color=never -lz -std=c++17
