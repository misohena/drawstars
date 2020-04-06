#!/usr/bin/env python
import os
import glob
import gzip
import struct

with open('ARCHIVE_DIR', 'r') as fd:
    archive_dir = fd.readline().rstrip('\n')

print('archive_dir=' + archive_dir);

dir = archive_dir + "/cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source/csv/";
files = glob.glob(os.path.join(dir, "*.csv.gz"));

# number of stars: 1,692,919,135 < 2^31 (int)
num_output_star = 0
num_invalid_star = 0

with open('gaia_ra_dec_g_bp_rp_teff.dat', 'wb') as outfile:
    for file_index, filename in enumerate(files):
        print(str(file_index) + "/" + str(len(files)) + ": " + os.path.basename(filename))

        with gzip.open(filename, "rt") as csv_file:
            column_names = csv_file.readline().rstrip("\n\r").split(",")
            index_ra = column_names.index("ra");
            index_dec = column_names.index("dec");
            index_g_mag = column_names.index("phot_g_mean_mag");
            index_bp_mag = column_names.index("phot_bp_mean_mag");
            index_rp_mag = column_names.index("phot_rp_mean_mag");
            index_teff = column_names.index("teff_val");
            #index_bp_rp = column_names.index("bp_rp");
            #index_bp_g = column_names.index("bp_g");

            for line in csv_file:
                columns = line.rstrip("\n\r").split(",")
                try:
                    ra = float(columns[index_ra])
                    dec = float(columns[index_dec])
                    g_mag = float(columns[index_g_mag])
                    bp_mag = float(columns[index_bp_mag])
                    rp_mag = float(columns[index_rp_mag])
                    teff = float(columns[index_teff])

                    outfile.write(
                        struct.pack('f', ra) +
                        struct.pack('f', dec) +
                        struct.pack('f', g_mag) +
                        struct.pack('f', bp_mag) +
                        struct.pack('f', rp_mag) +
                        struct.pack('f', teff));
                    num_output_star = num_output_star + 1
                except:
                    num_invalid_star = num_invalid_star + 1

print("output: " + str(num_output_star))
print("invalid: " + str(num_invalid_star))
