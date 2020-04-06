#!/usr/bin/env python

# Byte-by-byte Description of file: hip2.dat
# --------------------------------------------------------------------------------
#    Bytes Format Units    Label   Explanations
# --------------------------------------------------------------------------------
#    1-  6  I6    ---      HIP     Hipparcos identifier
#    8- 10  I3    ---      Sn      [0,159] Solution type new reduction (1)
#       12  I1    ---      So      [0,5] Solution type old reduction (2)
#       14  I1    ---      Nc      Number of components
#   16- 28 F13.10 rad      RArad   Right Ascension in ICRS, Ep=1991.25
#   30- 42 F13.10 rad      DErad   Declination in ICRS, Ep=1991.25
#   44- 50  F7.2  mas      Plx     Parallax
#   52- 59  F8.2  mas/yr   pmRA    Proper motion in Right Ascension
#   61- 68  F8.2  mas/yr   pmDE    Proper motion in Declination
#   70- 75  F6.2  mas    e_RArad   Formal error on RArad
#   77- 82  F6.2  mas    e_DErad   Formal error on DErad
#   84- 89  F6.2  mas    e_Plx     Formal error on Plx
#   91- 96  F6.2  mas/yr e_pmRA    Formal error on pmRA
#   98-103  F6.2  mas/yr e_pmDE    Formal error on pmDE
#  105-107  I3    ---      Ntr     Number of field transits used
#  109-113  F5.2  ---      F2      Goodness of fit
#  115-116  I2    %        F1      Percentage rejected data
#  118-123  F6.1  ---      var     Cosmic dispersion added (stochastic solution)
#  125-128  I4    ---      ic      Entry in one of the suppl.catalogues
#  130-136  F7.4  mag      Hpmag   Hipparcos magnitude
#  138-143  F6.4  mag    e_Hpmag   Error on mean Hpmag
#  145-149  F5.3  mag      sHp     Scatter of Hpmag
#      151  I1    ---      VA      [0,2] Reference to variability annex
#  153-158  F6.3  mag      B-V     Colour index
#  160-164  F5.3  mag    e_B-V     Formal error on colour index
#  166-171  F6.3  mag      V-I     V-I colour index
#  172-276 15F7.2 ---      UW      Upper-triangular weight matrix (G1)
# --------------------------------------------------------------------------------

#  32349   0 4 1  1.7678185359 -0.2916993748  379.21  -546.01 -1223.07   1.21   1.04   1.58   1.33   1.24  47  0.00  0    0.0    0 -1.0876 0.0024 0.040 0  0.009 0.007 -0.020   0.83  -0.04   0.96  -0.21   0.32   0.68   0.46  -0.12   0.13   0.90  -0.20   0.62  -0.14  -0.16   1.00

import gzip
from struct import *

with gzip.open('hip2.dat.gz', 'r') as infile:
    with open('hip2_ra_dec_mag_bv.dat', 'wb') as outfile:

        for line in infile:
            ra = float(line[15:28]) # rad
            dec = float(line[29:42]) # rad
            mag = float(line[129:136])
            bv = float(line[152:158])
            outfile.write(pack('f', ra)+pack('f', dec)+pack('f', mag)+pack('f', bv));

