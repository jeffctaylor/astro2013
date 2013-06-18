# I'm adding this comment from github in my browser.
# Let's see if it shows up on the Mac, and on my phone.

import astropy
import glob
import pyfits

# Grab all of the .fits and .fit files
all_files = glob.glob('/Users/jeff.c.taylor/Dropbox/ASTROINFORMATICs/RAWdata/RAW/*.fit*')

# Lists to store information
all_wavelengths = []
all_naxis1 = []
all_naxis2 = []

for i in all_files:
        wavelength = pyfits.getval(i, 'WAVELENG')
        all_wavelengths.append(wavelength)
        naxis1 = pyfits.getval(i, 'NAXIS1')
        all_naxis1.append(naxis1)
        naxis2 = pyfits.getval(i, 'NAXIS2')
        all_naxis2.append(naxis2)

for i in all_wavelengths:
        print("Wavelength only: " + str(i))

for i in all_naxis1:
        print("NAXIS1 only: " + str(i))

for i in all_naxis2:
        print("NAXIS2 only: " + str(i))

print("Maximums: ")
print(str(max(all_naxis1)))
print(str(max(all_naxis2)))
