# I'm adding this comment from github in my browser.
# Let's see if it shows up on the Mac, and on my phone.

# And now I'm making this comment on the Mac.
# Let's see if I can commit it and get it to show up on GitHub.

import astropy
import glob
from astropy.io import fits

# Grab all of the .fits and .fit files
all_files = glob.glob('/Users/jeff.c.taylor/Dropbox/ASTROINFORMATICs/RAWdata/RAW/*.fit*')

# Lists to store information
all_wavelengths = []
all_comments = []

for i in all_files:
        print("Processing " + i)
        hdulist = fits.open(i)
        wavelength = hdulist[0].header['waveleng']
        comment = hdulist[0].header.comments['waveleng']
        all_wavelengths.append(wavelength)
        all_comments.append(comment)

        print("\tWavelength: " + `wavelength`)
        print("\tComment: " + `comment`)

        hdulist.close()

#for i in all_wavelengths:
        #print("Wavelength value + comment: " + str(i))

#for i in all_comments:
        #print("Comment only: " + str(i))
