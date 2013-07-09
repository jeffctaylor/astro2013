#19th June 2013 - Image Registration: from IRAF  --->  PyRAF
#2013-06-25 (JCT): What I want to do is add some command line parameters to this script.
#   For starters: the location of the image (so that it's not hardcoded) and the ncols
#   and nlines values, since Sophia has said it's best that the user provides these
#   values.

# A change from branch test1
# A change in branch test2

#things to import regarding pyraf & iraf

import pyraf
from pyraf import iraf
#the following line is extremely important in order to get the tasks running properly without dependence on iraf's "login.cl" - it can be created anywhere and I just chose the current working directory
iraf.set(uparm='.')

from iraf import noao, images
from iraf import artdata, immatch, imcoords

# getopt gives us some powerful command line processing tools
import sys
import getopt

import glob

from astropy.io import fits
from astropy.nddata import NDData
from astropy import units as u
from astropy import constants
import numpy as np

# Simple function to determine whether s is a number or not
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# Display usage information in case of a command line error
def print_usage():
    print("Usage: " + sys.argv[0] + " --ncols <ncols> --nlines <nlines> --image <image file>")

# This function will convert the input wavelength units into microns so that we only
# have to deal with a single unit.
def wavelength_to_microns(wavelength, unit):
    if (unit in u.micron.names or unit in u.um.names):
        return_value = float(wavelength)
    elif (unit in u.angstrom.names):
        return_value = u.angstrom.to(u.micron, float(wavelength))
    # This is a placeholder default value for now - it is not intended to be used
    # for real!
    else:
        return_value = 1

    return return_value

# A function to determine which instrument was used. This is done by checking certain
# keywords in the FITS header.
def get_instrument(header):
    instrument = ''
    # Check for INSTRUME keyword first
    if ('INSTRUME' in header):
        instrument = header['INSTRUME']
    # Then check for 'galex' in the INF0001 keyword
    elif ('INF0001' in header):
        if ('galex' in header['INF0001']):
            instrument = 'GALEX'
    # Finally, check to see if the 'ORIGIN' keyword has a value of '2MASS'
    elif ('ORIGIN' in header):
        if (header['ORIGIN'] == '2MASS'):
            instrument = '2MASS'
    else:
        print("could not determine instrument")
        sys.exit()

    return instrument

# A function to obtain the factor that is necessary to convert an image's flux units
# to Jy/pixel.
def get_conversion_factor(header, instrument):
    # Give a default value that can't possibly be valid; if this is still the value
    # after running through all of the possible cases, then an error has occurred.
    conversion_factor = 0

    if (instrument == 'IRAC'):
        print("IRAC; wavelength: " + `header['WAVELENG']` + "; CHNLNUM: " + `header['CHNLNUM']`)
        pixelscale = header['PXSCAL1']
        print("Pixel scale: " + `pixelscale`)
        # This is a hardcoded value from what Sophia gave me.
        # I would like to see if we could also obtain this from units.
        # MJy/sr to Jy/pixel
        conversion_factor = (2.3504 * 10**(-5)) * (pixelscale**2)

    elif (instrument == 'MIPS'):
        print("MIPS; wavelength: " + `header['WAVELENG']` + "; CHNLNUM: " + `header['CHNLNUM']`)
        pixelscale = header['PLTSCALE']
        print("Pixel scale: " + `pixelscale`)
        conversion_factor = (2.3504 * 10**(-5)) * (pixelscale**2)

    elif (instrument == 'GALEX'):
        print("GALEX; wavelength: " + `header['WAVELENG']`)
        print("Speed of light: " + `constants.c.to('um/s').value`)
        f_lambda_con = 0
        # I am using a < comparison here to account for the possibility that the given
        # wavelength is not EXACTLY 1520 AA or 2310 AA
        if (header['WAVELENG'] < 0.2): 
            f_lambda_con = 1.40 * 10**(-15)
        else:
            f_lambda_con = 2.06 * 10**(-16)
        conversion_factor = (10**23 * f_lambda_con * header['WAVELENG']**2) / (constants.c.to('um/s').value)

    elif (instrument == '2MASS'):
        print("2MASS; wavelength: " + `header['WAVELENG']` + "; FILTER: " + `header['FILTER']`)
        print("MAGZP: " + `header['MAGZP']`)
        fvega = 0
        if (header['FILTER'] == 'j'):
            fvega = 1594
        elif (header['FILTER'] == 'h'):
            fvega = 1024
        elif (header['FILTER'] == 'k'):
            fvega = 666.7
        conversion_factor = fvega * 10**(-0.4 * header['MAGZP'])

    elif (instrument == 'PACS'):
        print("PACS; wavelength: " + `header['WAVELENG']`)
        # Confirm that the data is already in Jy/pixel by checking the BUNIT header
        # keyword
        if ('BUNIT' in header):
            if (header['BUNIT'].lower() != 'jy/pixel'):
                print("Instrument is PACS, but Jy/pixel is not being used in BUNIT.")
        conversion_factor = 1;

    elif (instrument == 'SPIRE'):
        print("SPIRE; wavelength: " + `header['WAVELENG']` + "; pixelscale: " + `header['CDELT2']`)
        pixelscale = header['CDELT2']
        wavelength = header['WAVELENG']
        if (wavelength == 250):
            conversion_factor = (pixelscale**2) / 423
        elif (wavelength == 350):
            conversion_factor = (pixelscale**2) / 751
        elif (wavelength == 500):
            conversion_factor = (pixelscale**2) / 1587
    
    return conversion_factor

# A function to allow wavelength values to be obtained from a number of different
# header keywords
def get_wavelength(header):
    wavelength = 0
    wavelength_units = ''
    if ('WAVELENG' in header):
        wavelength = header['WAVELENG']
        wavelength_units = header.comments['WAVELENG']
    elif ('WAVELNTH' in header):
        wavelength = header['WAVELNTH']
        wavelength_units = 'micron'
    elif ('FILTER' in header):
        wavelength = header['FILTER']
        wavelength_units = 'micron'

    return wavelength, wavelength_units

ncols_input = ""
nlines_input = ""
image_input = ""

# Parse the command line options
try:
    opts, args = getopt.getopt(sys.argv[1:], "", ["ncols=", "nlines=", "image="])
except getopt.GetoptError:
    print("An error occurred. Check your parameters and try again.")
    sys.exit(2)
for opt, arg in opts:
    if opt in ("--ncols"):
        ncols_input = float(arg)
    if opt in ("--nlines"):
        nlines_input = float(arg)
    if opt in ("--image"):
        image_input= arg

# Now make sure that the values we have just grabbed from the command line are valid.
# ncols and nlines should both be numbers...
if (not is_number(ncols_input)):
    print_usage()
    print("ncols must be a number.")
    sys.exit()
if (not is_number(nlines_input)):
    print_usage()
    print("nlines must be a number.")
    sys.exit()

# And the image should be a valid file.
try:
   with open(image_input): pass
except IOError:
    print_usage()
    print("The specifiied image is not a file.")
    sys.exit()

# Read the input image data and header into an NDData object
hdulist = fits.open(image_input)
d1 = NDData(hdulist[0].data, meta=hdulist[0].header)
hdulist.close()

#print("Sample header value: " + d1.meta['OBJECT'])
lngref_input = d1.meta['CRVAL1']
latref_input = d1.meta['CRVAL2']
xmag_input = d1.meta['CDELT1']
ymag_input = d1.meta['CDELT2']

#print("lngref: " + `lngref_input`)
#print("latref: " + `latref_input`)
#print("xmag: " + `xmag_input` + "; converted: " + `u.deg.to(u.arcsec, xmag_input)`)
#print("ymag: " + `ymag_input` + "; converted: " + `u.deg.to(u.arcsec, ymag_input)`)
xmag_input = u.deg.to(u.arcsec, xmag_input)
ymag_input = u.deg.to(u.arcsec, ymag_input)
#print("xmag: " + `xmag_input`)

# Grab all of the .fits and .fit files
all_files = glob.glob('/Users/jeff.c.taylor/Dropbox/ASTROINFORMATICs/RAWdata/RAWb/*.fit*')

# Lists to store information
# We may need to another list to store filenames. This would enable us to create
# output filenames based on the input filenames; e.g. for image1.fits, we could output
# files like image1_registered.fits, image1_convolved.fits, image1_resampled.fits, etc.
# First check to make sure that this information is not already in the header somewhere.
# Actually, it can be found using hdulist.filename(), but it's not actually in the
# header. So another list seems like it would be the best way to retain it.
image_data = []
headers = []

for i in all_files:
    hdulist = fits.open(i)
    image = hdulist[0].data
    header = hdulist[0].header
    hdulist.close()
    #wavelength = header['WAVELENG']
    wavelength, wavelength_units = get_wavelength(header)
    #wavelength_units = header.comments['WAVELENG']
    wavelength_microns = wavelength_to_microns(wavelength, wavelength_units)
    #print("Wavelength " + `wavelength_microns` + "; Sample image value: " + `image[30][30]`)
    header['WAVELENG'] = (wavelength_microns, 'micron')
    image_data.append(image)
    headers.append(header)

# Sort the lists by their WAVELENG value
images_with_headers_unsorted = zip(image_data, headers)
images_with_headers = sorted(images_with_headers_unsorted, key=lambda header: header[1]['WAVELENG'])

#print("----------\nAfter sorting\n----------")
for i in range(0, len(images_with_headers)):
    #print("Data: " + `image_data`)
    #print("Data shape" + `image_data[i].shape`)
    wavelength = images_with_headers[i][1]['WAVELENG']
    wavelength_units = images_with_headers[i][1].comments['WAVELENG']
    #print("Wavelength: " + `wavelength` + ' ' + `wavelength_units` + "; Sample image value: " + `images_with_headers[i][0][30][30]`)
    #print("Wavelength units: " + `wavelength_units`)
    #wavelength_microns = wavelength_to_microns(wavelength, wavelength_units)
    #print("Wavelengths in microns: " + `wavelength_microns`)
    
    instrument = get_instrument(images_with_headers[i][1])
    print("Instrument: " + instrument)

    conversion_factor = get_conversion_factor(images_with_headers[i][1], instrument)
    print("Conversion factor: " + `conversion_factor`)

    # Now try to save all of the image data and headers as a new FITS image.
    # This was just a test - no need to actually do it right now, or yet.
    #hdu = fits.PrimaryHDU(image_data[i], headers[i])
    #hdu.writeto(`i` + '.fits')

sys.exit()

#First we create an artificial fits image
# unlearn some iraf tasks

iraf.unlearn('mkpattern')
#create a fake image "apixelgrid.fits", to which we will register all fits images

artdata.mkpattern(input="apixelgrid.fits", output="apixelgrid.fits", pattern="constant", pixtype="double", ndim=2, ncols=ncols_input, nlines=nlines_input)
#note that in the exact above line, the "ncols" and "nlines" should be wisely chosen, depending on the input images - they provide the pixel-grid 
#for each input fits image, we will create the corresponding artificial one - therefore we can tune these values such that we cover, for instance, XXarcsecs of the target - so the best is that user provides us with such a value

#Then, we tag the desired WCS in this fake image:
# unlearn some iraf tasks

iraf.unlearn('ccsetwcs')
#tag the desired WCS in the fake image "apixel.fits"

iraf.ccsetwcs(images="apixelgrid.fits", database="", solution="", xref=ncols_input/2.0, yref=nlines_input/2.0, xmag=xmag_input, ymag=ymag_input, xrotati=0.,yrotati=0.,lngref=lngref_input, latref=latref_input, lngunit="degrees", latunit="degrees", transpo="no", project="tan", coosyst="j2000", update="yes", pixsyst="logical", verbose="yes")
#note that the "xref" and "yref" are actually half the above "ncols", "nlines", respectively, so that we center each image
#note also that "xmag" and "ymag" is the pixel-scale, which in the current step ought to be the same as the native pixel-scale of the input image, for each input image - so we check the corresponding header value in each image
#note that "lngref" and "latref" can be grabbed by the fits header, it is actually the center of the target (e.g. ngc1569)
#note that we should make sure that the coordinate system is in coosyst="j2000" by checking the header info, otherwise we need to adjust that
# As of 2013-07-02, xmag, yman, lngref, and latref are all being obtained from the
# header of the input image

# Then, register the fits file of interest to the WCS of the fake fits file
# unlearn some iraf tasks

iraf.unlearn('wregister')
#register the sciense fits image

iraf.wregister(input=image_input, reference="apixelgrid.fits", output="scitestout.fits", fluxconserve="no")

