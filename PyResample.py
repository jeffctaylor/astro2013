#11th July 2013 - Resampling using pyraf

#things to import regarding pyraf & iraf
#
import pyraf
from pyraf import iraf
#the following line is extremely important in order to get the tasks running properly without dependence on iraf's "login.cl" - it can be created anywhere and I just chose the current working directory
iraf.set(uparm='.')
#
from iraf import noao, images
from iraf import artdata, immatch, imcoords

# First we create an artificial fits image, 
# The difference with the registration step is that the artificial image is now created only once, and it is common for all the input_images_convolved (or imput_images_gaussian_convolved)
# unlearn some iraf tasks
#
iraf.unlearn('mkpattern')
#
#create a fake image "grid_final_resample.fits", to which we will register all fits images
#
artdata.mkpattern(input="grid_final_resample.fits", output="grid_final_resample.fits", pattern="constant", pixtype="double", ndim=2, ncols=parameter1, nlines=parameter2)
# parameter1 & parameter2 depend on the "fwhm" of the convolution step, and following the Nyquist sampling rate. 
# parameter1 = parameter2 = physical_zize / (fwhm / 3.3) 

#Then, we tag the desired WCS in this fake image:
# unlearn some iraf tasks
#
iraf.unlearn('ccsetwcs')
#tag the desired WCS in the fake image "apixel.fits"
#
iraf.ccsetwcs(images="grid_final_resample.fits", database="", solution="", xref=ncols/2, yref=nlines/2, xmag=fwhm/3.3, ymag=fwhm/3.3, xrotati=0.,yrotati=0.,lngref=same_as_we_put_in_the_registration, latref=same_as_we_put_in_the_registration, lngunit="hours", latunit="degrees", transpo="no", project="tan", coosyst="j2000", update="yes", pixsyst="logical", verbose="yes")

# Then, register the fits file of interest to the WCS of the fake fits file
# unlearn some iraf tasks
#
iraf.unlearn('wregister')
#register the sciense fits image
#
iraf.wregister(input="input_image_convolved.fits", reference="grid_final_resample.fits", output="input_image_final.fits", fluxconserve="yes")
#
