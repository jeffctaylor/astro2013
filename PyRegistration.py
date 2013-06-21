#19th June 2013 - Image Registration: from IRAF  --->  PyRAF

#things to import regarding pyraf & iraf
#
import pyraf
from pyraf import iraf
#the following line is extremely important in order to get the tasks running properly without dependence on iraf's "login.cl" - it can be created anywhere and I just chose the current working directory
iraf.set(uparm='.')
#
from iraf import noao, images
from iraf import artdata, immatch, imcoords

#First we create an artificial fits image
# unlearn some iraf tasks
#
iraf.unlearn('mkpattern')
#create a fake image "apixelgrid.fits", to which we will register all fits images
#
artdata.mkpattern(input="apixelgrid.fits", output="apixelgrid.fits", pattern="constant", pixtype="double", ndim=2, ncols=100, nlines=100)
#note that in the exact above line, the "ncols" and "nlines" should be wisely chosen, depending on the input images - they provide the pixel-grid 
#for each input fits image, we will create the corresponding artificial one - therefore we can tune these values such that we cover, for instance, XXarcsecs of the target - so the best is that user provides us with such a value

#Then, we tag the desired WCS in this fake image:
# unlearn some iraf tasks
#
iraf.unlearn('ccsetwcs')
#tag the desired WCS in the fake image "apixel.fits"
#
iraf.ccsetwcs(images="apixelgrid.fits", database="", solution="", xref=50., yref=50., xmag=1.5, ymag=1.5, xrotati=0.,yrotati=0.,lngref=4.5133861111111, latref=64.849611111111, lngunit="hours", latunit="degrees", transpo="no", project="tan", coosyst="j2000", update="yes", pixsyst="logical", verbose="yes")
#note that the "xref" and "yref" are actually half the above "ncols", "nlines", respectively, so that we center each image
#note also that "xmag" and "ymag" is the pixel-scale, which in the current step ought to be the same as the native pixel-scale of the input image, for each input image - so we check the corresponding header value in each image
#note that "lngref" and "latref" can be grabbed by the fits header, it is actually the center of the target (e.g. ngc1569)
#note that we should make sure that the coordinate system is in coosyst="j2000" by checking the header info, otherwise we need to adjust that

# Then, register the fits file of interest to the WCS of the fake fits file
# unlearn some iraf tasks
#
iraf.unlearn('wregister')
#register the sciense fits image
#
iraf.wregister(input="scitest.fits", reference="apixelgrid.fits", output="scitestout.fits", fluxconserve="no")
#