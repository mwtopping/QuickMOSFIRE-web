import astropy.io.fits
import numpy as np

# list of filenames to stack
images = ['12jan2019/co2_deep-H_262-419.fits','15mar2019/co2_deepH_mar2019_123-206.fits']
# list of weights for each images
weights = [5, 3]
# put a 1 if the images should be inverted
inverts = [0, 1]

#output file
output='combine.fits'
total = 0

# loop through each image
for invert, img, w in zip(inverts, images, weights):

    # read in and normalize each image
    data = astropy.io.fits.getdata(img)
    data = data / np.abs(np.median(data))

    # check if it is the first image read in
    if np.sum(total) == 0:
        total = w * data / np.sum(weights)

    else:
        # if you want it to be inverted then do that
        if invert:
            total -= w * data / np.sum(weights)
        else:
            total += w * data / np.sum(weights)

# write file
astropy.io.fits.writeto(output, total, overwrite=True) 
    


