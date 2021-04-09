from astropy.io import fits
import numpy as np




def quicklook(hdus, pattern, offset, prime):

    # these are the dither parameters converted into pixels
    pix_prime = int(prime/0.18)
    pix_offset = int(offset/0.18)

    # start with set number zero
    setno = 0
    ii = 0
    # the total number of dither patterns is the number of frames/4
    nsets = int(len(hdus)/4)

    # create a blank total image
    totalimg = np.zeros((2048-pix_prime-pix_offset,2048, nsets))

    # loop through each dither set
    for dither_set in range(nsets):

        # create blank images for the A and B dither positions
        A = np.zeros((2048-pix_prime, 2048))
        B = np.zeros((2048-pix_prime, 2048))

        # define the image numbers at the beginning and end of the set

        # loop through each image in the set
        for jj, i_img in enumerate(range(4)):

            # print out the file currently being read
            #  and then read in the data
            print('Processing {}'.format(hdus[ii][0]))
            imgdata = hdus[ii][0].data


            # define the sky frame based on which image we are reading in within
            #  the dither pattern
            print(ii, len(hdus))
            if jj != 0:
                sky1 = hdus[ii-1][0].data
            if jj != 3:
                sky2 = hdus[ii+1][0].data
            if jj == 0:
                sky1 = sky2
            if jj == 3:
                sky2 = sky1

            # subtract the sky image
            tmp_im = imgdata - (sky1+sky2)/2.

            # depending on where we are in the dither pattern, 
            #  set the sky-subtracted image to either A or B
            if even(i_img):
                A += tmp_im[:2048-pix_prime, :]
            else:
                B += tmp_im[:2048-pix_prime, :]

    
            ii += 1
        # save the total image from this one particluar dither pattern    
        totalimg[:,:,setno] = A[:2048-pix_prime-pix_offset, :]+B[pix_offset:2048-pix_prime, :]

        # move onto the next set
        setno += 1
    

    return np.median(totalimg, axis=2)
    # save the total image
#    hdu = fits.PrimaryHDU(np.median(totalimg, axis=2))
#    hdu.writeto('quick-m{}_{}-{}.fits'.format(utdate, ranges[0], ranges[1]))


# determine if the image number is even or not
def even(x):
    return x%2



if __name__=="__main__":
    quicklook([71, 142], 200311, 'abab', 2.7, 0.3)
