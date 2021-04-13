from astropy.io import fits
from scipy import ndimage
import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = './testdata/'


def quicklook(ranges, utdate, pattern, offset, prime):

    # these are the dither parameters converted into pixels
    pix_prime = int(prime/0.18)
    pix_offset = int(offset/0.18)

    # load all of the images
    images = {}
    headers = {}
    for imgno in np.arange(ranges[0], ranges[1]+1, 1):
        images['{}_{:04d}'.format(utdate, imgno)] = preprocess(utdate, imgno)
        headers['{}_{:04d}'.format(utdate, imgno)] = fits.getheader(DATA_DIR + 'm{}_{:04d}.fits'.format(utdate, imgno))

    # start with set number zero
    setno = 0
    # the total number of dither patterns is the number of frames/4
    nsets = int((ranges[1] - ranges[0] + 1)/4)

    # create a blank total image
    imgdata = preprocess(utdate, ranges[0])
    s = np.shape(imgdata)[1]
    print("This is happening here")

    totalimg = np.zeros((2048 - pix_prime - pix_offset, s, nsets))

    # loop through each dither set
    for dither_set in range(nsets):

        A = np.zeros((2048 - pix_prime, s))
        B = np.zeros((2048 - pix_prime, s))

        # define the image numbers at the beginning and end of the set
        set_start = ranges[0]+setno*4
        set_end = (ranges[0]+3)+setno*4

        # loop through each image in the set
        for ii, i_img in enumerate(np.arange(set_start, set_end+1, 1)):
            # print out the file currently being read
            #  and then read in the data
            print('Processing m{}_{:04d}.fits'.format(utdate, i_img))
#            imgdata = fits.getdata(DATA_DIR+'m{}_{:04d}.fits'.format(utdate, i_img))
#            imgdata = preprocess(utdate, i_img)
            imgdata = images['{}_{:04d}'.format(utdate, i_img)]
            #imgheader = fits.getheader(DATA_DIR+'m{}_{:04d}.fits'.format(utdate, i_img))

            # define the sky frame based on which image we are reading in within
            #  the dither pattern
            if i_img != set_start:
                #sky1 = fits.getdata(DATA_DIR+'m{}_{:04d}.fits'.format(utdate, i_img-1))
#                sky1 = preprocess(utdate, i_img-1)
                sky1 = images['{}_{:04d}'.format(utdate, i_img-1)]
            if i_img != set_end:
                #sky2 = fits.getdata(DATA_DIR+'m{}_{:04d}.fits'.format(utdate, i_img+1))
#                sky2 = preprocess(utdate, i_img+1)
                sky2 = images['{}_{:04d}'.format(utdate, i_img + 1)]
            if i_img == set_start:
                sky1 = sky2
            if i_img == set_end:
                sky2 = sky1

            # subtract the sky image
            tmp_im = imgdata - (sky1+sky2)/2.
            print("showing sky image")

            # depending on where we are in the dither pattern,
            #  set the sky-subtracted image to either A or B
            if even(i_img):
                A += tmp_im[:2048-pix_prime, :]
            else:
                B += tmp_im[:2048-pix_prime, :]

    
        # save the total image from this one particluar dither pattern    
        totalimg[:,:,setno] = A[:2048-pix_prime-pix_offset, :]+B[pix_offset:2048-pix_prime, :]

        # move onto the next set
        setno += 1

#    print(imgheader.keys)
    summedimage = np.median(totalimg, axis=2)

    fig, ax = plt.subplots(figsize=(12, 12))
    plt.imshow(summedimage, origin='lower', vmin=-30, vmax=30)
    print(summedimage)
    plt.show()
#    Nbars = 92
#    barnums = range(92)
#    barxpos = []
#    barypos = []
#    for ii, b in enumerate(barnums):
#        barypos.append((45-int(ii/2))*7.982/.18)
#        print((45-int(ii/2))*7.982/.18)
#        barxpos.append((367.2-imgheader['B{:02d}POS'.format(b+1)])/0.18)
#
#    # find the pixel of the left and right extremes of the bars
#    minxpix = np.min(barxpos)
#    maxxpix = np.max(barxpos)
#
#    # create a new image that can hold the shifted data
#    newimg = np.zeros((2048, 2048+np.int(maxxpix)-np.int(minxpix)))
#    new_width = np.shape(newimg)[1]
#
#    fig, ax = plt.subplots(figsize=(12, 12))
#    plt.imshow(summedimage, vmin=-30, vmax=30)
#
#    lastx = -10
#    # loop through all of the bars
#    for ii, (x, y) in enumerate(zip(barxpos, barypos)):
#        # only look at the right bars
#        if ii%2 and (not ii==0):
#            continue
#        print(x, lastx)
#        if np.abs(x-lastx) > 1/.18:
#            # this is a new bar
#            ax.axhline(y=int(y), color='black')
#        else:
#            ax.axhline(y=int(y), color='black', linestyle=':', linewidth=.8)
#        # get the sub image of this bar
#        # they are roughly
#        subimg = summedimage[int(y-pix_offset+pix_prime):int(y-pix_offset+pix_prime)+44+ii%2,:]
#        #subimg = summedimage[int(y-10):int(y-10)+44+ii%2,:]
#
#        sub_height = np.shape(subimg)[0]
#        sub_width = np.shape(subimg)[1]
#
#        dx = int(maxxpix)-int(x)
#
#
#        subimg = np.pad(subimg[::-1, :], ((2048-int(y)-sub_height, int(y)), (dx,new_width-sub_width-dx) ))
#        newimg += subimg
#        lastx = x
#
#    # save the total image
#    hdu = fits.PrimaryHDU(np.median(totalimg, axis=2))
#    hdu.writeto(DATA_DIR+'quick-m{}_{}-{}.fits'.format(utdate, ranges[0], ranges[1]), overwrite=True)
    hdu = fits.PrimaryHDU(summedimage)
    hdu.writeto(DATA_DIR+'quick-m{}_{}-{}_shifted.fits'.format(utdate, ranges[0], ranges[1]), overwrite=True)
#
#    plt.show()

# determine if the image number is even or not
def even(x):
    return x%2


def preprocess(utdate, i_img):
    imgdata = fits.getdata(DATA_DIR + 'm{}_{:04d}.fits'.format(utdate, i_img))
    imgheader = fits.getheader(DATA_DIR + 'm{}_{:04d}.fits'.format(utdate, i_img))

    Nbars = 92
    barnums = range(92)
    barxpos = []
    barypos = []
    for ii, b in enumerate(barnums):
        barypos.append((46-int((ii)/2))*7.982/.18)
        barxpos.append((367.2-imgheader['B{:02d}POS'.format(b+1)])/0.18)

    # refine the bar positions
    lastx = -10
    newbarxpos = []
    newbarypos = []
    for ii, (x, y) in enumerate(zip(barxpos, barypos)):
        if ii%2 and (not ii==0):
            continue
        if np.abs(x-lastx) > 1/.18:
            newbarxpos.append(x)
            newbarypos.append(y)
            # this is a new bar
        lastx=x

    # add in the last one
    newbarxpos.append(x)
    newbarypos.append(0)

    print(np.array(newbarxpos).astype(int))

#    for y in newbarypos:
#        print(y)
#        ax.axhline(y=y, color='white', linewidth=2)
    # find the pixel of the left and right extremes of the bars
    minxpix = np.min(barxpos)
    maxxpix = np.max(barxpos)

    # create a new image that can hold the shifted data
    newimg = np.zeros((2048, 2048+np.int(maxxpix)-np.int(minxpix)))
    new_width = np.shape(newimg)[1]

    # reverse the orders
    newbarxpos = newbarxpos[::-1]
    newbarypos = newbarypos[::-1]

    # loop through all of the bars
    for ii in range(len(newbarxpos)-1):

#        print(int(newbarypos[ii]),int(newbarypos[ii+1]))

        subimg = imgdata[int(newbarypos[ii]):int(newbarypos[ii+1]),:]

        sub_height = np.shape(subimg)[0]
        sub_width = np.shape(subimg)[1]

        dx = int(maxxpix)-int(newbarxpos[ii+1])


#        if ii == 0:
#            subimg = np.pad(subimg[::-1, :], ((2048-sub_height, 0),
#                                          (dx,new_width-sub_width-dx)))
#        else:
        subimg = np.pad(subimg[::-1, :], ((2048-int(newbarypos[ii])-sub_height, int(newbarypos[ii])),
                                          (dx,new_width-sub_width-dx)))
        newimg += subimg

#    fig, ax = plt.subplots(figsize=(12, 12))
#    plt.imshow(newimg[::-1,:], origin='lower', vmin=20, vmax=100)
#    plt.show()

    return newimg[::-1,:]

if __name__=="__main__":
    quicklook([236, 259], 210403, 'abab', 2.7, 0.3)
#    preprocess(210403, 256)
