from astropy.io import fits
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

# the location of the mosfire data relative to the directory of the script
# also the location the result of this script is saved in
DATA_DIR = './'  # if they are in the same directory
DATA_DIR = './testdata/'
DATA_DIR = './mosfire_apr18-2021/'

def linefunc(xarr, a, b, h, mu, sigma):
    return a + b*xarr + h * np.exp(-((xarr-mu)/sigma/2)**2)

def quicklook(imgfilename, utdate, pattern, offset, prime):
    filenames = np.genfromtxt(imgfilename, dtype=str)

    # these are the dither parameters converted into pixels
    pix_prime = int(prime/0.18)
    pix_offset = int(offset/0.18)

    # load all of the images
    images = {}
    headers = {}
    frames = []
    backgrounds = []
    sigmas = []
    areas = []
    centers = []
    colors = []
    colorarr = ['#e41a1c','#377eb8','#4daf4a','#e78ac3']
    initial_star_position = 0
    for ii, filename in enumerate(filenames):
        print('Processing {}'.format(filename))
        processed_image, (bkgd, sigma, area, center) = preprocess(filename)
        key = filename[:-5]
        frameno = int(filename[8:-5])
        images[key] = processed_image
        headers[key] = fits.getheader(DATA_DIR + filename)
        if ii == 0:
            initial_star_position = center
        if ii%4==0:
            dither_offset = 0
        if ii%4==1:
            dither_offset = offset/0.18
        if ii%4==2:
            dither_offset = prime/0.18
        if ii%4==3:
            dither_offset = (prime+offset)/0.18
        if area > 0:
            frames.append(frameno)
            backgrounds.append(bkgd)
            sigmas.append(sigma)
            areas.append(area)
            centers.append(center-dither_offset)
            colors.append(colorarr[ii%4])
        else:
            print("Cannot fit star for frame number {}".format(filename))


    # plot all of the observing summaries
    summaryfig, summaryaxs = plt.subplots(2, 2, figsize=(12, 10))
    summaryaxs[0,0].scatter(frames, (np.array(centers) - centers[0])*.18, c=colors, s=96, edgecolors='black', linewidths=0.5)
#    summaryaxs[0,0]
    summaryaxs[0,0].set_xlabel('Frame Number')
    summaryaxs[0,0].set_ylabel('Star drift [arcsec]')
    summaryaxs[0,0].text(0.05, 0.9, "A", transform=summaryaxs[0,0].transAxes, horizontalalignment='left', color=colorarr[0], fontsize=16, weight='bold')
    summaryaxs[0,0].text(0.05, 0.84, "B", transform=summaryaxs[0,0].transAxes, horizontalalignment='left', color=colorarr[1], fontsize=16, weight='bold')
    summaryaxs[0,0].text(0.05, 0.78, "A'", transform=summaryaxs[0,0].transAxes, horizontalalignment='left', color=colorarr[2], fontsize=16, weight='bold')
    summaryaxs[0,0].text(0.05, 0.72, "B'", transform=summaryaxs[0,0].transAxes, horizontalalignment='left', color=colorarr[3], fontsize=16, weight='bold')

    summaryaxs[0, 1].scatter(frames, sigmas, c=colors, s=96, edgecolors='black', linewidths=0.5)
    summaryaxs[0, 1].set_xlabel('Frame Number')
    summaryaxs[0, 1].set_ylabel('Star FWHM [arcsec]')

    summaryaxs[1,0].scatter(frames, areas, c=colors, s=96, edgecolors='black', linewidths=0.5)
    summaryaxs[1, 0].set_xlabel('Frame Number')
    summaryaxs[1, 0].set_ylabel('Star Flux [arbitrary]')

    summaryaxs[1,1].scatter(frames, backgrounds, c=colors, s=96, edgecolors='black', linewidths=0.5)
    summaryaxs[1, 1].set_xlabel('Frame Number')
    summaryaxs[1, 1].set_ylabel('Background level')

    # create a blank total image
    s = np.shape(processed_image)[1]

    totalimg = np.zeros((2048 - pix_prime - pix_offset, s, int(len(filenames)/4)))

    print("Creating Stack")
    # loop through each dither set

    for dither_set in range(int(len(filenames)/4)):

        A = np.zeros((2048 - pix_prime, s))
        B = np.zeros((2048 - pix_prime, s))

        # loop through each image in the set
        for ii in range(4):
            key = filenames[dither_set*4+ii][:-5]
            # print out the file currently being read
            #  and then read in the data
            imgdata = images[key]

            # define the sky frame based on which image we are reading in within
            #  the dither pattern
            if ii != 0:
                sky1 = images[filenames[dither_set*4+ii-1][:-5]]
            if ii != 3:
                sky2 = images[filenames[dither_set*4+ii+1][:-5]]
            if ii == 0:
                sky1 = sky2
            if ii == 3:
                sky2 = sky1

            # subtract the sky image
            tmp_im = imgdata - (sky1+sky2)/2.

            # depending on where we are in the dither pattern,
            #  set the sky-subtracted image to either A or B
            if even(dither_set*4+ii):
                A += tmp_im[:2048-pix_prime, :]
            else:
                B += tmp_im[:2048-pix_prime, :]


        # save the total image from this one particluar dither pattern
        totalimg[:,:,dither_set] = A[:2048-pix_prime-pix_offset, :]+B[pix_offset:2048-pix_prime, :]

#    print(imgheader.keys)
    summedimage = np.median(totalimg, axis=2)
    hdu = fits.PrimaryHDU(summedimage)
    hdu.writeto(DATA_DIR+'quick-m{}_{:.1f}hrs_H_shifted.fits'.format(utdate, len(filenames)*120/3600), overwrite=True)
    plt.savefig(DATA_DIR+'quick-m{}_{:.1f}hrs_H_summary.png'.format(utdate, len(filenames)*120/3600), dpi=200)
    print("Completed")

# determine if the image number is even or not
def even(x):
    return x%2


def preprocess(filename):
    imgdata = fits.getdata(DATA_DIR + filename)
    imgheader = fits.getheader(DATA_DIR + filename)

    # attempt to measure the star parameters
    try:
        # sum the image
        slice = np.median(imgdata, axis=1)
        slicesize=20
        slicemin = slice.argmax()-slicesize
        slicemax = slicemin + 2*slicesize

        # fit the line
        popt, pcov = curve_fit(linefunc, range(2*slicesize), slice[slicemin:slicemax],
                               p0=[1, 0, np.max(slice), slicesize, 3])


        height = popt[2]
        sigma = popt[4]
        center = popt[3]
        background = popt[0]
        area = height * sigma * np.sqrt(2*3.14)
    except:
        sigma = -1
        center = -1
        background = -1
        area = -1

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

        subimg = imgdata[int(newbarypos[ii]):int(newbarypos[ii+1]),:]

        sub_height = np.shape(subimg)[0]
        sub_width = np.shape(subimg)[1]

        dx = int(maxxpix)-int(newbarxpos[ii+1])


        subimg = np.pad(subimg[::-1, :], ((2048-int(newbarypos[ii])-sub_height, int(newbarypos[ii])),
                                          (dx,new_width-sub_width-dx)))
        newimg += subimg

    return newimg[::-1,:], (background, 2.355*.18*sigma, area, slicemin+center)

if __name__=="__main__":
    quicklook(DATA_DIR+'imlist.txt', 210409, 'abab', 2.7, 0.3)
