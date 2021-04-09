from flask import Flask, render_template, request
import matplotlib.pyplot as plt
import numpy as np
from io import BytesIO
import base64
from data import *
from forms import RegistrationForm, FileForm
import astropy.io.fits
from pyquicklook import *

app = Flask(__name__)

app.config['SECRET_KEY'] = open('secret_key.txt', 'r').readline()

@app.route('/', methods=['GET', 'POST'])
def main():
    filesubmit = FileForm()
    if filesubmit.validate_on_submit():
        files = filesubmit.fileinput.data
        hdus = []
        for uploaded_file in files:

            uploaded_file.seek(0)

            hdus.append(astropy.io.fits.open(uploaded_file, memmap=True, mode='update'))


        newimg = quicklook(hdus, 'abab', 3.0, 0.0)
        hdu = hdus[0]
        fig = plt.Figure()
        ax = fig.subplots()
        ax.imshow(newimg, vmin=0.2*np.median(newimg), vmax=5*np.median(newimg))
        buf = BytesIO()
        fig.savefig(buf, format='png', dpi=300)
        data = base64.b64encode(buf.getbuffer()).decode("ascii")
        
        img = data

        return render_template('index.html', img=img, fileform=filesubmit)


    arr.append(len(arr))
    fig = plt.Figure()
    ax = fig.subplots()
    buf = BytesIO()
    fig.savefig(buf, format="png")
    data = base64.b64encode(buf.getbuffer()).decode("ascii")
    img = data
    return render_template('index.html', img=img, fileform=filesubmit)
