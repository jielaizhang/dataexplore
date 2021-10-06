import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.visualization import LinearStretch
from astropy.visualization import ZScaleInterval
from astropy.visualization import ImageNormalize



def make_stamps(RA,DEC,fitsfiles_2Darray,outpng='stamp.png',labels=False,size=50,debug=False,axisticks=False):
    '''fitsfiles_2Darray should be a 2D array that specifies the number of rows and columns in output stamps.
    E.g. [[temp_1.fits,temp_2.fits],[sci_1.fits,sci2.fits],[sub_1.fits,sub_2.fits]]
    '''
    n_cols = np.shape(fitsfiles_2Darray)[0]
    n_rows = np.shape(fitsfiles_2Darray)[1]
    if debug:
        print('n_cols, n_rows:',n_cols,n_rows)
        
    # Get WCS 
    f   = fitsfiles_2Darray[0][0]
    _,h = fits.getdata(f,header=True)
    wcs = WCS(h)
    (X, Y)       = wcs.wcs_world2pix(RA, DEC, 0)
    X,Y          = [int(x) for x in [X,Y]]
    wcs_cutout   = wcs[X-size:X+size,Y-size:Y+size] 
    
    # Set up subplot
    if axisticks == True:
        fig, axs = plt.subplots(n_rows, n_cols,figsize=(n_cols*5,n_rows*5),
                                subplot_kw={'projection': wcs_cutout, 'label':"overlays"})
    else:
        fig, axs = plt.subplots(n_rows, n_cols,figsize=(n_cols*5,n_rows*5))
    fig.subplots_adjust(left=None, bottom=None, right=1., top=None, wspace=None, hspace=None)
    
    for row in range(n_rows):
        for col in range(n_cols):
            if debug:
                print(f"fitsfiles_2Darray[{col}][{row}] #[col][row]")
                print(np.shape(axs))
                
            # Read in data
            f   = fitsfiles_2Darray[col][row]
            d,h = fits.getdata(f,header=True)
            
            # Do cutout
            wcs = WCS(h)
            (X, Y)       = wcs.wcs_world2pix(RA, DEC, 0)
            X,Y          = [int(x) for x in [X,Y]]
            image_cutout = d[X-size:X+size,Y-size:Y+size]
            wcs_cutout   = wcs[X-size:X+size,Y-size:Y+size] 
            
            # pixel value normalisation
            norm = ImageNormalize(image_cutout, interval=ZScaleInterval(),stretch=LinearStretch())
            
            # imshow
            if n_rows == 1 or n_cols ==1:
                axs[max(row,col)].imshow(image_cutout, norm = norm, cmap='gray') 
            else:
                axs[row, col].imshow(image_cutout, norm = norm, cmap='gray') 
            
            # labels
            if labels:
                label=labels[col][row]
            else:
                label=''
                
            # axis ticks   
            if axisticks == True:
                if n_rows == 1 or n_cols ==1:
                    axs[max(row,col)].coords[1].set_major_formatter('d.ddd')
                    axs[max(row,col)].coords[0].set_major_formatter('d.ddd')
                    axs[max(row,col)].coords[0].set_axislabel(f'\nRight Ascension \n {label} \n', fontsize=10)
                    axs[max(row,col)].coords[1].set_axislabel('Declination',fontsize=10)  
                else:
                    axs[row, col].coords[1].set_major_formatter('d.ddd')
                    axs[row, col].coords[0].set_major_formatter('d.ddd')
                    axs[row, col].coords[0].set_axislabel(f'\nRight Ascension \n {label} \n', fontsize=10)
                    axs[row, col].coords[1].set_axislabel('Declination',fontsize=10)  
            else:
                if n_rows == 1 or n_cols ==1:
                    axs[max(row,col)].set_axis_off()
                    axs[max(row,col)].set_title(label, fontsize=10)
                else:
                    axs[row, col].set_axis_off()
                    axs[row, col].set_title(label, fontsize=10)

    return outpng
