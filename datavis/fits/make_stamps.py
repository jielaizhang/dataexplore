import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.visualization import LinearStretch
from astropy.visualization import ZScaleInterval
from astropy.visualization import ImageNormalize


from datavis.fits.create_cutouts import create_cutout_centre


def make_stamps(RA,DEC,fitsfiles_2Darray,output='stamp.png',labels=False,size=50,debug=False,
                axisticks=False,mark_target=False):
    '''fitsfiles_2Darray should be a 2D array that specifies the number of rows and columns in output stamps.
    E.g. [[temp_1.fits,temp_2.fits],[sci_1.fits,sci2.fits],[sub_1.fits,sub_2.fits]]
    '''
            
    ##################
    # Set up subplot #
    ##################
    n_cols = np.shape(fitsfiles_2Darray)[0]
    n_rows = np.shape(fitsfiles_2Darray)[1]

    if axisticks == True:
        # Start figure with RA,DEC axis
        cutout,cutout_h = create_cutout_centre(fitsfiles_2Darray[0][0],RA,DEC,size,verbose=debug,debug=debug)
        cutout_wcs = WCS(cutout_h)
        fig, axs = plt.subplots(n_rows, n_cols,figsize=(n_cols*5,n_rows*5),
                                subplot_kw={'projection': cutout_wcs, 'label':"overlays"})
        if debug:
            print(f"cutout shape: {np.shape(cutout)}")
    else:
        # start figure with pixel axis
        fig, axs = plt.subplots(n_rows, n_cols,figsize=(n_cols*5,n_rows*5))
    fig.subplots_adjust(left=None, bottom=None, right=1., top=None, wspace=None, hspace=None)
    
    ##############################
    # plot one subplot at a time #
    ##############################
    for row in range(n_rows):
        for col in range(n_cols):
                
            # Read in data, get cutout
            cutout,cutout_h = create_cutout_centre(fitsfiles_2Darray[col][row],
                                            RA,DEC,size,verbose=debug,debug=debug)
            
            # imshow, with pixel value normalisation
            norm = ImageNormalize(cutout, interval=ZScaleInterval(),stretch=LinearStretch())
            if n_rows == 1 or n_cols ==1:
                axs[max(row,col)].imshow(cutout, norm = norm, cmap='gray') 
            else:
                axs[row, col].imshow(cutout, norm = norm, cmap='gray') 
                
            # labels
            if labels:
                label=labels[col][row]
            else:
                label=''
                
            # axis ticks   
            if axisticks == True:
                if n_rows == 1 or n_cols ==1:
                    axs[max(row,col)].set_title(label, fontsize=10)
                    axs[max(row,col)].coords[1].set_major_formatter('d.ddd')
                    axs[max(row,col)].coords[0].set_major_formatter('d.ddd')
                    axs[max(row,col)].coords[0].set_axislabel(f'Right Ascension', fontsize=10)
                    axs[max(row,col)].coords[1].set_axislabel('Declination',fontsize=10)  
                else:
                    axs[row, col].set_title(label, fontsize=10)
                    axs[row, col].coords[1].set_major_formatter('d.ddd')
                    axs[row, col].coords[0].set_major_formatter('d.ddd')
                    axs[row, col].coords[0].set_axislabel(f'Right Ascension', fontsize=10)
                    axs[row, col].coords[1].set_axislabel('Declination',fontsize=10)  
            else:
                if n_rows == 1 or n_cols ==1:
                    axs[max(row,col)].set_axis_off()
                    axs[max(row,col)].set_title(label, fontsize=10)
                else:
                    axs[row, col].set_axis_off()
                    axs[row, col].set_title(label, fontsize=10)
                    
            # mark target
            if mark_target==True:
                cutout_wcs = WCS(cutout_h)
                (X, Y)     = cutout_wcs.wcs_world2pix(RA, DEC, 0)
                if n_rows == 1 or n_cols ==1:
                    axs[max(row,col)].scatter(X,Y,s=200,marker='x',color='green',alpha=0.2)
                else:
                    axs[row, col].scatter(X,Y,s=200,marker='x',color='green',alpha=0.2)
                    
    plt.savefig(output)

    return output
