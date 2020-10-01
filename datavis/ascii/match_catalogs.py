# ascii reading
from astropy.io import ascii
import pandas as pd

# Catalog matching
from astropy.coordinates import SkyCoord
from astropy import units as u

def match_catalogs(f_cat_ref, f_cat_sci,radius_threshold=2 * u.arcsec,
               ref_RA_KEY='X_WORLD',ref_DEC_KEY='Y_WORLD',
               sci_RA_KEY='X_WORLD',sci_DEC_KEY='Y_WORLD'):
    '''Input ascii file for ref and sci catalogs, output pandas DataFrame of 
    matched ref catalog entries and matched sci catalog entries. All columns in inputs are preserved.
    '''
    cat_ref = ascii.read(f_cat_ref)
    cat_sci = ascii.read(f_cat_sci)
    df_ref  = pd.DataFrame(cat_ref.as_array())
    df_sci  = pd.DataFrame(cat_sci.as_array())
    coords_ref = SkyCoord(ra=df_ref['X_WORLD'], dec=df_ref['Y_WORLD'],unit='deg')
    coords_sci = SkyCoord(ra=df_sci['X_WORLD'], dec=df_sci['Y_WORLD'],unit='deg')
    
    # index = nearest in csci that matches cref
    # index = for each cref, get nearest csci
    # len(idx) = len(df_ref)
    idx, d2d, d3d = coords_ref.match_to_catalog_3d(coords_sci)
    
    # Sep constraint
    sep_constraint = d2d < radius_threshold

    # Get entries in cat_ref with a match
    df_ref_matched = df_ref[sep_constraint]
    
    # Get matched entries in cat_sci
    df_sci_matched = df_sci.iloc[idx[sep_constraint]]
    
    return df_ref_matched, df_sci_matched
