# ascii reading
from astropy.io import ascii
import pandas as pd

# Catalog matching
from astropy.coordinates import SkyCoord
from astropy import units as u

# from astropy import units as u
# from datavis.ascii.match_catalogs import find_closest_in_cat

def match_catalogs(f_cat_ref, f_cat_sci,radius_threshold=2 * u.arcsec,
               ref_RA_KEY='X_WORLD',ref_DEC_KEY='Y_WORLD',
               sci_RA_KEY='X_WORLD',sci_DEC_KEY='Y_WORLD',reindex=True):
    '''Input ascii file for ref and sci catalogs, output pandas DataFrame of 
    matched ref catalog entries and matched sci catalog entries. All columns in inputs are preserved.
    '''
    cat_ref = ascii.read(f_cat_ref)
    cat_sci = ascii.read(f_cat_sci)
    df_ref  = pd.DataFrame(cat_ref.as_array())
    df_sci  = pd.DataFrame(cat_sci.as_array())
    coords_ref = SkyCoord(ra=df_ref[ref_RA_KEY], dec=df_ref[ref_DEC_KEY],unit='deg')
    coords_sci = SkyCoord(ra=df_sci[sci_RA_KEY], dec=df_sci[sci_DEC_KEY],unit='deg')
    
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
    
    # re-index to match two dfs
    if reindex:
        df_ref_matched_reindex=df_ref_matched.reset_index()
        df_sci_matched_reindex=df_sci_matched.reset_index()
        return df_ref_matched_reindex, df_sci_matched_reindex
    else:
        return df_ref_matched, df_sci_matched

def match_catalogs_df(df_ref, df_sci,radius_threshold=2 * u.arcsec,
               ref_RA_KEY='X_WORLD',ref_DEC_KEY='Y_WORLD',
               sci_RA_KEY='X_WORLD',sci_DEC_KEY='Y_WORLD',
               find_close=True, reindex=True):
    '''Input pandas data frame object for ref and sci catalogs, output pandas DataFrame of 
    matched ref catalog entries and matched sci catalog entries. All columns in inputs are preserved.
    '''
    coords_ref = SkyCoord(ra=df_ref[ref_RA_KEY], dec=df_ref[ref_DEC_KEY],unit='deg')
    coords_sci = SkyCoord(ra=df_sci[sci_RA_KEY], dec=df_sci[sci_DEC_KEY],unit='deg')
    
    # index = nearest in csci that matches cref
    # index = for each cref, get nearest csci
    # len(idx) = len(df_ref)
    idx, d2d, d3d = coords_ref.match_to_catalog_3d(coords_sci)
    
    # Sep constraint
    if find_close:
        sep_constraint = d2d <= radius_threshold
    else:
        sep_constraint = d2d > radius_threshold

    # Get entries in cat_ref with a match
    df_ref_matched = df_ref[sep_constraint]
    
    # Get matched entries in cat_sci
    df_sci_matched = df_sci.iloc[idx[sep_constraint]]

    # re-index to match two dfs
    if reindex:
        df_ref_matched_reindex=df_ref_matched.reset_index()
        df_sci_matched_reindex=df_sci_matched.reset_index()
        return df_ref_matched_reindex, df_sci_matched_reindex
    else:
        return df_ref_matched, df_sci_matched

def find_closest_in_cat(df_sci,RA,DEC,radius_threshold=2 * u.arcsec,
               sci_RA_KEY='X_WORLD',sci_DEC_KEY='Y_WORLD'):
    '''Input pandas DataFrame catalog, output line in DataFrame closest in input RA,DEC. 
    All columns in inputs are preserved.
    '''
    try:
        x=len(RA)
        coords_ref = SkyCoord(ra=RA, dec=DEC,unit='deg')
    except:
        coords_ref = SkyCoord(ra=[RA], dec=[DEC],unit='deg')
    coords_sci = SkyCoord(ra=df_sci[sci_RA_KEY], dec=df_sci[sci_DEC_KEY],unit='deg')
    
    # index = nearest in csci that matches input RA,DEC
    # index = for each input RA, DEC, get nearest csci
    # len(idx) = len(RA)
    idx, d2d, d3d = coords_ref.match_to_catalog_3d(coords_sci)
    
    # Sep constraint
    sep_constraint = d2d < radius_threshold
    
    # Get matched entries in cat_sci
    df_sci_matched = df_sci.iloc[idx[sep_constraint]]
    
    return df_sci_matched,d2d


def find_all_within(df_ref, RA,DEC,radius_threshold=2 * u.arcsec,
               ref_RA_KEY='X_WORLD',ref_DEC_KEY='Y_WORLD',
               find_close=True):
    '''Input pandas DataFrame object for ref, output all lines in pandas DataFrame within radius threshold of RA and DEC. 
    All columns in inputs are preserved.
    '''
    coords_ref = SkyCoord(ra=df_ref[ref_RA_KEY], dec=df_ref[ref_DEC_KEY],unit='deg')
    coords_sci = SkyCoord(ra=[RA], dec=[DEC],unit='deg')
    
    # index = nearest in csci that matches cref
    # index = for each cref, get nearest csci
    # len(idx) = len(df_ref)
    idx, d2d, d3d = coords_ref.match_to_catalog_3d(coords_sci)
    
    # Sep constraint
    if find_close:
        sep_constraint = d2d < radius_threshold
    else:
        sep_constraint = d2d > radius_threshold

    # Get entries in cat_ref with a match
    df_ref_matched = df_ref[sep_constraint]
    
    return df_ref_matched
