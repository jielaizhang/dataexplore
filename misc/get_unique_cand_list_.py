#from misc.get_unique_cand_list_ import get_unique_cand_list
#(unique_cand_RAs_1,
# unique_cand_DECs_1,
# cand_DETECTEDIN_Dict_1,
# d2d_list_1,
# unique_cand_EXPNAMEs,
# unique_cand_MAGSCIs,
# unique_CAND_MAGSUBs) = get_unique_cand_list(RA_all,DEC_all,EXPNAME_all,SCIMAG_all,SUBMAG_all,
#                                             same_cand_radii_threshold=1.0*u.arcsec) 



from astropy import units as u
from astropy.coordinates import SkyCoord

def get_unique_cand_list(RA_all,DEC_all,EXPNAME_all,MAGSCI_all,MAGSUB_all,same_cand_radii_threshold = 1*u.arcsec):

    unique_cand_RAs  = []
    unique_cand_DECs = []
    unique_cand_EXPNAMEs = []
    unique_cand_MAGSCIs = []
    unique_cand_MAGSUBs = []
    next_cand_number = 0
    d2d_list = []
    
    for RA,DEC,EXP_name,MAGSCI,MAGSUB in zip(RA_all,DEC_all,EXPNAME_all,MAGSCI_all,MAGSUB_all):

        #.................................................
        # Check if this is a unique candidate 
        if next_cand_number != 0:
            # Find closest candidate in existing cand list
            coords_thiscand      = SkyCoord(ra=[RA], dec=[DEC],unit='deg')
            coords_numberedcands = SkyCoord(ra=unique_cand_RAs, dec=unique_cand_DECs,unit='deg')
            # Find closest in numberedcands that matches thiscand
            idx, d2d, d3d = coords_thiscand.match_to_catalog_3d(coords_numberedcands)
            very_close = d2d<same_cand_radii_threshold

            # If not near any existing candidate in unique candidate, list, make new unique candidate
            if very_close[0]==False:
                unique_cand_RAs.append(RA)
                unique_cand_DECs.append(DEC)
                unique_cand_EXPNAMEs.append(EXP_name)
                unique_cand_MAGSCIs.append(str(MAGSCI))
                unique_cand_MAGSUBs.append(str(MAGSUB))
                next_cand_number += 1

            # If near existing canddiate, append detection exposure name to cand_DETECTEDIN_Dict
            elif very_close[0]==True:
                d2d_list.append(d2d[0])
                unique_cand_EXPNAMEs[idx[0]] = unique_cand_EXPNAMEs[idx[0]]+','+EXP_name
                unique_cand_MAGSCIs[idx[0]]  = unique_cand_MAGSCIs[idx[0]]+','+str(MAGSCI)
                unique_cand_MAGSUBs[idx[0]]  = unique_cand_MAGSUBs[idx[0]]+','+str(MAGSUB)

        #.................................................
        # Add candidate if first candidate
        if next_cand_number == 0:
            unique_cand_RAs.append(RA)
            unique_cand_DECs.append(DEC)
            unique_cand_EXPNAMEs.append(str(EXP_name))
            unique_cand_MAGSCIs.append(str(MAGSCI))
            unique_cand_MAGSUBs.append(str(MAGSUB))            
            next_cand_number += 1

    print('DONE')
    print('========')
    print(f'Total number of detected candidates: {len(RA_all)}')
    print(f'Number of unique candidates        : {len(unique_cand_RAs)}')
    
    return (unique_cand_RAs,unique_cand_DECs,
            d2d_list,
            unique_cand_EXPNAMEs,unique_cand_MAGSCIs,unique_cand_MAGSUBs)
