#!/usr/bin/env python
# For instructions on how to get results.csv and url.txt, see "How to X" slides from Jielai Zhang

#####################
### User Settings ###
#####################
f  = 'results.csv'# downloaded csv from NOAO archive after web search
f2 = 'url.txt' # contains single line with url where NOAO archive data is staged
desired_prod_type = 'image1' # choose what product type you want
f_out = 'wget.sh' # file where output wget commands will be written
#####################


import astropy.io.ascii as ascii
import pandas as pd
import numpy as np

# Read in csv downloaded from NOAO archive
d  = ascii.read(f)
filename_long = d['archive_filename']
filename      = np.array([x.split('/')[-1] for x in filename_long])
prod_type     = np.array(d['prod_type'])
exposure      = np.array(d['exposure'])
prop_id       = np.array(d['prop_id'])
release_date  = np.array(d['release_date'])
date_obs      = np.array(d['date_obs'])

# get filenames for product type image1
index = np.where(prod_type==desired_prod_type)
selectec_prod_type    = prod_type[index[0]]
selected_filename     = filename[index[0]]
selected_exposure     = exposure[index[0]]
selected_prop_id      = prop_id[index[0]]
selected_release_date = release_date[index[0]]
selected_date_obs     = date_obs[index[0]]
print('These images will be included in wget script output by this script.')
print('type      exp  prop_id   release    observe                         filename')
for tt,exp,ii,release,obs_date,nn in zip(selectec_prod_type,
                                        selected_exposure,
                                        selected_prop_id,
                                        selected_release_date,
                                        selected_date_obs,
                                        selected_filename):
    print(f'{tt} {exp:6} {ii} {release} {obs_date} {nn}')
print('')

# get url, strip new line char
myfile = open(f2, "r")
url = myfile.readline().strip()
print('This is the url used for wget from NOAO archive.')
print(url)
print('')
myfile.close()

# get wget commands, and write to file
writehere = open(f_out, "w")
print(f'These are the wget commands being written to an output file {f_out}')
for n in selected_filename:
    command = f'wget {url}/{n}'
    writehere.write(command+'\n')
    print(command)
writehere.close()
