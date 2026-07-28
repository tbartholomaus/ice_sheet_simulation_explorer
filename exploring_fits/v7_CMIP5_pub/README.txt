# ISMIP6-Greenland scalar model output
# Heiko Goelzer (h.goelzer@uu.nl) 2020
# Version v0 - November 28 2019
# Version v1 - November 30 2019
# Version v4 - December 22 2019 - Used for paper submissions
# Version v5 - Feb 13 2020 - including ctrl_exp and historical
# Version v6 - April 18 2020 
# Version v7 - May 15 2020 - All models, including Zwally basins
# Version v7_pub - July 15 2020 - same as v7 but renamed models to match publication

Data usage notice:
If you use any of these results, please acknowledge the work of the people involved in the process producing this data set. Acknowledgements should have language similar to the below (if you only use CMIP5 forcing, remove CMIP6 and vice versa).

“We thank the Climate and Cryosphere (CliC) effort, which provided support for ISMIP6 through sponsoring of workshops, hosting the ISMIP6 website and wiki, and promoted ISMIP6. We acknowledge the World Climate Research Programme, which, through it's Working Group on Coupled Modelling, coordinated and promoted CMIP5 and CMIP6. We thank the climate modeling groups for producing and making available their model output, the Earth System Grid Federation (ESGF) for archiving the CMIP data and providing access, the University at Buffalo for ISMIP6 data distribution and upload, and the multiple funding agencies who support CMIP5 and CMIP6 and ESGF. We thank the ISMIP6 steering committee, the ISMIP6 model selection group and ISMIP6 dataset preparation group for their continuous engagement in defining ISMIP6."

You should also refer to and cite the following papers:

Heiko Goelzer, Sophie Nowicki, Anthony Payne, Eric Larour, Helene Seroussi, William H. Lipscomb, Jonathan Gregory, Ayako Abe-Ouchi, Andy Shepherd, Erika Simon, Cecile Agosta, Patrick Alexander, Andy Aschwanden, Alice Barthel, Reinhard Calov, Christopher Chambers, Youngmin Choi, Joshua Cuzzone, Christophe Dumas, Tamsin Edwards, Denis Felikson, Xavier Fettweis, Nicholas R. Golledge, Ralf Greve, Angelika Humbert, Philippe Huybrechts, Sebastien Le clec'h, Victoria Lee, Gunter Leguy, Chris Little, Daniel P. Lowry, Mathieu Morlighem, Isabel Nias, Aurelien Quiquet, Martin Rückamp, Nicole-Jeanne Schlegel, Donald Slater, Robin Smith, Fiamma Straneo, Lev Tarasov, Roderik van de Wal, and Michiel van den Broeke: The future sea-level contribution of the Greenland ice sheet: a multi-model ensemble study of ISMIP6 , The Cryosphere, 2020. doi:10.5194/tc-2019-319

Sophie Nowicki, Antony Payne, Heiko Goelzer, Helene Seroussi, William Lipscomb, Ayako Abe-Ouchi, Cecile Agosta, Patrick Alexander, Xylar Asay-Davis, Alice Barthel, Thomas Bracegirdle, Richard Cullather, Denis Felikson, Xavier Fettweis, Jonathan Gregory, Tore Hatterman, Nicolas Jourdain, Peter Kuipers Munneke, Eric Larour, Christopher Little, Mathieu Morlinghem, Isabel Nias, Andrew Shepherd, Erika Simon, Donald Slater, Robin Smith, Fiammetta Straneo, Luke Trusel, Michiel van den Broeke, and Roderik van de Wal: 
Experimental protocol for sea level projections from ISMIP6 standalone ice sheet models, The Cryosphere, doi:10.5194/tc-2019-322, 2020.


About the data:
- The results are based on model output regridded conservatively to a 5x5 km regular ISMIP6 grid unless this is already the native grid. 
- The results are calculated over the ice-covered area of Greenland, map projection error corrected, ice sheet model specific densities taken into account.
- The contribution of peripheral glaciers and ice caps has been removed, by considering their area-coverage in each grid cell.
- The results for the projections 'exp*' are all calculated as differences to the control experiment ctrl_proj (suffix cr in filename for control removed).
- Results for ctrl_proj and historical are un-corrected (no suffix cr in filename).


Directory structure:
versionid
  groupname1
    modelname1
      expid
        scalars_mm_cr_GIS_groupname1_modelname1_expid.nc
        scalars_rm_cr_GIS_groupname1_modelname1_expid.nc
        scalars_zm_cr_GIS_groupname1_modelname1_expid.nc
...

Variables per file:

scalars_mm_cr_GIS ----------------- Greenland wide numbers 

oarea - assumed ocean area [m2]
rhof - model specific freshwater density [kg m-3]
rhoi - model specific ice density [kg m-3]
rhow - model specific ocean water density [kg m-3]

time - time, typically in days since X
iarea - Fraction of grid cell covered by land ice [1]
iareagr - Fraction of grid cell covered by grounded ice sheet
iareafl - Fraction of grid cell covered by ice sheet flowing over seawater

ivol - ice volume [m3]
ivolgr - grounded ice volume [m3]
ivolfl - floating ice volume [m3]
ivaf - ice volume above flotation [m3]

lim - ice mass [kg]
limgr - grounded ice mass [kg]
limfl - floating ice mass [kg]
limaf - ice mass above flotation [kg]

sle - sea-level equivalent mass [m] !! decreases with mass loss !! 
smb - spatially integrated surface mass balance anomaly [kg s-1]


scalars_rm_cr_GIS ----------------- IMBIE2-Rignot basins xx=[no,ne,se,sw,cw,nw]

oarea - assumed ocean area [m2]
rhof - model specific freshwater density [kg m-3]
rhoi - model specific ice density [kg m-3]
rhow - model specific ocean water density [kg m-3]

time - time, typically in days since X
ivaf_xx - ice volume above flotation [m3]
smb_xx - spatially integrated surface mass balance anomaly [kg s-1]
limaf_xx - ice mass above flotation [kg]
sle_xx - sea-level equivalent mass [m] !! decreases with mass loss !! 


scalars_zm_cr_GIS ----------------- IMBIE2-Zwally basins xx=[z11,z12,z13,z14,z21,z22,z31,z32,z33,z41,z42,z43,z50,z61,z62,z71,z72,z81,z82]

oarea - assumed ocean area [m2]
rhof - model specific freshwater density [kg m-3]
rhoi - model specific ice density [kg m-3]
rhow - model specific ocean water density [kg m-3]

time - time, typically in days since X
ivaf_xx - ice volume above flotation [m3]
smb_xx - spatially integrated surface mass balance anomaly [kg s-1]
limaf_xx - ice mass above flotation [kg]
sle_xx - sea-level equivalent mass [m] !! decreases with mass loss !! 

