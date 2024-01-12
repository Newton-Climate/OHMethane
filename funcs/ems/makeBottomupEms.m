function ems = makeBottomupEms(sYear, eYear)


ems = struct(); % allocate output struct

% read the biomass burning emissions from GFED
edgar_file = 'biomass_burning_gfed.txt';
fireEms = readFireEms(edgar_file, sYear, eYear);
ems.fires_nh = fireEms.nh_fires;
ems.fires_sh = fireEms.sh_fires;
% EDGAR anthropogenic emissions 
edgar_file = 'EDGAR_V7_Hemispheric.npy';
edgarEms = getEdgarEms(edgar_file, sYear, eYear);
anthro_sectors = fieldnames(edgarEms);
for i=1:length(anthro_sectors)
  sector = char(anthro_sectors(i));
  ems.(sector) = edgarEms.(sector);
end

% read wetland emissions from the wetcharts database 
wetland_file = 'wetland_ems.nc';
wetlandEms = readWetlandEms(wetland_file, sYear, eYear);
ems.wetlands_nh = wetlandEms.nh_wetlands;
ems.wetlands_sh = wetlandEms.sh_wetlands;
end

