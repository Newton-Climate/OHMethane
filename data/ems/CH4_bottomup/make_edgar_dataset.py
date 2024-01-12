# Read EDGARv7 for NEWTON
dir = 'Bottom_up/ANTHROP/EDGARv7.0/'
ff_list = ['CHE','ENE','FFF','IND','IRO','PRO',
          'RCO','REF_TRF','TNR_Aviation_CDS','TNR_Aviation_CRS',
          'TNR_Aviation_LTO','TNR_Other','TNR_Ship','TRO_noRES',] 
waste_list = ['SWD_INC', 'SWD_LDF', 'WWT']
rice_list = ['AGS']
animal_list = ['ENF','MNM']
biofuel_list = ['AWB']
Anthro = np.zeros((5, (2022-1970), 2)) 
for yyyy in range(1970, 2022):
    for var in ff_list:
        f = Dataset(dir+var+'_nc/'+'v7.0_FT2021_CH4_'+str(yyyy)+'_'+var+'.0.1x0.1.nc')
        va = f.variables['emi_ch4'][:] #kg m-2 s-1 
        Anthro[0, yyyy-1970, 0] += np.nansum( va[:900, :]*area01[:900, :])*3600*24*365*1.e-9     
        Anthro[0, yyyy-1970, 1] += np.nansum( va[900:, :]*area01[900:, :])*3600*24*365*1.e-9 
    for var in waste_list:
        f = Dataset(dir+var+'_nc/'+'v7.0_FT2021_CH4_'+str(yyyy)+'_'+var+'.0.1x0.1.nc')
        va = f.variables['emi_ch4'][:] #kg m-2 s-1 
        Anthro[1, yyyy-1970, 0] += np.nansum( va[:900, :]*area01[:900, :])*3600*24*365*1.e-9     
        Anthro[1, yyyy-1970, 1] += np.nansum( va[900:, :]*area01[900:, :])*3600*24*365*1.e-9  
        
    var = 'AGS'    
    f = Dataset(dir+var+'_nc/'+'v7.0_FT2021_CH4_'+str(yyyy)+'_'+var+'.0.1x0.1.nc')
    va = f.variables['emi_ch4'][:] #kg m-2 s-1 
    Anthro[2, yyyy-1970, 0] += np.nansum( va[:900, :]*area01[:900, :])*3600*24*365*1.e-9     
    Anthro[2, yyyy-1970, 1] += np.nansum( va[900:, :]*area01[900:, :])*3600*24*365*1.e-9  
        
    for var in animal_list:
        f = Dataset(dir+var+'_nc/'+'v7.0_FT2021_CH4_'+str(yyyy)+'_'+var+'.0.1x0.1.nc')
        va = f.variables['emi_ch4'][:] #kg m-2 s-1 
        Anthro[3, yyyy-1970, 0] += np.nansum( va[:900, :]*area01[:900, :])*3600*24*365*1.e-9     
        Anthro[3, yyyy-1970, 1] += np.nansum( va[900:, :]*area01[900:, :])*3600*24*365*1.e-9  
        
    var = 'AWB'
    f = Dataset(dir+var+'_nc/'+'v7.0_FT2021_CH4_'+str(yyyy)+'_'+var+'.0.1x0.1.nc')
    va = f.variables['emi_ch4'][:] #kg m-2 s-1 
    Anthro[4, yyyy-1970, 0] += np.nansum( va[:900, :]*area01[:900, :])*3600*24*365*1.e-9     
    Anthro[4, yyyy-1970, 1] += np.nansum( va[900:, :]*area01[900:, :])*3600*24*365*1.e-9     
    
np.save('EDGAR_V7_Hemispheric.npy', Anthro)





5:09