function obs = update_2021(obs, tRes)

  % data for BHD station in NH
  nh = [-47.470, -47.440, -47.520, -47.50, -47.480, NaN, NaN, -47.610, -47.610, -47.62, -47.65, -47.57]

  % Data for AHR station in SH
  sh = [-47.43, -47.45, NaN, -47.46, -47.53, -47.59, -47.60, -47.60, NaN, NaN, NaN, NaN]
  

    if strcmp(tRes,'year') || strcmp(tRes,'YEAR') || strcmp(tRes,'yearly')
      obs.nh_ch4c13 = nanmean(nh);
      obs.sh_ch4c13 = nanmean(sh);
obs.nh_ch4c13_err(end) = obs.nh_ch4c13_err(end-1);
obs.sh_ch4c13_err(end) = obs.sh_ch4c13_err(end-1);
 else
     obs.nh_ch4c13(end-11:end) = nh
  obs.sh_ch4c13(end-11:end) = sh
end       
  end
  
