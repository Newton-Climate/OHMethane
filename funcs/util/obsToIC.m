function IC = obsToIC(obs_baseline, t_baseline, params)

IC = zeros(28,1);
nh_oh = (1e6/params.n_air)*1d9;
sh_oh = (1e6/params.n_air)*1d9;
nh_ch4c13 = obs_baseline.nh_ch4(t_baseline)*(1+obs_baseline.nh_ch4c13(t_baseline)/1e3);
sh_ch4c13 = obs_baseline.sh_ch4(t_baseline)*(1+obs_baseline.sh_ch4c13(t_baseline)/1e3);


% set IC to the first year of observed values
% [OH] comes from MCF in case5.mat
IC(1:14) = [obs_baseline.nh_ch4(t_baseline), obs_baseline.sh_ch4(t_baseline), nh_ch4c13, sh_ch4c13, obs_baseline.nh_mcf(t_baseline), obs_baseline.sh_mcf(t_baseline), obs_baseline.nh_n2o(t_baseline), obs_baseline.sh_n2o(t_baseline), obs_baseline.nh_c2h6(t_baseline), obs_baseline.sh_c2h6(t_baseline), nh_oh, sh_oh, obs_baseline.nh_co(t_baseline), obs_baseline.sh_co(t_baseline)];

IC(21:22) = [obs_baseline.nh_n2o_strat(t_baseline), obs_baseline.sh_n2o_strat(t_baseline)];
end
