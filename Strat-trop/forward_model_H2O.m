%%% =======================================================================
%%% = DriverScript.m
%%% = Alex Turner
%%% = Originally created on 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES:
%%% =
%%% = This is the driver script for the 2-box model methane inversion.
%%% = There are currently two different inversions implemented: (1) a
%%% = linear or non-linear deterministic inversion following Rodgers (2000)
%%% = and (2) an inversion using the non-linear Covariance Matrix Adaptation
%%% = Evolution Strategy (CMA-ES).  Case (1) requires us to compute gradients
%%% = and only allows for Gaussian errors.  Case (2) is  a stochastic method
%%% = that automatically tunes the proposal distribution to improve the sampling,
%%% = however it does not provide error statistics that are consistent with
%%% = the distributions.  Case (2) also allows one to specify non-analytic
%%% = distributions (e.g., bounded Gaussians or uniform distributions).
%%% =======================================================================


%%
%%% =======================================================================
%%% 1. Initialize
%%% =======================================================================
profile off
%%% Clear the MatLab space
clf
clear all
close all
clc
diary run1.txt
diary on

%%% Header
fprintf('\n ***********************************\n')
fprintf(' *** STARTING GLOBAL 2-BOX MODEL ***\n')
fprintf(' ***********************************\n')

%%% Define the directories
baseDir = pwd;
utilDir = sprintf('%s/funcs/', baseDir);
dataDir = sprintf('%s/data/',  baseDir);
outDir  = sprintf('%s/output/',baseDir);

%%% Add the utility functions
addpath(utilDir);
addpath(sprintf('%s/obs',               utilDir));
addpath(sprintf('%s/ems',               utilDir));
addpath(sprintf('%s/model',             utilDir));
addpath(sprintf('%s/util',              utilDir));
addpath(sprintf('%s/plot',              utilDir));
addpath(sprintf('%s/inv',               utilDir));
addpath(sprintf('%s/inv/deterministic', utilDir));
addpath(sprintf('%s/inv/stochastic',    utilDir));

%%% Define the time period
sYear = 1980;
eYear = 2080;
%eYear = 2100;
tRes  = 'year';     % Can be 'year' or 'month' (year preferred)
tAvg  = 'year';     % Smooth the observations
St    = getTime(sYear,eYear,tRes); % Time vector
nT    = length(St);

%%% Export variables to mat file
export_data = true; % do we want to export data to data_filename.mat?
data_filename  = 'synthetic_inversion_H2O';

disp('** Beginning the Synthetic Inversion Test **');



%%% Execute in parallel?
run_parallel = false;
if run_parallel
    nWorkers     = 4;
    setupParallel(run_parallel,nWorkers);
end


%%% What kind of inversions do we want to do?
do_deterministic = true;    % Rodgers (2000)
do_cmaes         = false;    % Covariance Matrix Adaptation Evolution Strategy

%%% For reading the observations
% Do we want to reread the raw data?
reread.flag  = false;
% Other flags for re-reading
reread.sYear = sYear;
reread.eYear = eYear;
reread.tRes  = tRes;
reread.tAvg  = tAvg;
reread.dir   = dataDir;

%%% Other options and flags
% Use globals for some flags
global fixedCH4 fixedOH onlyCH4 onlyMCF schaefer          % Linear inversion
global k_mcf_flag smooth_MCF set_MCF_EMS MCF_EMS_val      % Methyl Chloroform
global k_co_flag use_strat interactive_OH use_other_sinks ignoreCO % Other
global fit_strat_exchange no_temporal_correlation % temporal correlation flag
% Plotting flags
ftype           = 'pdf';    % Type of plots to make? (eps, pdf, tif, or png)
plot_prior      = false;     % Plot the prior?
plot_raw        = false;    % Plot the raw observations?
plot_old_cmaes  = false;    % Plot an old CMA-ES solution (false means run a new one)
% General flags
use_strat       = true;     % Use a stratosphere?
fit_strat_exchange = false;
interactive_OH  = true;     % Allow OH feedbacks?
use_other_sinks = false;     % Use non-OH sinks?
% Linear inversion flags
use_other_sinks = false;     % Use non-OH sinks?
% Linear inversion flags
det_linear      = false;     % Use a linear deterministic inversion?
fixedCH4        = false;    % Use fixed methane emissions
fixedOH         = false;    % Use fixed OH sources
onlyCH4         = false;    % Only invert for methane emissions
ignoreCO = true; % keep CO emissions fixed
onlyMCF         = false;    % Only invert for MCF emissions
schaefer        = false;    % Case that is most similar to Schaefer et al.
% MCF sensitivity test flags
k_co_flag       = true;     % Use k_CO that AJT derived
k_mcf_flag      = true;     % Use k_MCF that AJT derived
smooth_MCF      = false;    % Smooth the MCF emissions with a 5-year filter?
set_MCF_EMS     = false;    % Set post-2000 emissions to a fixed value?
MCF_EMS_val     = 0.0;      % Fixed post-2000 MCF emissions value (Gg/yr)
reduce_MCFerr   = false;    % Reduce the errors in MCF observations?
MCF_ERR_val     = 2.0;      % Error in MCF observations (ppt)
% Flags for other tests to run
use_OH_stratMLO = false;    % Use the OH derived from MLO strat ozone?
use_Ed          = false;    % Use Ed Dlugokencky's hemispheric averages?
no_temporal_correlation = true; % do we want to get rid of temporal correlations?


%%% Set the seed for repeatability
rng('default');


%%
%%% =======================================================================
%%% 2. Load the obs
%%% =======================================================================

%%% Diagnostic
fprintf('\n *** LOADING THE OBSERVATIONS *** \n');

%%% Load the observations
% Structures with with three fields:
% - "obs":  Observations from each NOAA site (ppb)
% - "tim":  Julian date for the observation
% - "lat":  Latitude of the NOAA site
try % Add a try-catch statement in case the user hasn't downloaded the data
    ch4_obs     = getCH4(dataDir,reread);      % CH4 observations (ppb)
    ch4c13_obs  = getCH4C13(dataDir,reread);   % delta13C observations (permil)
    ch4h2_obs   = getCH4H2(dataDir,reread);    % deltaD observations (permil)
    mcf_obs     = getMCF(dataDir,reread);      % Methylchloroform observations (ppt)
    n2o_obs     = getN2O(dataDir,reread);      % N2O observations (ppb)
    c2h6_obs    = getC2H6(dataDir,reread);     % Ethane observations (ppt)
    co_obs      = getCO(dataDir,reread);       % carbon monoxide observations (ppb)
    o3strat_obs = getO3strat(dataDir,reread);  % Stratospheric ozone observations (DU)
catch % Some data is missing
    try % See if ethane is the only problem
        fprintf(' * SOME DATA IS MISSING\n');
        ch4_obs     = getCH4(dataDir,reread);
        ch4c13_obs  = getCH4C13(dataDir,reread);
        ch4h2_obs   = getCH4H2(dataDir,reread);
        mcf_obs     = getMCF(dataDir,reread);
        n2o_obs     = getN2O(dataDir,reread);
        co_obs      = getCO(dataDir,reread);
        o3strat_obs = getO3strat(dataDir,reread);
        c2h6_obs    = NaN;
    catch % Otherwise, set the observation structures to NaN
        fprintf(' * UNABLE TO READ OBSERVATIONS!\n');
        ch4_obs     = NaN;
        ch4c13_obs  = NaN;
        ch4h2_obs   = NaN;
        mcf_obs     = NaN;
        n2o_obs     = NaN;
        c2h6_obs    = NaN;
        co_obs      = NaN;
        o3strat_obs = NaN;
    end
end

%%% Make the observation structure
% Structure with 12 fields:
% - NH/SH CH4    obs & err (ppb)
% - NH/SH CH4C13 obs & err (permil)
% - NH/SH MCF    obs & err (ppt)
% - NH/SH N2O    obs & err (ppb)
% - NH/SH C2H6   obs & err (ppt)
% - NH/SH CO     obs & err (ppb)
%obs = makeObs(St,tAvg,ch4_obs,ch4c13_obs,mcf_obs,n2o_obs,c2h6_obs,co_obs,dataDir,reread);
%
% blow up CO error:
%obs.nh_co_err(:)=500;
%obs.sh_co_err(:)=500;
%%% Use Ed Dlugokencky's obs? (sensitivity test)
if use_Ed
    ajt_obs = obs;
    ed_obs  = getEdObs(dataDir,ajt_obs,St,tAvg,reread);
    obs     = ed_obs;
    plotEdObs(St,ajt_obs,ed_obs,sprintf('%s/%s/raw_EdObs.%s',outDir,tRes,ftype))
end

%%% Reduce MCF errors?  (sensitivity test)
if reduce_MCFerr
    obs.nh_mcf_err = min([obs.nh_mcf_err,MCF_ERR_val*ones(size(obs.nh_mcf_err))],[],2);
    obs.sh_mcf_err = min([obs.sh_mcf_err,MCF_ERR_val*ones(size(obs.sh_mcf_err))],[],2);
end

%%% Diagnostics (check the raw data)
if plot_raw
    deseasonalize  = true;
    plot_all_sites = false;
    plotAllObs(St,obs,ch4_obs,   tAvg, 'ch4' ,sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,ch4c13_obs,tAvg, 'd13C',sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    %plotAllObs(St,obs,ch4h2_obs, tAvg, 'dD',  sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,mcf_obs,   tAvg, 'mcf' ,sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,n2o_obs,   tAvg, 'n2o' ,sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,c2h6_obs,  tAvg, 'c2h6',sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,co_obs,    tAvg, 'co'  ,sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
end


%%
%%% =======================================================================
%%% 3. Load the emissions (all will be arrays with a length matching "St")
%%% =======================================================================

%%% Diagnostic
fprintf('\n *** LOADING THE EMISSIONS *** \n');

%%% Get the CH4 emissions
% Stucture with two fields
% - "nh": CH4 emissions from the Northern hemisphere (Tg/yr)
% - "sh": CH4 emissions from the Southern hemisphere (Tg/yr)
ch4_ems = getCH4ems(St,tRes,dataDir);

%%% Get the delta13C composition for NH/SH CH4 emissions
% Stucture with two fields
% - "nh": delta13C composition from the Northern hemisphere (permil)
% - "sh": delta13C composition from the Southern hemisphere (permil)
ch4c13_ems = getCH4C13ems(St,tRes,dataDir);

%%% Get the MCF emissions (assumed to be in NH only)
% Stucture with two fields
% - "prinn":     MCF emissions from Prinn (Gg/yr)
% - "mcculloch": MCF emissions from McCulloch (Gg/yr)
mcf_ems = getMCFems(St,tRes,dataDir);

%%% Get the N2O emissions
% Stucture with two fields
% - "nh": N2O emissions from the Northern hemisphere (Tg/yr)
% - "sh": N2O emissions from the Southern hemisphere (Tg/yr)
n2o_ems = getN2Oems(St,tRes,dataDir);

%%% Get the C2H6 emissions
% Stucture with two fields
% - "nh": C2H6 emissions from the Northern hemisphere (Tg/yr)
% - "sh": C2H6 emissions from the Southern hemisphere (Tg/yr)
c2h6_ems = getC2H6ems(St,tRes,dataDir);

%%% Get the OH emissions
% Stucture with two fields
% - "nh": OH emissions from the Northern hemisphere (Tg/yr)
% - "sh": OH emissions from the Southern hemisphere (Tg/yr)
oh_ems = getOHems(St,tRes,dataDir);

%%% Get the CO emissions
% Stucture with two fields
% - "nh": CO emissions from the Northern hemisphere (Tg/yr)
% - "sh": CO emissions from the Southern hemisphere (Tg/yr)
co_ems = getCOems(St,tRes,dataDir);
% somehow, the NH emissions are way too low:
%co_ems.nh(:) = 1400;

%%
%%% =======================================================================
%%% 4. Initialize the 2-box model
%%% =======================================================================

%%% Diagnostic
fprintf('\n *** RUN THE 2-BOX MODEL WITH PRIOR FLUXES *** \n');

%%% OH scaling factor
oh_scale.nh = ones(nT,1);
oh_scale.sh = ones(nT,1);
% Derive OH from the stratospheric ozone?
if use_OH_stratMLO
    OH_sensitivity  = 4.2;          % a 1% increase in strat O3 leads to a 4.2% decrease in OH (Murray et al., 2013)
    O3_site         = 'mlo_NOAA';   % which site to use?
    fDays           = 365.25*2;     % How long of a smoothing?
    tO3             = o3strat_obs.tim.(O3_site);
    yO3             = o3strat_obs.obs.(O3_site);
    yO3             = DeseasonalizeData(tO3,yO3,fDays);
    [tO3, yO3, ~]   = BlockAverage_AltError(tO3,yO3,ones(size(tO3)),365.25);
    oh_change       = yO3 / nanmean(yO3); % Convert strat O3 to OH change
    oh_change       = 1 ./ ((oh_change - 1) * OH_sensitivity + 1);
    yOH             = interp1(tO3,oh_change,St);
    yOH(isnan(yOH)) = 1;
    % Store this OH
    oh_scale.nh = yOH;
    oh_scale.sh = yOH;
end

%%% Strat-trop exchange
tau_TS = 5.0 * ones(nT,1); % years
if ~use_strat
    % Set this to something high, Inf results in trouble:
    %tau_TS(:) = Inf; % No exchange with stratosphere
    tau_TS(:) = 1e4;
end

%%% ArbitrÏary reactions with OH
% CF Needed to adapt NH as there would otherwise be a rather large IH
% difference in OH
f = 2.07;
%kX_NH = 1.0*ones(nT,1); % s^-1
%kX_SH = 1.3*ones(nT,1); % s^-1
kX_NH = 0.99*ones(nT,1); % s^-1 for 6600 tg/yr OH source
kX_SH = 1.23*ones(nT,1); % s^-1

%%% Structure of sources with 17 fields:
% - NH CH4 emissions
% - SH CH4 emissions
% - NH CH4C13 composition
% - SH CH4C13 composition
% - NH MCF emissions
% - SH MCF emissions
% - NH N2O emissions
% - SH N2O emissions
% - NH C2H6 emissions
% - SH C2H6 emissions
% - NH OH emissions
% - SH OH emissions
% - NH CO emissions
% - SH CO emissions
% - Strat-trop exchange
% - NH arbitrary OH reaction rate
% - SH arbitrary OH reaction rate
ems.nh_ch4    = ch4_ems.nh;
ems.sh_ch4    = ch4_ems.sh;
ems.nh_ch4c13 = ch4c13_ems.nh;
ems.sh_ch4c13 = ch4c13_ems.sh;
ems.nh_mcf    = mcf_ems.nh;
ems.sh_mcf    = mcf_ems.sh;
ems.nh_n2o    = n2o_ems.nh;
ems.sh_n2o    = n2o_ems.sh;
ems.nh_c2h6   = c2h6_ems.nh;
ems.sh_c2h6   = c2h6_ems.sh;
ems.nh_oh     = oh_ems.nh;
ems.sh_oh     = oh_ems.sh;
ems.nh_co     = co_ems.nh;
ems.sh_co     = co_ems.sh;
ems.tau_TS    = tau_TS;
ems.kX_NH     = kX_NH;
ems.kX_SH     = kX_SH;

% Convert the structure to a matrix
ems = assembleEms(ems);
ems_no_strat = ems;
%ems_no_strat(:,15) = 1e4 * ones(nT,1);
% Make the perturbations for emissions 
ems_pert = ems;
% fix strat-trop exchange variability to lower strat H2O variability
H2O_strat = read_strat_obs('global_strat_obs.nc', {'combinedanomfillanomh2oq', 'combinedeqfillh2oq'});
H2O_var_percent = H2O_strat{1} ./ mean(H2O_strat{2}, 'omitnan') + 1;
%H2O_var_percent = H2O_var_percent(end-18:end);



ems(:,15) = 3*ones(nT,1);
after_pert = 10;
start = nT-length(H2O_var_percent)+1 - after_pert;
stop = start + length(H2O_var_percent) - 1;
ems_pert(end-length(H2O_var_percent)+1:end, 15) = H2O_var_percent .* ems(end-length(H2O_var_percent)+1:end, 15);

%%% Run the box model
params = getParameters(St); % Only need to do this once
IC     = params.IC;         % Guess for the inital conditions
interactive_OH  = true;     % Allow OH feedbacks?

% Run to get steady state for IC for 2nd box model run (synthetic data)
out    = boxModel_wrapper(St,ems,IC,params);
IC(1:14) = [out.nh_ch4(end), out.sh_ch4(end), out.nh_ch4(end)*(1+out.nh_ch4c13(end)/1e3), out.sh_ch4(end)*(1+out.sh_ch4c13(end)/1e3), out.nh_mcf(end), out.sh_mcf(end), out.nh_n2o(end), out.sh_n2o(end), out.nh_c2h6(end), out.sh_c2h6(end), (out.nh_oh(end)/params.n_air)*1d9, (out.sh_oh(end)/params.n_air)*1d9, out.nh_co(end), out.sh_co(end)];
IC(15:end) = [out.nh_ch4_strat(end), out.sh_ch4_strat(end), out.nh_ch4_strat(end)*(1+out.nh_ch4c13_strat(end)/1e3), out.sh_ch4_strat(end)*(1+out.sh_ch4c13_strat(end)/1e3), out.nh_mcf_strat(end), out.sh_mcf_strat(end), out.nh_n2o_strat(end), out.sh_n2o_strat(end), out.nh_c2h6_strat(end), out.sh_c2h6_strat(end), (out.nh_oh_strat(end)/params.n_air)*1d9, (out.sh_oh_strat(end)/params.n_air)*1d9, out.nh_co_strat(end), out.sh_co_strat(end)];


% Run to get synthetic data 
out    = boxModel_wrapper(St,ems_pert,IC,params);

% assign new observations and errors
obs = out;

synthetic_error = 0.01; %  make error 1 percent of observation because we know the observations 

obs.nh_ch4_err = obs.nh_ch4*synthetic_error;
obs.sh_ch4_err = obs.sh_ch4*synthetic_error;
obs.nh_ch4c13_err = obs.nh_ch4c13 * synthetic_error;
obs.sh_ch4c13_err = obs.sh_ch4c13*synthetic_error;
obs.nh_mcf_err = obs.nh_mcf*synthetic_error;
obs.sh_mcf_err = obs.sh_mcf*synthetic_error;
obs.nh_n2o_err = obs.nh_n2o * synthetic_error;
obs.sh_n2o_err = obs.sh_n2o * synthetic_error;
obs.nh_c2h6_err = obs.nh_c2h6*synthetic_error;
obs.sh_c2h6_err = obs.sh_c2h6*synthetic_error;
obs.nh_co_err = obs.nh_co * synthetic_error;
obs.sh_co_err = obs.sh_co * synthetic_error;

if plot_prior
    plotNewObs(St,out,obs,sprintf('%s/%s/prior_%%s.%s',outDir,tRes,ftype));
    %writeData(St,obs,out,ems,IC,sprintf('%s/%s/prior_%%s.csv',outDir,tRes));
    %plotObs(St,out,obs,sprintf('%s/%s/prior_%%s.%s',outDir,tRes,ftype));
    %plotDrivers(St,ems,NaN*ems,sprintf('%s/%s/prior_%%s.%s',outDir,tRes,ftype),dataDir);
end


%%
%%% =======================================================================
%%% 5. Deterministic inversion (Rodgers, 2000)
%%% =======================================================================

if do_deterministic
    
    %%% Diagnostic
    fprintf('\n *** DETERMINISTIC INVERSION *** \n');
    
    %%% Invert

interactive_OH  = true;     % Allow OH feedbacks?



% Create an arbitrary guess for initial ems

% NN: run inversion with strat 
    [anal_soln,jacobian_ems,jacobian_IC,reltol,abstol, mati] = invert_methane(St,obs,ems_pert,IC,params,det_linear,run_parallel);
out_strat    = boxModel_wrapper(St,anal_soln{1},IC,params);



% NN: Let's now run our inversion without the stratosphere 
interactive_OH = true;
% model spin-up
out    = boxModel_wrapper(St,ems_no_strat,IC,params);
IC(1:14) = [out.nh_ch4(end), out.sh_ch4(end), out.nh_ch4(end)*(1+out.nh_ch4c13(end)/1e3), out.sh_ch4(end)*(1+out.sh_ch4c13(end)/1e3), out.nh_mcf(end), out.sh_mcf(end), out.nh_n2o(end), out.sh_n2o(end), out.nh_c2h6(end), out.sh_c2h6(end), (out.nh_oh(end)/params.n_air)*1d9, (out.sh_oh(end)/params.n_air)*1d9, out.nh_co(end), out.sh_co(end)];

det_linear = false;


    [anal_soln2,jacobian_ems2,jacobian_IC2,reltol2,abstol2, mati2] = invert_methane(St,obs,ems_no_strat,IC,params,det_linear,run_parallel);
out_no_strat    = boxModel_wrapper(St,anal_soln2{1},IC,params);

out_no_strat    = boxModel_wrapper(St,anal_soln2{1},IC,params);
strat_ems_nh = anal_soln{1}(:,1);
strat_ems_sh = anal_soln{1}(:,2);
strat_ems = strat_ems_nh + strat_ems_sh

no_strat_ems_nh = anal_soln2{1}(:,1);
no_strat_ems_sh = anal_soln2{1}(:,2);
no_strat_ems = no_strat_ems_nh + no_strat_ems_sh

true_ems_nh = ems_pert(:,1);
true_ems_sh = ems_pert(:,2);
true_ems = true_ems_nh + true_ems_sh

% get global average for concentrations 
out_strat.ch4 = (out_strat.nh_ch4 + out_strat.sh_ch4)/2;
out_strat.oh = (out_strat.nh_oh + out_strat.sh_oh)/2;

out_no_strat.ch4 = (out_no_strat.nh_ch4 + out_no_strat.sh_ch4)/2;
out_no_strat.oh = (out_no_strat.nh_oh + out_no_strat.sh_oh)/2;
obs.ch4 = (obs.nh_ch4 + obs.sh_ch4)/2;
obs.oh = (obs.nh_oh + obs.sh_oh)/2;

if export_data
    fprintf('Exporting all variables in this run to %s \n', data_filename)
    save(data_filename);
end


% plot the fitted concentrations 
clf
close all
figure(1)
subplot(211)
time = [1: length(St)];
time = time(5:end-5);
plot(time, obs.ch4(5:end-5), 'r-', time, out_noninteractive.ch4(5: end-5), 'g-', time, out_interactive.ch4(5 : end - 5), 'b-')
xlabel('years')
ylabel('ppb')
title('CH4 concentrations')

subplot(212)
plot(time, obs.oh(5:end-5), 'r-', time, out_noninteractive.oh(5: end-5), 'g-', time, out_interactive.oh(5 : end - 5), 'b-')
xlabel('years')
ylabel('molec/cm^3')
title('OH concentrations')



%NN: Let's stop the script here
return

%%% NN Plotting the results of the inversion here
fig1 =figure(1);
subplot(1,2,1);
time = [1980:2040];
plot(time , anal_soln{1}(:,1), 'p')
hold on
plot(time, anal_soln2{1}(:,1), 'b');
hold on
plot(time, ems(:,1) , 'p-')
hold on
xlabel('years');
ylabel('Tg');
title('NH Synthetic Emissions Test');
legend('Interactive Chemistry', 'Noninteractive Chemistry', 'True Emissions')

subplot(1,2,2);
diff1 = anal_soln{1}(:,1) - ems(:,1);
diff2 = anal_soln{1}(:,2) - ems(:,2);

plot(time , anal_soln{1}(:,1), 'p')
hold on
plot(time, anal_soln2{1}(:,1), 'b');
hold on
plot(time, ems(:,1) , 'p-')
hold on
title('SH Synthetic Emissions Test')
xlabel('years')
ylabel('Tg')
saveas(figure(1), 'synthetic_emissions_tests.png', 'png')
return

    
    %%% Plot the Jacobians
    %[jacobian_ems,jacobian_IC] = define_Jacobian( St, ems, IC, params, run_parallel );


    plotJacobian(St,jacobian_ems,tRes,sprintf('%s/%s/jacobian_%%s.%s',outDir,tRes,ftype));
    
    %%% Try plotting the solution
    ems_anal = anal_soln{1};
    IC_anal  = anal_soln{2};
    % Comment this out for now:
    out_anal1 = boxModel_wrapper(St,ems_anal,IC_anal,params);
    plotNewObs(St,out_anal,obs,sprintf('%s/%s/anal_%%s.%s',outDir,tRes,ftype));
    %writeData(St,obs,out_anal,ems_anal,IC_anal,sprintf('%s/%s/anal_%%s.csv',outDir,tRes));
    %plotObs(St,out_anal,obs,sprintf('%s/%s/anal_%%s.%s',outDir,tRes,ftype));
    %plotDrivers(St,ems_anal,ems,sprintf('%s/%s/anal_%%s.%s',outDir,tRes,ftype),dataDir);
    
end


%%
%%% =======================================================================
%%% 6. Invert with CMAES
%%% =======================================================================

if do_cmaes
    
    %%% Diagnostics
    fprintf('\n *** STARTING CMA-ES INVERSION *** \n');
    
    %%% Use a log-likelihood and get matrix dimensions
    use_log = true;
    nE      = size(ems,2);
    nI      = length(IC);
    
    %%% Set the function parameters for the box model
    fun_param.run_parallel = run_parallel;
    fun_param.p_prior      = @(ems,IC,input_param) define_prior(ems,IC,input_param);
    fun_param.p_like       = @(ems,IC,input_param) define_likelihood(ems,IC,input_param);
    fun_param.use_log      = use_log;
    fun_param.params       = params;
    fun_param.IC           = params.IC;
    fun_param.St           = St;
    fun_param.ch4_ems      = ch4_ems;
    fun_param.ch4c13_ems   = ch4c13_ems;
    fun_param.mcf_ems      = mcf_ems;
    fun_param.n2o_ems      = n2o_ems;
    fun_param.c2h6_ems     = c2h6_ems;
    fun_param.oh_scale     = oh_scale;
    fun_param.tau_TS       = tau_TS;
    fun_param.nT           = nT;
    fun_param.nE           = nE;
    fun_param.nI           = nI;
    fun_param.obs          = obs;
    
    %%% Set the options for CMAES
    CMAES_opts                   = cmaes;
    CMAES_opts.EvalParallel      = 'yes';
    CMAES_opts.SaveFilename      = sprintf('output/%s/cmaes_dat/variablescmaes.mat',tRes);
    CMAES_opts.LogFilenamePrefix = sprintf('output/%s/cmaes_dat/outcmaes',tRes);
    CMAES_opts.Resume            = 'no'; % Default is 'no'
    CMAES_opts.CMA.active        = 2;
    CMAES_opts.DiagonalOnly      = 100;
    %     % Full
    %     CMAES_opts.StopFunEvals      = 5000000;
    %     CMAES_opts.MaxIter           = 20000;
    %     CMAES_opts.Restarts          = 10;
    %     % Medium
    %     CMAES_opts.StopFunEvals      = 1000000;
    %     CMAES_opts.MaxIter           = 5000;
    %     CMAES_opts.Restarts          = 6;
    % Small
    CMAES_opts.StopFunEvals      = 10000;
    CMAES_opts.MaxIter           = 200;
    CMAES_opts.Restarts          = 2;
    
    %%% Get the starting point and standard deviations
    % Starting point
    if do_deterministic % Use the deterministic inversion as a starting point
        xstart = ems_anal;
        xstart(xstart(:,1)<0,1)   = 0+eps; % Enforce positive CH4 emissions
        xstart(xstart(:,2)<0,2)   = 0+eps; % Enforce positive CH4 emissions
        xstart(xstart(:,1)<0,3)   = 0+eps; % Enforce positive MCF emissions
        xstart(xstart(:,2)<0,4)   = 0+eps; % Enforce positive MCF emissions
        xstart(xstart(:,7)<0,7)   = 0+eps; % Enforce positive N2O emissions
        xstart(xstart(:,8)<0,8)   = 0+eps; % Enforce positive N2O emissions
        xstart(xstart(:,9)<0,9)   = 0+eps; % Enforce positive C2H6 emissions
        xstart(xstart(:,10)<0,10) = 0+eps; % Enforce positive C2H6 emissions
        xstart = assembleStateVector(xstart,IC_anal);
    else
        xstart = ems;
        xstart = assembleStateVector(xstart,IC);
    end
    % Standard deviations
    xsigma.ems       = ones(nT,nE);
    sig.ch4          =   1.000;     % Tg/yr
    sig.ch4c13       =   0.050;     % permil
    sig.mcf_nh       =   5.000;     % Gg/yr
    sig.mcf_sh       =   5.000;     % Gg/yr
    sig.n2o          =   0.500;     % Tg/yr
    sig.c2h6         = 100.000;     % Gg/yr
    sig.oh           =   0.005;     % OH scale factor
    sig.tau          =   0.500;     % 1/yr
    xsigma.IC        = [0.05,0.05,0.01,0.01,0.05,0.05,0.05,0.05,0.05,0.05,...
        0.05,0.05,0.01,0.01,0.05,0.05,0.05,0.05,0.05,0.05]; % Initial conditions
    xsigma.ems(:,1)  = xsigma.ems(:,1)  * sig.ch4;
    xsigma.ems(:,2)  = xsigma.ems(:,2)  * sig.ch4;
    xsigma.ems(:,3)  = xsigma.ems(:,3)  * sig.ch4c13;
    xsigma.ems(:,4)  = xsigma.ems(:,4)  * sig.ch4c13;
    xsigma.ems(:,5)  = xsigma.ems(:,5)  * sig.mcf_nh;
    xsigma.ems(:,6)  = xsigma.ems(:,6)  * sig.mcf_sh;
    xsigma.ems(:,7)  = xsigma.ems(:,7)  * sig.n2o;
    xsigma.ems(:,8)  = xsigma.ems(:,8)  * sig.n2o;
    xsigma.ems(:,9)  = xsigma.ems(:,9)  * sig.c2h6;
    xsigma.ems(:,10) = xsigma.ems(:,10) * sig.c2h6;
    xsigma.ems(:,11) = xsigma.ems(:,11) * sig.oh;
    xsigma.ems(:,12) = xsigma.ems(:,12) * sig.oh;
    xsigma.ems(:,13) = xsigma.ems(:,13) * sig.tau./params.YrToDay;
    xsigma           = assembleStateVector(xsigma.ems,xsigma.IC);
    
    %%% Use an old estimate?
    use_old = false;
    if use_old
        fprintf('   * PLOTTING THE OLD CMA-ES RESULTS\n')
        fname = sprintf('./%s',CMAES_opts.SaveFilename);
        if (exist(fname,'file') == 2)
            load(fname);
            xstart = bestever.x;
        end
    end
    
    % Run the CMA-ES inversion or just plot an old one?
    if plot_old_cmaes
        [cmaes_res.ems,cmaes_res.IC] = disassembleStateVector(bestever.x,nT,nE,nI);
    else
        %%% Invert with CMAES
        % INPUTS
        %  - fitfun:  name of objective/fitness function
        %  - xstart:  objective variables initial point, determines N
        %  - insigma: initial coordinate wise standard deviation(s)
        %  - inopts:  options struct
        % OUTPUTS
        %  - xmin:      minimum search point of last iteration
        %  - fmin:      function value of xmin
        %  - counteval: number of function evaluations done
        %  - stopflag:  stop criterion reached
        %  - outCMAES:  struct with various histories and solutions
        %  - bestever:  struct containing overall best solution (for convenience)
        [xmin,fmin,counteval,stopflag,outCMAES,bestever] = cmaes('cmaes_fun_eval', xstart, xsigma, CMAES_opts, fun_param);
        [cmaes_res.ems,cmaes_res.IC]                     = disassembleStateVector(bestever.x,nT,nE,nI);
        %[cmaes_res.ems,cmaes_res.IC]                     = disassembleStateVector(xmin,nT,nE,nI);
    end
    
    %%% Plot the best one
    ems_best = cmaes_res.ems;
    IC_best  = cmaes_res.IC;
    out_best = boxModel_wrapper(St,ems_best,IC_best,params);
    plotNewObs(St,out_best,obs,sprintf('%s/%s/cmaes_%%s.%s',outDir,tRes,ftype));
    %writeData(St,obs,out_best,ems_best,IC_best,sprintf('%s/%s/cmaes_%%s.csv',outDir,tRes));
    %plotObs(St,out_best,obs,sprintf('%s/%s/cmaes_%%s.%s',outDir,tRes,ftype));
    %plotDrivers(St,ems_best,ems,sprintf('%s/%s/cmaes_%%s.%s',outDir,tRes,ftype),dataDir);
    
end

if export_data
    fprintf('Exporting all variables in this run to %s \n', data_filename)
    save(data_filename);
end


%%
%%% Finished simulation
fprintf('\n ***********************************\n')
fprintf(' ***            DONE!            ***\n')
fprintf(' ***********************************\n\n')



