
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
sYear = 2000;
eYear = 2100;
%eYear = 2100;
tRes  = 'year';     % Can be 'year' or 'month' (year preferred)
tAvg  = 'month';     % Smooth the observations
St    = getTime(sYear,eYear,tRes); % Time vector
nT    = length(St);

%%% Execute in parallel?
run_parallel = false;
nWorkers     = 4;
setupParallel(run_parallel,nWorkers);

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
global k_co_flag use_strat interactive_OH use_other_sinks % Other
global use_N2O_strat
% Plotting flags
ftype           = 'pdf';    % Type of plots to make? (eps, pdf, tif, or png)
plot_prior      = false;     % Plot the prior?
plot_raw        = false;    % Plot the raw observations?
plot_old_cmaes  = false;    % Plot an old CMA-ES solution (false means run a new one)
% General flags
use_strat       = true;     % Use a stratosphere?
interactive_OH  = true;     % Allow OH feedbacks?
use_other_sinks = false;     % Use non-OH sinks?
% Linear inversion flags
det_linear      = true;     % Use a linear deterministic inversion?
fixedCH4        = false;    % Use fixed methane emissions
fixedOH         = true;    % Use fixed OH anomalies
onlyCH4         = false;    % Only invert for methane emissions
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
tau_TS = 2.5 * ones(nT,1); % years
if ~use_strat
    % Set this to something high, Inf results in trouble:
    %tau_TS(:) = Inf; % No exchange with stratosphere
    tau_TS(:) = 1e4;
end

%%% Arbitrary reactions with OH
% CF Needed to adapt NH as there would otherwise be a rather large IH
% difference in OH
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

ems_no_strat = ems
% Convert the structure to a matrix
ems = assembleEms(ems);

%%% Run the box model
params = getParameters(St); % Only need to do this once
%params.tau_NS = 0.1;

% set the interhemispheric exchange time to infinite so we deal with only 1 box
%params.tau_NS_strat=1e99;
%%params.odeOpts.MinStep=1;

% NN: Load observations for IC
%load('obs');
IC = params.IC;
%IC(1) = 1910;
%IC(2) = 1860;


% CO
%IC(13) = 87;
%IC(14) = 53;

%%% Define new emissions 
ems_1x = ems;

use_strat = true; 
interactive_OH = true;
ems_1x(:,15) = 3*ones(nT,1);


out_strat = boxModel_wrapper(St,ems_1x,IC,params);



ems_no_strat = ems;
%%% Strat-trop exchange
ems_no_strat(:,15) = 3.0* ones(nT,1); % years
if ~use_strat
    % Set this to something high, Inf results in trouble:
    %tau_TS(:) = Inf; % No exchange with stratosphere
    ems_no_strat(:,15) = 1e4*ones(nT,1);
end




out_no_strat = boxModel_wrapper(St,ems_no_strat,IC,params);
%save('steady_state_test')
time = [1:length(St)];
  lw = 3;
 
fig1 =figure(1);
ax1 = subplot(4,1,1);
plot(time , out_strat.nh_ch4, 'p', 'linewidth',lw)
%xlabel('years');
ylabel('ppb');
title('CH4 NH Trop');


ax2 = subplot(4,1,2)
plot(time, out_strat.sh_ch4, 'p', 'linewidth',lw);
title('SH Trop')
%xlabel('years')
ylabel('ppb')


ax3 = subplot(4,1,3)
plot(time , out_strat.nh_ch4_strat, 'p', 'linewidth',lw)
%xlabel('years');
ylabel('ppb');
title('NH Strat');

ax4 = subplot(4,1,4)
plot(time, out_strat.sh_ch4_strat, 'p', 'linewidth',lw);
title('SH Strat')
xlabel('years')
ylabel('ppb')
linkaxes([ax1 ax2 ax3 ax4],'x');
saveas(figure(1), 'strat_forward_modeltest_CH4.pdf', 'pdf')



fig2 =figure();
 ax1 = subplot(4,1,1);
plot(time , out_strat.nh_n2o, 'p', 'linewidth',lw)
%xlabel('years');
ylabel('ppb');
title('N2O NH Trop');


ax2 = subplot(4,1,2)
plot(time, out_strat.sh_n2o, 'p', 'linewidth',lw);
title('SH Trop')
%xlabel('years')
ylabel('ppb')


ax3 = subplot(4,1,3)
plot(time , out_strat.nh_n2o_strat, 'p', 'linewidth',lw)
%xlabel('years');
ylabel('ppb');
title('NH Strat');

ax4 = subplot(4,1,4)
plot(time, out_strat.sh_n2o_strat, 'p', 'linewidth',lw);
title('SH Strat')
xlabel('years')
ylabel('ppb')
linkaxes([ax1 ax2 ax3 ax4],'x');
saveas(figure(2), 'strat_forward_modeltest_N2O.pdf', 'pdf')
 
