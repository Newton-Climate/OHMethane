%%% Header
fprintf('\n ***********************************\n')
fprintf(' *** STARTING GLOBAL 2-BOX MODEL ***\n')
fprintf(' ***********************************\n')

%%% Define the directories
baseDir = fullfile(pwd, '../');
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
sYear = 1984;
eYear = 2022;
tRes  = 'year';     % Can be 'year' or 'month' (year preferred)
tAvg  = 'month';     % Smooth the observations
St    = getTime(sYear,eYear,tRes); % Time vector
nT    = length(St);

%%% Export variables to mat file
export_data = true; % do we want to export data to data_filename.mat?
data_filename  = 'monthly_obs.csv';


%%% For reading the observations
% Do we want to reread the raw data?
reread.flag  = false;
% Other flags for re-reading
reread.sYear = sYear;
reread.eYear = eYear;
reread.tRes  = tRes;
reread.tAvg  = tAvg;
reread.dir   = dataDir;


%%% grab the data and put them into a struct 
    ch4_obs     = getCH4(dataDir,reread);      % CH4 observations (ppb)
    ch4c13_obs  = getCH4C13(dataDir,reread);   % delta13C observations (permil)
    ch4h2_obs   = getCH4H2(dataDir,reread);    % deltaD observations (permil)
    mcf_obs     = getMCF(dataDir,reread);      % Methylchloroform observations (ppt)
    n2o_obs     = getN2O(dataDir,reread);      % N2O observations (ppb)
    c2h6_obs    = getC2H6(dataDir,reread);     % Ethane observations (ppt)
    co_obs      = getCO(dataDir,reread);       % carbon monoxide observations (ppb)
    o3strat_obs = getO3strat(dataDir,reread);  % Stratospheric ozone observations (DU)



%%% Diagnostics (check the raw data)
plot_raw = false;
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

% combine obs into one struct and perform hemispheric averaging 
obs = makeObs(St,tAvg,ch4_obs,ch4c13_obs,mcf_obs,n2o_obs,c2h6_obs,co_obs,dataDir,reread);

%%% for H2O strat data 
    stratH2O_global_fh = fullfile(dataDir, 'obs/h2o_strat/global_strat_obs_monthly.nc');
    H2O_global_strat = read_strat_obs(stratH2O_global_fh, St, tRes);
    stratH2O_tropical_fh = fullfile(dataDir, 'obs/h2o_strat/tropical_strat_obs_monthly.nc');
    H2O_tropical_strat = read_strat_obs(stratH2O_tropical_fh, St, tRes);
obs = makeObs(St,tAvg,ch4_obs,ch4c13_obs,mcf_obs,n2o_obs,c2h6_obs,co_obs,dataDir,reread);
obs.h2o_global_strat = H2O_global_strat.obs;
obs.h2o_tropical_strat = H2O_tropical_strat.obs;

%%% read N2O strat data 
obs = readN2OStrat(fullfile(dataDir, 'obs/N2O_strat/N2O_lowerstrat_obs.nc'), St, obs);



% convert to table and export to csv
T = struct2table(obs); % Convert struct to table
T.timestamp = datetime( double(St), 'ConvertFrom', 'datenum');    

% make timestamp first column in csv file 
T = movevars(T, 'timestamp', 'Before', T.Properties.VariableNames{1});

writetable(T, data_filename); % Write table to CSV file