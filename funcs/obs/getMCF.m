%%% =======================================================================
%%% = getMCF.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Reads in the d13C observations.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): dataDir -- Directory containing the data.
%%% =  ( 2): reread  -- Structure that says if we'll re-read the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- A structure containing the observation information.
%%% =======================================================================

function [ out ] = getMCF( dataDir, reread )

%%% Diagnostic
fprintf('   * MCF\n');


%%% =======================================================================
%%% HAVE WE READ THIS DATA BEFORE?
%%% =======================================================================

%%% Build the filename
OutName = sprintf('%sobs/StoredData/mcf_%4i-%4i_%s-%s.mat',...
                  reread.dir,reread.sYear,reread.eYear,reread.tRes,reread.tAvg);

%%% Load pre-existing data file
if ~reread.flag
    % Check if a file exists
    if exist(OutName, 'file') == 2
        fprintf('   * LOADING OLD OBS STRUCTURE\n');
        load(OutName);
        return % Don't need to read the data
    end
end


%%% =======================================================================
%%% READ DATA
%%% =======================================================================

%%% Create the output structure
out = struct;
out.obs = struct;
out.tim = struct;
out.lat = struct;


%%% =======================================================================
%%% NOAA
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/mcf/NOAA/combined/',dataDir);

%%% Define the site names, header lengths, and latitudes
% Filename structure
fName = 'HATS_global_MC.txt';
nHDR  = 85;
% Specify the site names and latitudes
sNames = {'alt', 'sum', 'brw', 'mhd', 'thd', 'nwr', 'kum', 'mlo', 'smo', 'cgo', 'psa', 'spo'};
sLat   = [82.5,   72.6,   71.3,   53,    41,    40.052, 19.5,  19.5,  -14.3, -40.7, -64.6, -90];
sCol   = [  9,     11,    13,     15,    17,    19,     21,    23,    25,    27,    29,    31];

%%% Load the data
dat   = importdata(sprintf('%s%s',dataDirU,fName),' ',nHDR);
dat   = dat.data;
tDatO = datenum(dat(:,1),dat(:,2),ones(size(dat(:,1))));

%%% Read the data
for i = 1:length(sNames)
    % Get the column for this site the data
    yDat = dat(:,sCol(i));
    tDat = tDatO;
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_NOAA',sNames{i})) = yDat;
        out.tim.(sprintf('%s_NOAA',sNames{i})) = tDat;
        out.lat.(sprintf('%s_NOAA',sNames{i})) = sLat(i);
    end
end


%%% =======================================================================
%%% GAGE
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/mcf/GAGE/',dataDir);

%%% Define the site names, header lengths, and latitudes
% Filename structure
fNameS = '%s-gage.mon';
sNames = { 'CGO', 'MHD', 'RPB', 'SMO', 'ORG'};
sLat   = [-40.68, 53.33, 13.17,-14.23, 45.00];
nHDR   = [     6,     6,     6,     6,     6];

%%% Read the data
for i = 1:length(sNames)
    % Current filename
    fName = sprintf('%s%s',dataDirU,sprintf(fNameS,sNames{i}))
    % Load the data
    dat   = importdata(fName,' ',nHDR(i));
    dat   = dat.data;
    tDat  = datenum(dat(:,3),dat(:,2),ones(size(dat(:,1))));
    yDat  = dat(:,14);
    yDat(yDat == 0) = NaN;
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_GAGE',sNames{i})) = yDat;
        out.tim.(sprintf('%s_GAGE',sNames{i})) = tDat;
        out.lat.(sprintf('%s_GAGE',sNames{i})) = sLat(i);
    end
end


% <-- AJT (2019/05/01): Currently reading GAGE data as ALE data, commented for now. 
% %%% =======================================================================
% %%% GAGE
% %%% =======================================================================
% 
% %%% Append the directory onto the dataDir
% dataDirU = sprintf('%sobs/mcf/GAGE/',dataDir);
% 
% %%% Define the site names, header lengths, and latitudes
% % Filename structure
% fNameS = '%s-gage.mon';
% sNames = {  'CGO', 'RPB', 'SMO', 'ORG'};
% sLat   = [ -40.68, 13.17,-14.23, 45.00];
% nHDR   = [      6,     6,     6,     6];
% 
% %%% Read the data
% for i = 1:length(sNames)
%     % Current filename
%     fName = sprintf('%s%s',dataDirU,sprintf(fNameS,sNames{i}));
%     % Load the data
%     dat   = importdata(fName,' ',nHDR(i));
%     dat   = dat.data;
%     tDat  = datenum(dat(:,3),dat(:,2),ones(size(dat(:,1))));
%     yDat  = dat(:,14);
%     yDat(yDat == 0) = NaN;
%     % Remove NaNs
%     ind  = ~isnan(yDat) & ~isnan(tDat);
%     if sum(ind) > 0
%         tDat = tDat(ind);
%         yDat = yDat(ind);
%         % Put the data in a structure
%         out.obs.(sprintf('%s_ALE',sNames{i})) = yDat;
%         out.tim.(sprintf('%s_ALE',sNames{i})) = tDat;
%         out.lat.(sprintf('%s_ALE',sNames{i})) = sLat(i);
%     end
% end


%%% =======================================================================
%%% AGAGE
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/mcf/AGAGE/',dataDir);

%%% Define the site names, header lengths, and latitudes
% Filename structure for CH3CCl3 data
fNameS = 'AGAGE-GCMD_%s_ch3ccl3_mon.txt';  % Updated pattern
sNames = {'CGO', 'MHD', 'RPB', 'SMO', 'THD'};
sLat   = [-40.68, 53.33, 13.17, -14.23, 41.05];
nHDR   = [33, 33, 33, 33, 33];

%%% Read the data
for i = 1:length(sNames)
    % Current filename
    fName = fullfile(dataDirU, sprintf(fNameS, sNames{i}));
    % Check if the file exists
    if ~isfile(fName)
        warning('File does not exist: %s', fName);
        continue;
    end

    % Load the data
    dat   = importdata(fName, ' ', nHDR(i));
    dat   = dat.data;
    tDat  = datenum(dat(:,2), dat(:,3), ones(size(dat(:,1))));
    yDat  = dat(:,4);
    yDat(yDat == 0) = NaN;
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_AGAGE', sNames{i})) = yDat;
        out.tim.(sprintf('%s_AGAGE', sNames{i})) = tDat;
        out.lat.(sprintf('%s_AGAGE', sNames{i})) = sLat(i);
    end
end


%%% =======================================================================
%%% SAVE THIS OBSERVATION FILE
%%% =======================================================================

%%% Save the structure
fprintf('   * SAVING OBS STRUCTURE\n');
if exist(OutName, 'file') == 2
    delete(OutName);
end
save(OutName,'out');

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================