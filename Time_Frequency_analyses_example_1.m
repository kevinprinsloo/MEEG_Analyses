%   FT_Time-frequency analysis
%
%   Summary:
%   In FieldTrip, the time-frequency decomposition is computed by the function ft_freqanalysis,
%   specifying one of the supported algorithms in cfg.method. The method mtmconvol implements a STFT,
%   whereas the method wavelet is based on Morlet wavelets (a Morlet wavelet is a complex sinusoid
%   tapered with a Gaussian function).
%
%   Status:
%   Under development
%
%   Notes:
%
%
%   References:
%      [1] https://www.fieldtriptoolbox.org/tutorial/timefrequencyanalysis/
%      [2] Tallon-Baudry C, Bertrand O, Delpuech C, Permier J. Oscillatory gamma-band (30-70 Hz)
%          activity induced by a visual search task in humans. J Neurosci. 1997 Jan 15;17(2):722-34.
%      [3] Bruns A, Fourier-, Hilbert-and wavelet-based signal analysis: are they really different approaches?.
%         J.Neurosci.Meth. 2004 Aug;137(2):321-332
%
%   Author(s): Kevin Prinsloo <kevinvago@gmail.com>
%              Cognitive Neurophysiology Laboratory,
%              University of Rochester, NY
%              May 2021; Last Revision: 23-May-2021
%
%   Copyright (c) 2021 Kevin Prinsloo 

%% Initialise Path Variables

dataPath = 'C:\Users\'; % It depends on where the
                                                   % wikiDatasets folder is
                                                   % stored
fileNameDescription = 'Task.set';
 
subjects = {'S01' 'S02' 'S03' ...
            'S04' 'S05'}; % the variable subjects contains the list of 
                          % subject-specific directories
nSubjs = numel(subjects);

%% Time-frequency domain analysis

% Standard tones
for i = 1:nSubjs
    fileName = 'dataStd.mat';
    dataSet = fullfile(dataPath,subjects{i},fileName);
    load(dataSet); % loading the variable 'data'
    
    % Computing time-frequency representation for each trial
    cfg = [];
    cfg.method = 'wavelet';
    cfg.foi = 1:1:40; % 1 = frequency resolution = 1/trial duration
    cfg.toi = -0.05:0.05:0.5;
    cfg.width = 3; % width of the wavelet ( = number of cycles). 
                   % Small values increase the temporal resolution at 
                   % the expense of frequency resolution
    cfg.keeptrials = 'yes';
    cfg.channel = {'all'};
    cfg.outputfile = fullfile(dataPath,subjects{i},'timeFreqStd');
    ft_freqanalysis(cfg,data);
    clear data
end
 
timeFreqAvgStd = cell(1,nSubjs);
for i = 1:nSubjs
    fileName = 'timeFreqStd.mat';
    dataSet = fullfile(dataPath,subjects{i},fileName);
    load(dataSet); % loading the variable 'freq'
    
    % Computing the average over trials
    cfg = [];
    timeFreqAvgStd{i} = ft_freqdescriptives(cfg,freq);
    clear freq
end
 
bslnTimeFreqAvgStd = cell(1,nSubjs);
for i = 1:nSubjs
    fileName = 'timeFreqStd.mat';
    dataSet = fullfile(dataPath,subjects{i},fileName);
    load(dataSet); % loading the variable 'freq'
    
    % Computing the change relative to a baseline
    cfg = [];
    cfg.baselinetype = 'relchange'; % 'absolute','relative','relchange' or 'db'
    cfg.baseline = [-0.05 0];
    bslnTimeFreqAvgStd{i} = ft_freqbaseline(cfg,freq);
    
    % Computing the average over trials
    cfg = [];
    bslnTimeFreqAvgStd{i} = ft_freqdescriptives(cfg,bslnTimeFreqAvgStd{i});
    clear freq
end
 
% Deviant tones
for i = 1:nSubjs
    fileName = 'dataDev.mat';
    dataSet = fullfile(dataPath,subjects{i},fileName);
    load(dataSet); % loading the variable 'data'
    
    % Computing time-frequency representation for each trial
    cfg = [];
    cfg.method = 'wavelet';
    cfg.foi = 1:1:40; % 1 = frequency resolution = 1/trial duration
    cfg.toi = -0.05:0.05:0.5;
    cfg.width = 3; % width of the wavelet ( = number of cycles). 
                   % Small values increase the temporal resolution at 
                   % the expense of frequency resolution
    cfg.keeptrials = 'yes';
    cfg.channel = {'all'};
    cfg.outputfile = fullfile(dataPath,subjects{i},'timeFreqDev');
    ft_freqanalysis(cfg,data);
    clear data
end
 
timeFreqAvgDev = cell(1,nSubjs);
for i = 1:nSubjs
    fileName = 'timeFreqDev.mat';
    dataSet = fullfile(dataPath,subjects{i},fileName);
    load(dataSet); % loading the variable 'freq'
    
    % Computing the average over trials
    cfg = [];
    timeFreqAvgDev{i} = ft_freqdescriptives(cfg,freq);
    clear freq
end
 
bslnTimeFreqAvgDev = cell(1,nSubjs);
for i = 1:nSubjs
    fileName = 'timeFreqDev.mat';
    dataSet = fullfile(dataPath,subjects{i},fileName);
    load(dataSet); % loading the variable 'freq'
    
    % Computing the change relative to a baseline
    cfg = [];
    cfg.baselinetype = 'relchange'; % 'absolute','relative','relchange' or 'db'
    cfg.baseline = [-0.05 0];
    bslnTimeFreqAvgDev{i} = ft_freqbaseline(cfg,freq);
    
    % Computing the average over trials
    cfg = [];
    bslnTimeFreqAvgDev{i} = ft_freqdescriptives(cfg,bslnTimeFreqAvgDev{i});
    clear freq
end

%% Grand average time-frequency data
% Standard tones
cfg = [];
gAvgTimeFreqStd = ft_freqgrandaverage(cfg,timeFreqAvgStd{:});
gAvgBslnTimeFreqStd = ft_freqgrandaverage(cfg,bslnTimeFreqAvgStd{:});
 
% Deviant tones
cfg = [];
gAvgTimeFreqDev = ft_freqgrandaverage(cfg,timeFreqAvgDev{:});
gAvgBslnTimeFreqDev = ft_freqgrandaverage(cfg,bslnTimeFreqAvgDev{:});

%% Plotting time-frequency analysis
% Defining channel layout
templateLayout = 'EEG1005.lay'; % one of the template layouts included
                                % in FieldTrip
cfg = [];
cfg.layout = which(templateLayout); 
layout = ft_prepare_layout(cfg);
 
% Increasing layout width and height
if strcmp(templateLayout,'EEG1005.lay')
    layout.width = 0.07 * ones(length(layout.width),1);
    layout.height = 0.04 * ones(length(layout.height),1);
end
 
% Multi plot
cfg = [];
cfg.layout = layout; %'easycapM1.mat';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
cfg.showlabels = 'yes';
cfg.fontsize = 10;
cfg.comment = newline; % New Implimentation on Matlab sprintf('\n')
cfg.xlim = [-0.05 0.5];
cfg.ylim = [8 40]; % we cannot analyze delta and theta band because 
                   % the trial lenght is too short
cfg.zlim = [-6 6]; 
cfg.colorbar = 'yes'; % or 'southoutside'
cfg.colormap = jet; % 'parula' or 'jet'
figure
ft_multiplotTFR(cfg,gAvgBslnTimeFreqStd);
title('Grand average - Standard (relchange)');
figure
cfg.zlim = [-6 6];
ft_multiplotTFR(cfg,gAvgBslnTimeFreqDev);
title('Grand average - Deviant (relchange)');
 
% Plotting single average across electrodes of interest
cfg = [];
cfg.ylim = [8 40];
cfg.zlim = [-6 6];
cfg.channel = {'Fz','FCz','Cz'};
cfg.colormap = jet;
figure;
ft_singleplotTFR(cfg,gAvgBslnTimeFreqStd);
set(gca,'Fontsize',20);
title('Mean over fronto-central electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('frequency (Hz)');
c = colorbar;
c.LineWidth = 1;
c.FontSize  = 18;
title(c,'Rel. change');

%% Statistical analysis of time-frequency representations at the group level
% Creating a neighbourhood structure to define how spatial clusters are formed
fileName = 'dataStd.mat';
dataSet = fullfile(dataPath,subjects{1},fileName);
load(dataSet); % loading the variable 'data'
 
templateElec = 'standard_1005.elc';
elec = ft_read_sens(which(templateElec));
 
idx = ismember(elec.label,data.label);
elec.chanpos = elec.chanpos(idx,:);
elec.chantype = elec.chantype(idx);
elec.chanunit = elec.chanunit(idx);
elec.elecpos = elec.elecpos(idx,:);
elec.label = elec.label(idx);
 
cfg = [];
cfg.method = 'distance'; %'template' or 'distance' or 'triangulation'
cfg.neighbourdist = 80; % maximum distance between neighbouring sensors
                        % (only for 'distance')
cfg.feedback = 'yes'; % show a neighbour plot
neighbours = ft_prepare_neighbours(cfg,elec);
 
% Preparing the design matrix for the statistical evaluation
% For within-subjects analysis, the design matrix contains two rows
design = [1:nSubjs 1:nSubjs; ones(1,nSubjs) ones(1,nSubjs)*2]; 
 
% Test the null-hypothesis of exchangeability between conditions
cfg = [];
cfg.channel = {'all'};
cfg.avgoverfreq = 'no';
cfg.method = 'montecarlo';
cfg.clusterthreshold = 'nonparametric_common';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbours;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 1000;
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.design = design;
cfg.uvar = 1; % row in design indicating subjects, repeated measure
cfg.ivar = 2; % row in design indicating condition for contrast
statTimeFreq = ft_freqstatistics(cfg,bslnTimeFreqAvgStd{:},bslnTimeFreqAvgDev{:});
 
% Plotting stat output over fronto-central electrodes
cfg = [];
cfg.layout = layout;
cfg.channel = {'Fz','FCz','Cz'};
cfg.parameter = 'stat';
cfg.colormap = jet;
cfg.ylim = [8 40];
cfg.marker = 'off';
cfg.style = 'fill';
cfg.comment = 'off';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.colorbar = 'yes';
cfg.zlim = [-5 5];
figure;
ft_singleplotTFR(cfg,statTimeFreq);
set(gca, 'Fontsize',20);
title ('Mean over fronto-central electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('frequency (Hz)');
c = colorbar;
c.LineWidth = 1;
c.FontSize  = 18;
title(c,'t-val');
