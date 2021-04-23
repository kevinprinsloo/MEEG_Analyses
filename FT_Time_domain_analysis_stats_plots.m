%   FT_Time_domain_analysis_stats_plots
%
%   Summary:
%   In FieldTrip, the time-frequency decomposition is computed by the function ft_freqanalysis,
%   specifying one of the supported algorithms in cfg.method. The method mtmconvol implements a STFT,
%   whereas the method wavelet is based on Morlet wavelets (a Morlet wavelet is a complex sinusoid
%   tapered with a Gaussian function).
%
%   Status:
%   Under developement
%
%   Notes:
%
%
%   References:
%      [1] https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/
%      [2] Maris E., Oostenveld R. Nonparametric statistical testing of EEG- and MEG-data.
%          J Neurosci Methods. 2007 Apr 10;
%      [3] Maris E., Schoffelen J.M., Fries P. Nonparametric statistical testing of coherence differences.
%
%   Author(s): Kevin Prinsloo <kevinvago@gmail.com>
%              Cognitive Neurophysiology Laboratory,
%              University of Rochester, NY
%              May 2021; Last Revision: 23-May-2021
%
%   Copyright (c) 2021 Kevin Prinsloo

% Looking at the different trigger codes present in the dataset
cfg = [];
cfg.dataset = 'C:\Users\...\S01\S01_Task.vhdr';
% file path
cfg.trialfun = 'ft_trialfun_general'; % the default trial function is fine
% when trials are defined according
% to specified trigger values
cfg.trialdef.eventtype = '?';
ft_definetrial(cfg);


% Defining trials using trigger values – Standard tones
cfg = [];
cfg.dataset = 'C:\Users\...\S01\S01_Task.set';
% file path. You may need to edit the path
cfg.trialfun = 'ft_trialfun_general'; % the default trial function is fine
% when trials are defined according
% to specified trigger values
cfg.trialdef.eventtype = 'Stimulus';
cfg.trialdef.eventvalue = {'S 103','S 104','S 105','S 106', ...
    'S 107','S 108','S 109','S 113','S 118'}; % standard tones
cfg.trialdef.prestim = 0.25; % in seconds
cfg.trialdef.poststim = 0.75; % in seconds
cfg = ft_definetrial(cfg);


% Removing trials containing discontinuities
trialsToReject = zeros(length(cfg.trl),1);
for eventId = 1:length(cfg.event)
    if isempty(cfg.event(eventId).type)
        tmpSample = cfg.event(eventId).sample;
        for trialInfoId = 1:length(cfg.trl)
            if (tmpSample - cfg.trl(trialInfoId,1)) * ...
                    (tmpSample - cfg.trl(trialInfoId,2)) < 0
                trialsToReject(trialInfoId) = trialInfoId;
            end
        end
    end
end
endtrialsToReject = trialsToReject(trialsToReject ~= 0);
cfg.trl(trialsToReject,:) = [];

% Reading trials and applying baseline correction – Standard tones
cfg.demean = 'yes';
cfg.baselinewindow = [-0.1 0]; % in seconds
dataStd = ft_preprocessing(cfg);

% Plotting a single trial from one channel
figure;
plot(dataStd.time{1}, dataStd.trial{1}(5,:)) % trial 1 from channel 5 (Fz)

% Defining trials using trigger values – Deviant tones
cfg = [];
cfg.dataset = 'C:\Users\...\Task.set';
% file path
cfg.trialfun = 'ft_trialfun_general'; % the default trial function is fine
% when trials are defined according
% to specified trigger values
cfg.trialdef.eventtype = 'Stimulus';
cfg.trialdef.eventvalue = {'S 203','S 204','S 205','S 206','S 207',...
    'S 208','S 209','S 213','S 218'}; % deviant tones
cfg.trialdef.prestim = 0.25; % in seconds
cfg.trialdef.poststim = 0.75; % in seconds
cfg = ft_definetrial(cfg);

% Removing trials containing discontinuities
trialsToReject = zeros(length(cfg.trl),1);
for eventId = 1:length(cfg.event)
    if isempty(cfg.event(eventId).type)
        tmpSample = cfg.event(eventId).sample;
        for trialInfoId = 1:length(cfg.trl)
            if (tmpSample - cfg.trl(trialInfoId,1)) * ...
                    (tmpSample - cfg.trl(trialInfoId,2)) < 0
                trialsToReject(trialInfoId) = trialInfoId;
            end
        end
    end
end

trialsToReject = trialsToReject(trialsToReject ~= 0);
cfg.trl(trialsToReject,:) = [];

% Reading trials and applying baseline correction – Deviant tones
cfg.demean = 'yes';
cfg.baselinewindow = [-0.1 0]; % in seconds
dataDev = ft_preprocessing(cfg);

%% Initializations
dataPath = 'C:\Users...'; % It depends on where
% the wikiDatasets folder
% is stored
fileNameDescription = 'Task.set';

subjects = {'S01' 'S02' 'S03' ...
    'S04' 'S05'}; % the variable subjects contains the list
% of subject-specific directories
nSubjs = numel(subjects);


%% Defining trials (or epochs) using trigger values
for i = 1:nSubjs
    fileName = strcat(subjects{i},fileNameDescription);
    dataSet = fullfile(dataPath,subjects{i},fileName);
    
    % Defining trials using trigger values - Standard tones
    cfg = [];
    cfg.dataset = dataSet; % file path
    cfg.trialfun = 'ft_trialfun_general'; % the default trial function is fine
    % when trials are defined according
    % to specified trigger values
    cfg.trialdef.eventtype = 'Stimulus';
    cfg.trialdef.eventvalue = {'S 103','S 104','S 105','S 106','S 107',...
        'S 108','S 109','S 113','S 118'}; % standard tones
    cfg.trialdef.prestim = 0.25; % in seconds
    cfg.trialdef.poststim = 0.75; % in seconds
    cfg = ft_definetrial(cfg);
    
    % Removing trials containing discontinuities
    trialsToReject = zeros(length(cfg.trl),1);
    for eventId = 1:length(cfg.event)
        if isempty(cfg.event(eventId).type)
            tmpSample = cfg.event(eventId).sample;
            for trialInfoId = 1:length(cfg.trl)
                if (tmpSample - cfg.trl(trialInfoId,1)) * ...
                        (tmpSample - cfg.trl(trialInfoId,2)) < 0
                    trialsToReject(trialInfoId) = trialInfoId;
                end
            end
        end
    end
    
    trialsToReject = trialsToReject(trialsToReject ~= 0);
    cfg.trl(trialsToReject,:) = [];
    
    % Reading trials and applying baseline correction - Standard tones
    cfg.demean = 'yes';
    cfg.baselinewindow = [-0.1 0]; % in seconds
    cfg.outputfile = fullfile(dataPath,subjects{i},'dataStd');
    ft_preprocessing(cfg);
    
    % Defining trials using trigger values - Deviant tones
    cfg = [];
    cfg.dataset = dataSet; % file path
    cfg.trialfun = 'ft_trialfun_general'; % the default trial function is fine
    % when trials are defined according
    % to specified trigger values
    cfg.trialdef.eventtype = 'Stimulus';
    cfg.trialdef.eventvalue = {'S 203','S 204','S 205','S 206','S 207',...
        'S 208','S 209','S 213','S 218'}; % deviant tones
    cfg.trialdef.prestim = 0.25; % in seconds
    cfg.trialdef.poststim = 0.75; % in seconds
    cfg = ft_definetrial(cfg);
    
    % Removing trials containing discontinuities
    trialsToReject = zeros(length(cfg.trl),1);
    for eventId = 1:length(cfg.event)
        if isempty(cfg.event(eventId).type)
            tmpSample = cfg.event(eventId).sample;
            for trialInfoId = 1:length(cfg.trl)
                if (tmpSample - cfg.trl(trialInfoId,1)) * ...
                        (tmpSample - cfg.trl(trialInfoId,2)) < 0
                    trialsToReject(trialInfoId) = trialInfoId;
                end
            end
        end
    end
    
    trialsToReject = trialsToReject(trialsToReject ~= 0);
    cfg.trl(trialsToReject,:) = [];
    
    % Reading trials and applying baseline correction - Deviant tones
    cfg.demean = 'yes';
    cfg.baselinewindow = [-0.1 0]; % in seconds
    cfg.outputfile = fullfile(dataPath,subjects{i},'dataDev');
    ft_preprocessing(cfg);
end

%% Computing ERP for each subject

% Standard tones
for i = 1:nSubjs
    fileName = 'dataStd.mat';
    dataSet = fullfile(dataPath,subjects{i},fileName);
    load(dataSet); % loading the variable 'data'
    
    % Averaging across trials
    cfg = [];
    cfg.channel = {'all'};
    cfg.outputfile = fullfile(dataPath,subjects{i},'erpStd');
    ft_timelockanalysis(cfg,data);
    clear data
end

% Deviant tones
for i = 1:nSubjs
    fileName = 'dataDev.mat';
    dataSet = fullfile(dataPath,subjects{i},fileName);
    load(dataSet); % loading the variable 'data'
    
    % Averaging across trials
    cfg = [];
    cfg.channel = {'all'};
    cfg.outputfile = fullfile(dataPath,subjects{i},'erpDev');
    ft_timelockanalysis(cfg,data);
    clear data
end

%% Computing group average ERP

% Standard tones
% Loading single subject ERPs in a single cell array
allErpStd = cell(1,nSubjs);
for i = 1:nSubjs
    fileName = 'erpStd.mat';
    dataSet = fullfile(dataPath,subjects{i},fileName);
    load(dataSet); % loading the variable 'timelock'
    allErpStd{i} = timelock;
    clear timelock
end

% Computing grand average across subjects
cfg = [];
cfg.latency = [-0.1 0.6];
cfg.keepindividual = 'yes'; % it is sometimes useful to keep individual...
% representations in the output
cfg.outputfile = fullfile(dataPath,'grandAvgErpStd');
ft_timelockgrandaverage(cfg,allErpStd{:});

% Deviant tones
% Loading single subject ERPs in a single cell array
allErpDev = cell(1,nSubjs);
for i = 1:nSubjs
    fileName = 'erpDev.mat';
    dataSet = fullfile(dataPath,subjects{i},fileName);
    load(dataSet); % loading the variable 'timelock'
    allErpDev{i} = timelock;
    clear timelock
end

% Computing grand average across subjects
cfg = [];
cfg.latency = [-0.1 0.6];
cfg.keepindividual = 'yes'; % it is sometimes useful to keep individual...
% representations in the output
cfg.outputfile = fullfile(dataPath,'grandAvgErpDev');
ft_timelockgrandaverage(cfg,allErpDev{:});

%% Defining channel layout

templateLayout = 'EEG128.lay'; % one of the template layouts included
% in FieldTrip
cfg = [];
cfg.layout = which(templateLayout);
layout = ft_prepare_layout(cfg);

figure
ft_plot_lay(layout);
title(templateLayout)

% Increasing layout width and height
layout.width = 0.0615 * ones(length(layout.width),1);
layout.height = 0.0356 * ones(length(layout.height),1);

%% Visualization of the group average ERP

% Plotting a channel layout with the ERP time courses of standard and...
% deviant tones at each electrode
fileName = 'grandAvgErpStd';
dataSet = fullfile(dataPath,fileName);
standard = load(dataSet); % loading the variable 'grandavg'

fileName = 'grandAvgErpDev';
dataSet = fullfile(dataPath,fileName);
deviant = load(dataSet); % loading the variable 'grandavg'

cfg = [];
cfg.layout = layout;
cfg.graphcolor = 'br';
cfg.showcomment = 'no';
cfg.showscale = 'no';
%cfg.showlabels = 'yes';
figure;
ft_multiplotER(cfg,standard.grandavg,deviant.grandavg);

% Plotting mean ERPs over fronto-central electrodes
cfg = [];
%cfg.ylim = [-1 1];
cfg.channel = {'Fz','FCz','Cz'};
cfg.linewidth = 2;
cfg.graphcolor = 'br';
figure;
ft_singleplotER(cfg,standard.grandavg,deviant.grandavg);
legend({'Standard','Deviant'})

% Matlab commands to make the figure intelligible
set(gca,'Fontsize',14);
%title('Mean over central electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('amplitude (\muV)');

% Plotting mean ERP of deviant tones over fronto-central electrodes...
% toghether with standard deviation
cfg = [];
cfg.channel = {'Fz','FCz','Cz'};
cfg.avgoverchan = 'yes'; % averaging over channels
grandavgSel = ft_selectdata(cfg,deviant.grandavg);
x = grandavgSel.time'; % defining x axis
y = mean(squeeze(grandavgSel.individual),1)'; % defining y axis
e = std(squeeze(grandavgSel.individual),1)'; % standard deviation
low = y - e; % lower bound
high = y + e; % upper bound
figure;
hp = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'r');
hold on;
hl = line(x,y);
set(hp,'facecolor',[1 0.8 0.8],'edgecolor','none');
set(hl,'color','r','linewidth',2);
set(gca,'FontSize',12);
title ('Mean and Standard deviation over fronto-central electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('amplitude (\muV)');
legend({'Deviant (StD)','Deviant (Mean)'})

% Plotting topography of the N100 component
cfg = [];
cfg.layout = layout;
cfg.xlim = [0.1 0.12];
cfg.zlim = [-2 2];
cfg.colormap = parula;
cfg.marker = 'off';
cfg.style = 'fill';
cfg.comment = ' ';
cfg.colorbar = 'yes';

figure;
subplot(2,1,1);ft_topoplotER(cfg,standard.grandavg);
title('Standard');
c = colorbar;
c.LineWidth = 1;
c.FontSize = 14;
title(c,'\muV');

subplot(2,1,2);ft_topoplotER(cfg,deviant.grandavg);
title('Deviant');
c = colorbar;
c.LineWidth = 1;
c.FontSize = 14;
title(c,'\muV');
suptitle('N100 topography');

%---------------------
%% FT Cluster Stats
%---------------------

%% Creating a neighbourhood structure to define how spatial clusters...
% are formed
fileName = 'dataStd.mat';
dataSet = fullfile(dataPath,subjects{1},fileName);
load(dataSet); % loading the variable 'data'

templateElec = 'standard_128.elc';
elec = ft_read_sens(which(templateElec));

idx = ismember(elec.label,data.label);
elec.chanpos = elec.chanpos(idx,:);
elec.chantype = elec.chantype(idx);
elec.chanunit = elec.chanunit(idx);
elec.elecpos = elec.elecpos(idx,:);
elec.label = elec.label(idx);

cfg = [];
cfg.method = 'distance'; %'template' or 'distance' or 'triangulation'
cfg.neighbourdist = 80; % maximum distance between neighbouring sensors...
% (only for 'distance')
cfg.feedback = 'yes'; % show a neighbour plot
neighbours = ft_prepare_neighbours(cfg,elec);

%% Statistical analysis of time-domain data at the group level

% Creating a FieldTrip design matrix
design = zeros(1,2*nSubjs); % for within-subjects analysis,...
% the design matrix contains two rows
for i = 1:nSubjs
    design(1,i) = i;
end
for i = 1:nSubjs
    design(1,nSubjs+i) = i;
end
design(2,1:nSubjs) = 1;
design(2,nSubjs+1:2*nSubjs) = 2;

% Computing the difference as a T-statistic and running an inferential test
cfg = [];
cfg.channel = {'all'};
%cfg.latency = [0.05 0.15];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.clusterthreshold = 'nonparametric_common';
cfg.correctm = 'cluster';
cfg.alpha = 0.025; % 0.05/2 because two-sided
cfg.numrandomization = 2000;
cfg.neighbours = neighbours;
cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;
statErp = ft_timelockstatistics(cfg,allErpStd{:},allErpDev{:});

% Plotting stat output
cfg = [];
cfg.layout = layout;
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.graphcolor = 'r';
cfg.showcomment = 'no';
cfg.showscale = 'no';
figure;
ft_multiplotER(cfg,statErp);

% Plotting stat difference ERP over fronto-central electrodes
cfg = [];
cfg.ylim = [-6 16];
cfg.channel = {'Fz','FCz','Cz'};
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.linewidth = 2;
cfg.graphcolor = 'r';
figure;
ft_singleplotER(cfg,statErp);
set(gca,'Fontsize',14);
title ('Mean over fronto-central electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('std vs. dev (t-val)');

% Plotting N100 topography in units of t-value
cfg = [];
cfg.layout = layout;
cfg.xlim = [0.1 0.12];
cfg.parameter = 'stat';
cfg.zlim = [-5 5];
cfg.colormap = parula;
cfg.marker = 'off';
cfg.style = 'fill';
cfg.comment = ' ';
cfg.colorbar = 'yes';
figure;
ft_topoplotER(cfg,statErp);
title('N100 topography')
c = colorbar;
c.LineWidth = 1;
c.FontSize = 14;
title(c,'t-val')

%% Single-subject analysis of time-domain data

subjId = 1; % subject that we will analyze

% Loading the trials with standard tones
fileName = 'dataStd.mat';
dataSet = fullfile(dataPath,subjects{subjId},fileName);
load(dataSet); % loading the variable 'data'

% Averaging across trials
cfg = [];
cfg.channel = {'all'};
cfg.keeptrials = 'yes';
timeLockStd = ft_timelockanalysis(cfg,data);

% Loading the trials with deviant tones
fileName = 'dataDev.mat';
dataSet = fullfile(dataPath,subjects{subjId},fileName);
load(dataSet); % loading the variable 'data'

% Averaging across trials
cfg = [];
cfg.channel = {'all'};
cfg.keeptrials = 'yes';
timeLockDev = ft_timelockanalysis(cfg,data);

% Creating a FieldTrip design matrix
design = zeros(1,size(timeLockStd.trial,1) + size(timeLockDev.trial,1)); % for between-trials analysis, the design matrix contains only one row
design(1,1:size(timeLockStd.trial,1)) = 1;
design(1,(size(timeLockStd.trial,1)+1):(size(timeLockStd.trial,1) + size(timeLockDev.trial,1)))= 2;

% Computing the difference as a T-statistic and running an inferential test
cfg = [];
cfg.channel = {'all'};
%cfg.latency = [0.05 0.15];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesT'; % independent samples T-statistic for a between-trials analysis
cfg.clusterthreshold = 'nonparametric_common';
cfg.correctm = 'cluster';
cfg.alpha = 0.025; % 0.05/2 because two-sided
cfg.numrandomization = 1000;
cfg.neighbours = neighbours;
cfg.design = design;
cfg.ivar = 1;
statErp2 = ft_timelockstatistics(cfg,timeLockStd,timeLockDev);

% Plotting stat output
cfg = [];
cfg.layout = layout;
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.graphcolor = 'r';
cfg.showcomment = 'no';
cfg.showscale = 'no';
figure;
ft_multiplotER(cfg,statErp2);


%% Group analysis of time-domain data (activity vs. baseline)

% Splitting data into activity and baseline
act = cell(1,nSubjs);
bl = cell(1,nSubjs);
for i = 1:nSubjs
    cfg = [];
    cfg.latency = [-0.1 0];
    bl{i} = ft_selectdata(cfg,allErpStd{i});
    cfg.latency = [0.05 0.15];
    act{i} = ft_selectdata(cfg,allErpStd{i});
    bl{i}.time = act{i}.time; % FieldTrip allows for a direct comparison only if the time-axes are equal
end

% Creating a FieldTrip design matrix, conditions: 1 = act, 2 = bl
design = zeros(2,2*nSubjs);
for i = 1:nSubjs
    design(1,i) = i;
end
for i = 1:nSubjs
    design(1,nSubjs+i) = i;
end
design(2,1:nSubjs) = 1;
design(2,nSubjs+1:2*nSubjs) = 2;

% Computing the difference as a T-statistic, and running an inferential test
cfg = [];
cfg.channel = {'all'};
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT'; % dependent samples T-statistics
cfg.clusterthreshold = 'nonparametric_common';
cfg.correctm = 'cluster';
cfg.alpha = 0.025;% .05/2 because two-sided
cfg.numrandomization = 1000;
cfg.neighbours = neighbours;
cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;
statErp3 = ft_timelockstatistics(cfg,act{:},bl{:});

% Plotting stat output
cfg = [];
cfg.layout = layout;
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.graphcolor = 'r';
cfg.showcomment = 'no';
cfg.showscale = 'no';
figure;
ft_multiplotER(cfg,statErp3);
