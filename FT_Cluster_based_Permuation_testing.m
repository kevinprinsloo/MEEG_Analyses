%   FT_Cluster_based_Permuation_testing.
%
%   Summary:
%   Fiedltrip implimentation to perform base cluster-based permuation
%   testing on event related fields. Additional code plotting of results
%
%   Status:
%   Incomplete
%
%   Notes:
%   This code required to data to be in a Fiedltrip format data strucure.
%   If you data is not already in this structure, section 1 allows you to
%   convert EEGLAB data structure to Fieldtrip strcuture
%   M/EEG-data are two-dimensional: channels (the spatial dimension) X time-points (the temporal dimension).
%
%
%   References:
%      [1] https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/
%      [2] Maris E., Oostenveld R. Nonparametric statistical testing of EEG- and MEG-data.
%          J Neurosci Methods. 2007 Apr 10;
%      [3] Maris E., Schoffelen J.M., Fries P. Nonparametric statistical testing of coherence differences.
%          J Neurosci Methods. 2007 Jun 15;163(1):161-75.
%
%   Author(s): Kevin Prinsloo <kevinvago@gmail.com>
%              Cognitive Neurophysiology Laboratory,
%              University of Rochester, NY
%              May 2021; Last Revision: 23-May-2021
%
%   Copyright (c) 2021 Kevin Prinsloo 

%% Initialise Path Variables
close all
clearvars
eeglab
close all
ft_defaults

% Converting between EEGLAB and FieldTrip data structures
% Use either method 1 or 2 to convert EEGLAB to Fiedltrip data structure
EEG = fieldtrip2eeglab(hdr, dat, events); % Method 1
EEG = pop_fileio(hdr, dat, events); % Method 2

% To import EEGLAB datasets in FieldTrip is to read EEGLAB datasets
% using the file File-IO interface either using the standard FILE-IO
% interface or using the low-level EEGLAB reading function of FILE-IO as below.
hdr = ft_read_header( EEGLABFILE );
data = ft_read_data( EEGLABFILE, 'header', hdr );
events = ft_read_event( EEGLABFILE, 'header', hdr );

% In this example we will compare two groups of subjects:
% Healthy controls and Patients

% Timelock fata
cfg = [];
cfg.keepindividual = 'yes';
Gavg_dat_1 = ft_timelockgrandaverage(cfg,  Gavg_1{:});
Gavg_dat_2 = ft_timelockgrandaverage(cfg,  Gavg_2{:});

%---------------------------
%% Run cluster Statistics
%---------------------------

% The configuration settings:
% Initialise Varibels 
channels_number_cephalic = 64; % 64 channel Biosemi layout

% Define EEG channel montage layout - Neighbourhood structure
cfg_neighb        = [];
cfg_neighb.method = 'triangulation';
cfg.layout = ['biosemi',num2str(channels_number_cephalic),'.lay'];
cfg_neighb.senstype = 'eeg';
neighbours        = ft_prepare_neighbours(cfg_neighb, Gavg_dat_1);
cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with
                                 % which other sensors it can form clusters
                                 
% Now run the Permutation Cluster Stats 
cfg=[];
cfg.parameter         = 'individual';
%cfg.parameter         = 'avg';
cfg.method            = 'montecarlo';    % use the Monte Carlo Method to calculate the significance probability
cfg.statistic         = 'indepsamplesT'; % use the independent samples T-statistic as a measure to
                                         % evaluate the effect at the sample level
cfg.neighbours        = neighbours;
cfg.clusteralpha      = 0.05;
cfg.correctm          = 'cluster';       % alpha level of the sample-specific test statistic that
                                         % will be used for thresholding
cfg.clusterstatistic  = 'maxsum';        % test statistic that will be evaluated under the permutation distribution.
cfg.alpha             = 0.05;
cfg.tail              = 0;               % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.correcttail       = 0;
cfg.alpha             = 0.025;           % alpha level of the permutation test
cfg.numrandomization  = 2000;            % number of draws from the permutation distribution
cfg.minnbchan         = 2;               % minimum number of neighborhood channels that is
                                         % required for a selected sample to be included
                                         % in the clustering algorithm (default=0).
isub = size(Gavg_dat_1.individual,1);
cfg.design(1,1:2*isub)=[ones(1,isub) 2*ones(1,isub)];
cfg.design(2,1:2*isub)=[1:isub 1:isub];
cfg.ivar=1; % the 1st row in cfg.design contains the independent variable
cfg.uvar=2; % the 2nd row in cfg.design contains the subject number

% All the restuls are stored in 'stat'
stat=ft_timelockstatistics(cfg,Gavg_dat_1, Gavg_dat_2);


%% Now have a quick look at the results of the analyses [channel x time]
figure;
imagesc(stat,'AlphaData', stat.mask, 'XData', [-0.1 0.8]); % Latency of epoched data -> [-0.1 0.8]
xticklabels({' ','0',' ','0.2',' ','0.4',' ','0.6',' ','0.8'})
xlabel('time (s)');
set(gca, 'Fontsize',14, 'FontName',  'HelvLight', 'box','off','LineWidth',1);
set(gca,'color','none');
colorbar;
ylabel(h, 'Perm \itt\rm-value')
colormap parula

%% Plot ERPs of the two subjects groups with subplot containing stats

% Choose two colours for you subject groups
clear c1 cl c3
[cb] =  cbrewer('seq','PuBuGn',16,'pchip');
c3(1,:) = cb(13,:); % 11
[cb] = cbrewer('seq','YlOrRd',16,'pchip');
cl(1,:) = cb(16,:); % 

% Select channel to plot - each channel will be a new plot
sel_chans = {'CPz','FCz','Cz'};

% Plot ERPS with subplot with stats
% Loop over our channels and make a new plot for each channel 
for ichan = 1:length(sel_chans)
    
    figure
    hhh = axes('Position', [0.12 0.38 0.85 0.6]);  % [left bottom width height]
    hold on
            
    cfg = [];
    cfg.channel = sel_chans{ichan};
    chan_num_sel = find(strcmp(chan_name_plot,label_names_list));   
    cfg.avgoverchan = 'yes'; % averaging over channels if you wat but here we are looking only at 1 channel
    
    % >> Now deal with subject group 1
    grandavgSel_1 = ft_selectdata(cfg,grandavgSel_11);
    x = grandavgSel_1.time'; % defining x axis
    y = nanmean(squeeze(grandavgSel_1.individual),1)'; % defining y axis    
    %e = std(squeeze(grandavgSel_1.individual),1)'; % (STD) standard deviation
    e = std(squeeze(grandavgSel_1.individual),1)'/sqrt(length(goodsubs)); % (SEM) standard error of the mean
    low = y - e; % lower bound
    high = y + e; % upper bound
    hp = patch([x; x(end:-1:1); x(1)], [low; high(end:-1:1); low(1)], 'r');
    hold on;
    h1 = line(x,y);
    set(hp,'facecolor',c3(1,:),'edgecolor','none','FaceAlpha', 0.2);
    set(h1,'color',c3(1,:),'linewidth',2.4);
    set(gca,'FontSize',12);
    alpha(0.2)
    hold on
    
    % >> Now deal with subject group 2
    grandavgSel_2 = ft_selectdata(cfg,grandavgSel_22);
    x2 = grandavgSel_2.time'; % defining x axis
    y2 = mean(squeeze(grandavgSel_2.individual),1)'; % defining y axis
    %e2 = std(squeeze(grandavgSel_2.individual),1)'; % (STD) standard deviation
    e2 = std(squeeze(grandavgSel_2.individual),1)'/sqrt(length(goodsubs)); % (SEM) standard error of the mean
    low = y2 - e2; % lower bound
    high = y2 + e2; % upper bound
    hp2 = patch([x2; x2(end:-1:1); x2(1)], [low; high(end:-1:1); low(1)], 'b');
    hold on;
    h2 = line(x2,y2);
    set(hp2,'facecolor',cl(1,:),'edgecolor','none','FaceAlpha', 0.8);
    set(h2,'color',cl(1,:),'linewidth',2.4);
    set(gca,'FontSize',16);
    alpha(0.2)
        
    % deal with some plotting extras
    hold on
    t5 = xline(0,'linestyle','--'); t5.Color = [.4 .4 .4];
    t5 = yline(0,'linestyle','--'); t5.Color = [.4 .4 .4];    
    set(gca,'box','on');
    xlabel(' ');
    ylabel('amplitude ({\ituV})','Fontsize',16,'FontName', 'HelvLight'); % , 'Interpreter', 'none'
    
    % Set limits
    xlim([-0.1 0.8]) % time axis
    ylim([-7 3.5])   % Amplitude voltage axis

    % Set legend info   
    leg = legend([hp hp2],{'Group 1','Group 2'});
    set(leg,'interpreter','none');    
    set(leg,'color','none','linewidth',1,'location','NorthEast')
    legend('BoxOff');
    set(gca, 'box','off','LineWidth',1);
    set(gca,'color','none');
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    ax = gca;
    ax.XColor = 'none';
    axPos = hhh.Position;
    ax2 = axes('position', (axPos .* [1 1 1 1e-3]) + [0 -0.05 0 0], 'color', 'none','Fontsize',16, 'FontName', 'HelvLight');
    ax2.XColor = 'none';
    ax2.XTickLabel = [];
    ax.XColor = 'none';
    
    % >> Now Plot subplot with stats    
    hh1 = axes('Position', [0.12 0.24 0.85 0.1]); % [left bottom width height]
    t5 = xline(0,'linestyle','--'); t5.Color = [.4 .4 .4];
    t5 = yline(0,'linestyle','--'); t5.Color = [.4 .4 .4];
    set(gca,'xlim',[-0.1 0.8]);
    
    % Deal with stats line for ourt selected channel
    maskLine = stat.mask(chan_num_sel,:);
    thinkLineTM = stat.stat(chan_num_sel,:);
    
    % Create figure and line thickness - if the line is thick is mean those
    % times are significant - think line is non-significant
    hold on
    thick = 1;
    t = grandavgSel_1.time;
    % Loop through logicals
    for k = 2:length(maskLine)
        if maskLine(k) == true
            plot(t((k-1):k),thinkLineTM((k-1):k),'LineWidth',4,'Color',[0.6 0.6 0.9])
            thick = thick + 1;
        else
            plot(t((k-1):k),thinkLineTM((k-1):k),'LineWidth',1.4,'Color',[0.5 0.5 0.5])
            thick = thick - 1;
        end        
        % Prevent line from dissapearing
        if thick <= 1
            thick = 1;
        end
    end
    ax3 = gca;
    ax3.XColor = 'none';
    axPos = hh1.Position;
    ax2 = axes('position', (axPos .* [1 1 1 1e-3]) + [0 -0.05 0 0], 'color', 'none','Fontsize',12, 'FontName', 'HelvLight');
    ax3.Color = 'none';
    ax3.XLabel.String = 'Time ({\itms})';
    ax.Color = 'none';
    set(gca, 'Fontsize',16, 'FontName',  'HelvLight', 'box','off','LineWidth',1);
    set(gca,'color','none');
    
    % Now label the pl;ot
    ylabel('\itt','Fontsize',14,'FontName', 'HelvLight')
    xlabel('Time ({\itms})','Fontsize',14,'FontName', 'HelvLight')
    ax3.YLabel.String = '\itt\rm-value';    
    hh1.YAxis.FontSize = 10;
    ylabel('\itt\rm-value','Fontsize',12,'FontName', 'HelvLight')
    set(gca,'xlim',[-100 800]);
    
end

%------------------
%% Topo Plots
%------------------

% Next we can plot a few topoplots with significant clusters
% Another way of looking at the data but now in topos averaged over time

% Take the difference of the averages (ERPs across out two groups to get difference waves) using ft_math
cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffectFICvsFC = ft_math(cfg, Gavg_1, Gavg_2);
pos_cluster_pvals = [stat.posclusters(:).prob];

% Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
% cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
% respectively
pos_clust = find(pos_cluster_pvals < 0.5);
pos       = ismember(stat.posclusterslabelmat, pos_clust);

% and now for the negative clusters...
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_clust         = find(neg_cluster_pvals < 0.5);
neg               = ismember(stat.negclusterslabelmat, neg_clust);

timestep      = 0.05; % timestep between time windows for each subplot (in seconds)
sampling_rate = 256; % Data has a temporal resolution of 300 Hz
sample_count  = length(stat.time);
% number of temporal samples in the statistics object
j = [0:timestep:1]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in M/EEG samples

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffectFICvsFC.label, stat.label);

figure
for k = 1:17
   subplot(4,5,k);
   cfg = [];
   cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
   %cfg.zlim = [-2.5e-13 2.5e-13];
   % If a channel is in a to-be-plotted cluster, then
   % the element of pos_int with an index equal to that channel
   % number will be set to 1 (otherwise 0).

   % Next, check which channels are in the clusters over the
   % entire time interval of interest.
   pos_int = zeros(numel(raweffectFICvsFC.label),1);
   neg_int = zeros(numel(raweffectFICvsFC.label),1);
   pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
   neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);

   cfg.highlight   = 'on';
   cfg.highlightsymbol = '.';
   cfg.highlightsize = 18;
   cfg.highlightcolor=[1 0 0];

   % Get the index of the to-be-highlighted channel
   cfg.highlightchannel = find(pos_int | neg_int);
   cfg.comment     = 'xlim';
   cfg.commentpos  = 'title';
   cfg.layout = ['biosemi',num2str(channels_number_cephalic),'.lay'];
   cfg.interactive = 'no';
   ft_topoplotER(cfg, raweffectFICvsFC);
end

% But now only plot two time points of interest
for k =  [6 10]
    
    figure
    cfg = [];
    cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
    %cfg.zlim = [-2.5e-13 2.5e-13];
    % If a channel is in a to-be-plotted cluster, then
    % the element of pos_int with an index equal to that channel
    % number will be set to 1 (otherwise 0).
    
    % Next, check which channels are in the clusters over the
    % entire time interval of interest.
    pos_int = zeros(numel(raweffectFICvsFC.label),1);
    neg_int = zeros(numel(raweffectFICvsFC.label),1);
    pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
    neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);
    
    cfg.highlight   = 'on';
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 26;
    cfg.highlightcolor=[1 1 1];
    
    % Get the index of the to-be-highlighted channel
    cfg.highlightchannel = find(pos_int | neg_int);
    cfg.layout = ['biosemi',num2str(channels_number_cephalic),'.lay'];
    cfg.interactive = 'no';
    cfg.comment = 'no';
    cfg.marker = 'off';
    cfg.zlim = [-60 60];
     cfg.gridscale = 128;
    ft_topoplotER(cfg, raweffectFICvsFC);
    colormap parula
    set(gca, 'Fontsize',24, 'FontName',  'HelvLight', 'box','off','LineWidth',1);
    %colorbar;
end

%------------------------------------------------
%% Box plots of single subjects across groups
%------------------------------------------------

% Now you can plot two Box plots comparing between the groups
% of EEG averaged over a time window

cfg = [];
cfg.channel =  {'CPz'};
cfg.latency = [0.2 0.3]; % Average data over 200-300 ms
cfg.avgovertime = 'yes';
grandavgSel_1 = ft_selectdata(cfg,Gavg_dat_1);
grandavgSel_2 = ft_selectdata(cfg,Gavg_dat_1);
data_1 = squeeze(grandavgSel_1.individual);
data_2 = squeeze(grandavgSel_2.individual);

% Now plot data
figure
ha = tight_subplot(1,2,[.1 0.01],[.12 .08],[0.15 0.01]);
hold on
ha(1).XLabel.String = 'N1 Amp'; ha(1).YColor = [0 0 0];
set(ha(1),'color','none','ycolor','none');
set(ha(1), 'Fontsize',12, 'FontName',  'HelvLight', 'box','off','LineWidth',1);
% Defien the space between plots
box_dodge_amount_2 = .26;
% Plot the rest
hold on
axes(ha(1));
h1 = raincloud_plot('X', data_1 , 'box_on', 1, 'color', c3(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .1, 'dot_dodge_amount', .1,...
    'box_col_match', 0,'line_width', 1.2, 'sig_value', sig_value(1));
h2= raincloud_plot('X',data_2, 'box_on', 1, 'color', cl(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', box_dodge_amount_2, 'dot_dodge_amount', box_dodge_amount_2, 'box_col_match', 0,'line_width', 1.2, 'sig_value', sig_value(2));
box off
set(gca,'color','none');
set(gca, 'Fontsize',16, 'FontName',  'HelvLight', 'box','off','LineWidth',1);
view([90 -90]);
xlabel('amplitude (\ituV)');
% Plot legend
leg = legend([h2{1} h1{1}],{'Group_1','Group_2'});
set(leg,'color','none','linewidth',1,'location','NorthEast')
set(leg,'Interpreter','none')
legend('BoxOff');
% Finish of plotting
ha = tight_subplot(1,2,[.1 0.01],[.12 .08],[0.15 0.01]);
hold on
ha(1).XLabel.String = 'N200'; ha(1).YColor = [0 0 0];
set(ha(1),'color','none','ycolor','none');
set(ha(1), 'Fontsize',16, 'FontName',  'HelvLight', 'box','off','LineWidth',1);


%------------------
%% Cluster Plots
%------------------

AL64 = {'Fp1' 'AF7' 'AF3'}; % 1-3
AC64 = {'Fpz' 'Afz'};       % 4-5 
AR64 = {'Fp2' 'AF8' 'AF4'}; % 6-8
FL64 = {'F7' 'FT7' 'F5' 'FC5' 'F3' 'FC3'}; % 9-14
FC64 = {'F1' 'FC1' 'Fz' 'FCz' 'F2' 'FC2'}; % 15-20
FR64 = {'F4' 'FC4' 'F6' 'FC6' 'F8' 'FT8'}; % 21-26
CL64 = {'T7' 'C5' 'C3'}; % 27-29
CC64 = {'C1' 'Cz' 'C2'}; % 30-32
CR64 = {'C4' 'C6' 'T8'}; % 33-35
PL64 = {'TP7' 'P7' 'CP5' 'P5' 'CP3' 'P3'}; % 36-41
PC64 = {'CP1' 'P1' 'CPz' 'Pz' 'CP2' 'P2'}; % 42-47
PR64 = {'CP4' 'P4' 'CP6' 'P6' 'TP8' 'P8'}; % 48-53
OL64 = {'P9' 'PO7' 'PO3' 'O1'};  % 54-57
OC64 = {'POz' 'Oz' 'Iz'};        % 58-60
OR64 = {'O2' 'PO4' 'PO8' 'P10'}; % 61-64

chnLst = {'Fp1' 'AF7' 'AF3' 'Fpz' 'AFz' 'Fp2' 'AF8' 'AF4' 'F7' 'FT7' 'F5' 'FC5' 'F3' 'FC3' 'F1' 'FC1' 'Fz' 'FCz' 'F2' 'FC2' ,...
    'F4' 'FC4' 'F6' 'FC6' 'F8' 'FT8' 'T7' 'C5' 'C3' 'C1' 'Cz' 'C2' 'C4' 'C6' 'T8' 'TP7' 'P7' 'CP5' 'P5' 'CP3' 'P3',...
    'CP1' 'P1' 'CPz' 'Pz' 'CP2' 'P2' 'CP4' 'P4' 'CP6' 'P6' 'TP8' 'P8' 'P9' 'PO7' 'PO3' 'O1' 'POz' 'Oz' 'Iz' 'O2' 'PO4' 'PO8' 'P10'};

clear chan_num_sel
for k = 1:64
    chan_num_sel(k) = find(strcmp(chnLst{k},label_names_list));
end

clear smask pval statT LogicalMask
cnt = 1;
for k = 1:64 
    chIdx = chan_num_sel(k);
    statT(cnt,:) = stat.stat(chIdx,:);
    LogicalMask(cnt,:) = stat.mask(chIdx,:);
    cnt=cnt+1;
end

figure;
h = imagesc(statT,'AlphaData', LogicalMask, 'XData', [-0.1 0.8]);
xticklabels({' ','0',' ','0.2',' ','0.4',' ','0.6',' ','0.8'})
xlabel('time (s)');
set(gca, 'Fontsize',14, 'FontName',  'HelvLight', 'box','off','LineWidth',1);
set(gca,'color','none');
h = colorbar;
ylabel(h, 'Perm \itt\rm-value')
colormap inferno
 ax = gca;
 spatialdivs = [anterofrontalreg(1) frontalreg(1) centralreg(1) parietalreg(1) occreg(1)];
 ax.YTick = spatialdivs;
 ax.YTickLabel = {'AnteroFrontal' 'Frontal' 'Central' 'Parietal' 'Occipital'};
t5 = xline(37,'linestyle','--'); t5.Color = [.4 .4 .4];


[cb1] = cbrewer('seq','YlOrRd',16,'pchip');
ALcolor = cb1(16,:);
ACcolor = cb1(13,:); 
ARcolor = cb1(10,:);

[cb2] = cbrewer('seq','PuBu',16,'pchip');
FLcolor = cb2(16,:); 
FCcolor = cb2(13,:); 
FRcolor = cb2(10,:);

[cb3] = cbrewer('seq','RdPu',16,'pchip');
CLcolor = cb3(16,:);
CCcolor = cb3(13,:); 
CRcolor = cb3(10,:); 

[cb4] = cbrewer('seq','Greens',16,'pchip');
PLcolor = cb4(16,:); 
PCcolor = cb4(13,:); 
PRcolor = cb4(10,:); 

[cb5] = cbrewer('seq','Purples',16,'pchip');
OLcolor = cb5(16,:); 
OCcolor = cb5(13,:); 
ORcolor = cb5(10,:) ;

finalorder = [AL64 AC64 AR64 FL64 FC64 FR64 CL64 CC64 CR64 PL64 PC64 PR64 OL64 OC64 OR64];

ALreg = [0.5 find(strcmp(AL64{end},finalorder))+0.5];
ACreg = [ALreg(2) find(strcmp(AC64{end},finalorder))+0.5];
ARreg = [ACreg(2) find(strcmp(AR64{end},finalorder))+0.5];

FLreg = [ARreg(2) find(strcmp(FL64{end},finalorder))+0.5];
FCreg = [FLreg(2) find(strcmp(FC64{end},finalorder))+0.5];
FRreg = [FCreg(2) find(strcmp(FR64{end},finalorder))+0.5];

CLreg = [FRreg(2) find(strcmp(CL64{end},finalorder))+0.5];
CCreg = [CLreg(2) find(strcmp(CC64{end},finalorder))+0.5];
CRreg = [CCreg(2) find(strcmp(CR64{end},finalorder))+0.5];

PLreg = [CRreg(2) find(strcmp(PL64{end},finalorder))+0.5];
PCreg = [PLreg(2) find(strcmp(PC64{end},finalorder))+0.5];
PRreg = [PCreg(2) find(strcmp(PR64{end},finalorder))+0.5];

OLreg = [PRreg(2) find(strcmp(OL64{end},finalorder))+0.5];
OCreg = [OLreg(2) find(strcmp(OC64{end},finalorder))+0.5];
ORreg = [OCreg(2) find(strcmp(OR64{end},finalorder))+0.5];

anterofrontalreg = [0.5 find(strcmp('AF4',finalorder))+0.5];
frontalreg = [anterofrontalreg(2) find(strcmp('FT8',finalorder))+0.5];
centralreg = [frontalreg(2) find(strcmp('T8',finalorder))+0.5];
parietalreg = [centralreg(2) find(strcmp('P8',finalorder))+0.5];
occreg = [parietalreg(2) length(finalorder)+0.5];

a1 = -0.075; %% ++++++++++++ SCALE FIXED TO INCLUDE BASELINE
b1 = -0.075;
d1 = -0.12;
c1 = -0.12;

patch([c1 a1 b1 d1],[ALreg(1) ALreg(1) ALreg(2) ALreg(2)],ALcolor,'FaceAlpha',1,'EdgeColor','none')
patch([c1 a1 b1 d1],[ACreg(1) ACreg(1) ACreg(2) ACreg(2)],ACcolor,'FaceAlpha',1,'EdgeColor','none')
patch([c1 a1 b1 d1],[ARreg(1) ARreg(1) ARreg(2) ARreg(2)],ARcolor,'FaceAlpha',1,'EdgeColor','none')

patch([c1 a1 b1 d1],[FLreg(1) FLreg(1) FLreg(2) FLreg(2)],FLcolor,'FaceAlpha',1,'EdgeColor','none')
patch([c1 a1 b1 d1],[FCreg(1) FCreg(1) FCreg(2) FCreg(2)],FCcolor,'FaceAlpha',1,'EdgeColor','none')
patch([c1 a1 b1 d1],[FRreg(1) FRreg(1) FRreg(2) FRreg(2)],FRcolor,'FaceAlpha',1,'EdgeColor','none')

patch([c1 a1 b1 d1],[CLreg(1) CLreg(1) CLreg(2) CLreg(2)],CLcolor,'FaceAlpha',1,'EdgeColor','none')
patch([c1 a1 b1 d1],[CCreg(1) CCreg(1) CCreg(2) CCreg(2)],CCcolor,'FaceAlpha',1,'EdgeColor','none')
patch([c1 a1 b1 d1],[CRreg(1) CRreg(1) CRreg(2) CRreg(2)],CRcolor,'FaceAlpha',1,'EdgeColor','none')

patch([c1 a1 b1 d1],[PLreg(1) PLreg(1) PLreg(2) PLreg(2)],PLcolor,'FaceAlpha',1,'EdgeColor','none')
patch([c1 a1 b1 d1],[PCreg(1) PCreg(1) PCreg(2) PCreg(2)],PCcolor,'FaceAlpha',1,'EdgeColor','none')
patch([c1 a1 b1 d1],[PRreg(1) PRreg(1) PRreg(2) PRreg(2)],PRcolor,'FaceAlpha',1,'EdgeColor','none')

patch([c1 a1 b1 d1],[OLreg(1) OLreg(1) OLreg(2) OLreg(2)],OLcolor,'FaceAlpha',1,'EdgeColor','none')
patch([c1 a1 b1 d1],[OCreg(1) OCreg(1) OCreg(2) OCreg(2)],OCcolor,'FaceAlpha',1,'EdgeColor','none')
patch([c1 a1 b1 d1],[ORreg(1) ORreg(1) ORreg(2) ORreg(2)],ORcolor,'FaceAlpha',1,'EdgeColor','none')

patch([0 512 512 0],[anterofrontalreg(1) anterofrontalreg(1) anterofrontalreg(2) anterofrontalreg(2)],[0.85 0.85 0.85],'FaceAlpha',0.25,'EdgeColor','none')
patch([0 512 512 0],[frontalreg(1) frontalreg(1) frontalreg(2) frontalreg(2)],[0.95 0.95 0.95],'FaceAlpha',0.25,'EdgeColor','none') 
patch([0 512 512 0],[centralreg(1) centralreg(1) centralreg(2) centralreg(2)],[0.85 0.85 0.85],'FaceAlpha',0.25,'EdgeColor','none') 
patch([0 512 512 0],[parietalreg(1) parietalreg(1) parietalreg(2) parietalreg(2)],[0.95 0.95 0.95],'FaceAlpha',0.25,'EdgeColor','none') 
patch([0 512 512 0],[occreg(1) occreg(1) occreg(2) occreg(2)],[0.85 0.85 0.85],'FaceAlpha',0.25,'EdgeColor','none') %occipital

hold on
h = imagesc(statT,'AlphaData', LogicalMask, 'XData', [-0.1 0.8]);
xticklabels({' ','0',' ','200',' ','400',' ','600',' ','800'})
xlabel('Time ({\itms})','Fontsize',14,'FontName', 'HelvLight')
set(gca, 'Fontsize',14, 'FontName',  'HelvLight', 'box','off','LineWidth',1);
set(gca,'color','none');
h = colorbar;
ylabel(h, 'Perm \itt\rm-value')
colormap inferno
caxis([-5 5])
ax = gca;
spatialdivs = [anterofrontalreg(1) frontalreg(1) centralreg(1) parietalreg(1) occreg(1)];
ax.YTick = spatialdivs;
ax.YTickLabel = {'AnteroFrontal' 'Frontal' 'Central' 'Parietal' 'Occipital'};
t5 = xline(37,'linestyle','--'); t5.Color = [.4 .4 .4];
