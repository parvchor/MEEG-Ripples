clc
clear 
close all

addpath(genpath('space/seh10/5/halgdev/projects/swilson/training/code/RippleDetection_clean'))

addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
% addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple'))

addpath(genpath('/home/jgarret/mat')) %mask2bounds
addpath(genpath('/home/jgarret/ArtifactDetection'))

addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/eeglab14_1_2b/'))
rmpath(genpath('/space/seh8/1/halgdev/projects/cdickey/eeglab14_1_2b/functions/octavefunc'))
 
%% Paths / Inputs

lfpDirec = '/space/seh10/5/halgdev/projects/swilson/training/data/lfp_data';
subj_list_full = {'SUBJ1'};
         
sleepfiles_set = {1,1,1,1,1,1,1,1,1,1,1,1};
subjID = 1:length(subj_list_full);
winSec = 2.000; % full window size (ripple peak +- winSec)
rippleWindowSec = 0.100; %0.100;-human % zoomed-in window (ripple peak +- rippleWindowSec)
% those two should match the values previously set in RippleSelection.m

species = 'human'; % rodent, human
recordingState = 'wake'; %sleep or wake
tag = recordingState;
locations = {'ALL'}; %NC, HC, TH or AMY
modifier = '';

exportFolder = '/space/seh10/5/halgdev/projects/swilson/training/data/ripples/Figures';
if ~isfolder(exportFolder); mkdir(exportFolder); end   

preProcExportFolder = '/space/seh10/5/halgdev/projects/swilson/training/data/preproc';
              
rippleStatsFolder = '/space/seh10/5/halgdev/projects/swilson/training/data/rippleout';  

%% Load rippleStats and compute bandpasses
for isbj = 1:length(subjID)
    subj = subjID(isbj);
    subject = subj_list_full{subj};

    for l = 1:length(locations)
        location = locations{l};
        load(fullfile(rippleStatsFolder,[subject,'_ripple_stats_',recordingState,'_',location,'_',modifier,'.mat']));

        data = [];
        switch recordingState
            case 'sleep'
                sleepfiles = sleepfiles_set{subj};


                for si = 1:length(sleepfiles)
                    s = sleepfiles(si);
                    if contains(subject,'CC')
                        if strcmp(location,'NC')
                                            
                            load(fullfile(preProcExportFolder,[subject,'_',num2str(s),'sleep','_polCheck.mat']), 'polCheck')

                            LFPdir = sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/segment_indices/',subject);
                            temp = load(fullfile(LFPdir,sprintf('sleep_%i_N23_seg_NC.mat', s)),'sleepdata','chan_labels');
                        elseif strcmp(location,'HC')
                            load(fullfile(preProcExportFolder,[subject,'_',num2str(s),'sleep_',location,'_polCheck.mat']), 'polCheck')

                            LFPdir = sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/noWaveletReject/HC',subject);
                            temp = load(fullfile(LFPdir,sprintf('sleep_sleep%i_N23.mat', s)),'sleepdata','chan_labels');
                        end
                    else
                        load(fullfile(preProcExportFolder,[subject,'_',num2str(s),'_sleep_',location,'_',modifier,'_polCheck.mat']), 'polCheck')
                
                        temp = load(fullfile(lfpDirec,sprintf('%s', subject)),'data','channel');
                        if ~isa(temp.data, 'double')
                            temp.data = double(temp.data);
                        end

                        if size(temp.data,1) > size(temp.data,2); temp.data = temp.data'; end

                        if ~exist('temp.sleepdata','var')
%                             temp.sleepdata = temp.datNREM';
                            temp.sleepdata = temp.data;
                            temp.data = [];
                        end
                    end


                    if isempty(data)
                        data = temp.sleepdata .* polCheck;
                    else
                        data = [data, temp.sleepdata .* polCheck];
                    end
                end

            case 'wake'
                if contains(subject,'CC')
                    LFPdir = sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/waking/',subject);           

                    if strcmp(location,'NC')
                        load(fullfile(LFPdir,'NC','waking_1.mat'),'data');
                    elseif strcmp(location,'HC')
                        load(fullfile(LFPdir,'HC','waking_1.mat'),'data');
                    end
                else
                    LFPdir = lfpDirec;
                    temp = load(fullfile(LFPdir, sprintf('%s', subject)), 'data', 'label');
                    if ~exist('temp.sleepdata','var')
                        temp.sleepdata = temp.data;
%                         temp.sleepdata = temp.data;
                        temp.data = [];
                    end
                end

                data = temp.sleepdata';

        end
    %     
    %     if contains(subject,'CC') && strcmp(location, 'NC')
    %         load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/segment_indices/sleep_%01d_N23_seg_NC_polchecked_newRipp.mat', subject,  sleep_ID))
    %         load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/noWaveletReject/segment_indices/sleep_%01d_N23_polchecked_newRipp.mat', subject,  sleep_ID))
    %     else
    %         load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/%s/%s/%s_%s',subject,recordingState,location,recordingState, '1_polchecked_newRipp.mat'))
    %         load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/%s/%s/%s_%s',subject,recordingState,location,recordingState, '1.mat'));
    %     end
        
        try 
            load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/%s/nan_edge_mask.mat', subject, recordingState),'nan_edge_mask') 
        catch 
            warning('could not load epoch edge mask');
            nan_edge_mask = false(1,size(data,2));
        end

        if isfield(rippleStats,'ephys_ch')
            data = data(rippleStats.ephys_ch,:);
        end

    %     polCheck = checkPolarity(data,rippleStats.fs); %check polarity by comparing slow oscillations w HG
    %     data = polCheck .* data;

        win = round(winSec * rippleStats.fs);
        rippleWindow = round(rippleWindowSec * rippleStats.fs);

        chan_labels = rippleStats.chanLabels;
        rippleband = rippleStats.RB;

        for ch = 1:size(data,1)

            [rippleSleepSet,~] = AnalyzeRipple(species, data, rippleband, win, rippleWindow, 0, 0, nan_edge_mask, 0, ...
                                                           ch, rippleStats.locs{ch}, [], rippleStats.fs, chan_labels, 0,0);
            LFPbandpasses = [];                                           
            fn = fieldnames(rippleSleepSet);
            %      rippleAllSleepSet = rippleSleepSet;
            for f = 1:numel(fn)
                LFPbandpasses(1).(fn{f}){ch} = [];
            end

            fn = fieldnames(rippleSleepSet);
            if numel(rippleSleepSet.rippleAmp) > 5
                for f = 1:numel(fn)
                    fDat = rippleSleepSet.(fn{f});
                    DIM = size(fDat,1);
                    if DIM > 1 || size(fDat,2) == 4001
                        LFPbandpasses.(fn{f}){ch} = [LFPbandpasses.(fn{f}){ch}; fDat];
                    elseif DIM == 1
                        LFPbandpasses.(fn{f}){ch} = [LFPbandpasses.(fn{f}){ch}, fDat];
                    end

                end
            end

            LFPbandpasses.recordingLength = size(data,2);

            %% Plotting
            fprintf(sprintf('exporting plots for channel %s\n', chan_labels(ch)))

            summaryFolder = fullfile(exportFolder, subject, location, ['ChannelSummary_', tag, '_',modifier]);
            singleTrialFolder = fullfile(exportFolder, subject, location, ['SingleTrials_',tag,'_', modifier], num2str(chan_labels(ch)));
            TFexportFolder = fullfile(exportFolder, subject, location, ['TF_', tag,'_',modifier]);

            if ~isfolder(fullfile(summaryFolder,'PNG')); mkdir(fullfile(summaryFolder, 'PNG')); end
            if ~isfolder(fullfile(summaryFolder,'PDF')); mkdir(fullfile(summaryFolder, 'PDF')); end
            if ~isfolder(fullfile(summaryFolder,'Fig')); mkdir(fullfile(summaryFolder, 'Fig')); end

            if ~isfolder(fullfile(singleTrialFolder,'PNG')); mkdir(fullfile(singleTrialFolder, 'PNG')); end
            if ~isfolder(fullfile(singleTrialFolder,'PDF')); mkdir(fullfile(singleTrialFolder, 'PDF')); end
            if ~isfolder(fullfile(singleTrialFolder,'Fig')); mkdir(fullfile(singleTrialFolder, 'Fig')); end

            if ~isfolder(fullfile(TFexportFolder,'PNG')); mkdir(fullfile(TFexportFolder, 'PNG')); end
            if ~isfolder(fullfile(TFexportFolder,'PDF')); mkdir(fullfile(TFexportFolder, 'PDF')); end
            if ~isfolder(fullfile(TFexportFolder,'Fig')); mkdir(fullfile(TFexportFolder, 'Fig')); end
 

            if numel(LFPbandpasses.rippleAmp{ch}) >= 20
               meanLFP = mean(LFPbandpasses.raw{ch});

                for rip_num = 1:size(LFPbandpasses.raw{ch},1)
                   rippleCOV = cov(LFPbandpasses.raw{ch}(rip_num,:), meanLFP);
                   LFPbandpasses.cov{ch}(rip_num) = rippleCOV(1,2);     
                end

                hsummary = chanOverviewPlot(species, LFPbandpasses, rippleStats, ch, win, rippleStats.fs, subject, chan_labels);
                savepdf(hsummary, sprintf('%s/PDF/%s_%s_summary.pdf', summaryFolder, subject, num2str(chan_labels(ch))))
                saveas(hsummary,  sprintf('%s/PNG/%s_%s_summary.png', summaryFolder, subject, num2str(chan_labels(ch))))
                saveas(hsummary,  sprintf('%s/Fig/%s_%s_summary.fig', summaryFolder, subject, num2str(chan_labels(ch))))

                hTF = chanTFplot(LFPbandpasses, ch, win, rippleStats.fs, subject, chan_labels,'TF');
                savepdf(hTF, sprintf('%s/PDF/%s_%s_TimeFreq.pdf', TFexportFolder, subject, num2str(chan_labels(ch))))
                saveas(hTF,  sprintf('%s/PNG/%s_%s_TimeFreq.png', TFexportFolder, subject, num2str(chan_labels(ch))))
                saveas(hTF,  sprintf('%s/Fig/%s_%s_TimeFreq.fig', TFexportFolder, subject, num2str(chan_labels(ch))))

                hTF = chanTFplot(LFPbandpasses, ch, win, rippleStats.fs, subject, chan_labels,'LFP');
                savepdf(hTF, sprintf('%s/PDF/%s_%s_MeanLFP.pdf', TFexportFolder, subject, num2str(chan_labels(ch))))
                saveas(hTF,  sprintf('%s/PNG/%s_%s_MeanLFP.png', TFexportFolder, subject, num2str(chan_labels(ch))))
                saveas(hTF,  sprintf('%s/Fig/%s_%s_MeanLFP.fig', TFexportFolder, subject, num2str(chan_labels(ch))))

              %  Plot by rippleBand amp 
%                 [~,I] = sort(LFPbandpasses.rippleAmp{ch}, 'ascend'); % plot the worst ripples (those with lowest amplitudes)
                [~,I] = sort(LFPbandpasses.rippleAmp{ch}, 'descend'); 
                for f = 1:4
                    h3 = figure('Position', [3 30 1908 885], 'Visible', 'off');
                    for i = 1:4

                        rip_num = I(4*(f-1) + i); 
%                         rip_num = randi(length(I),1); % random ripple plotting
                        rippleFullPlot(species, LFPbandpasses, i, ch, rip_num, rippleStats.fs, rippleband)

                    end
                    orient(h3,'landscape')
                    savepdf(h3, sprintf('%s/PDF/%s_%s_exampleRipple_%i.pdf', singleTrialFolder, subject, num2str(chan_labels(ch)), f))
                    saveas(h3,  sprintf('%s/PNG/%s_%s_exampleRipple_%i.png', singleTrialFolder, subject, num2str(chan_labels(ch)), f))
                    saveas(h3,  sprintf('%s/Fig/%s_%s_exampleRipple_%i.fig', singleTrialFolder, subject, num2str(chan_labels(ch)), f))
                    close
                end

                close all
            end
        end
    end
end