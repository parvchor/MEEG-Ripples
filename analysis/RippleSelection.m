clc
clear 
close all

addpath(genpath('/space/seh10/5/halgdev/projects/swilson/training/code/RippleDetection_clean'))

addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))

addpath(genpath('/home/jgarret/mat')) %mask2bounds
addpath(genpath('/home/jgarret/ArtifactDetection'))

% addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/eeglab14_1_2b/'))
% rmpath(genpath('/space/seh8/1/halgdev/projects/cdickey/eeglab14_1_2b/functions/octavefunc'))

masterFolder = '/space/seh10/5/halgdev/projects/swilson/training/data/ripples/Figures';
if ~isfolder(masterFolder); mkdir(masterFolder); end

matExportFolder = '/space/seh10/5/halgdev/projects/swilson/training/data/rippleout';
if ~isfolder(matExportFolder); mkdir(matExportFolder); end

preProcExportFolder = '/space/seh10/5/halgdev/projects/swilson/training/data/preproc';
if ~isfolder(preProcExportFolder); mkdir(preProcExportFolder); end

%% Inputs
lfpDirec = '/space/seh10/5/halgdev/projects/swilson/training/data/lfp_data';

subj_list_filtered = {'SUBJ1'};
runList = 1:length(subj_list_filtered); 

recordingState = 'wake'; % for cc wake, use sleep established polarity check
locations = {'ALL'}; %{'NC', 'TH'};%,'HC','TH'}; %{'NC','HC','TH','AMY'};
modifier = '';

arrayName = 'lateral';
if strcmp(recordingState,'sleep') || strcmp(recordingState,'REM')
    isSleep = 1;
elseif strcmp(recordingState,'wake')
    isSleep = 0;
else
    error("incorrect RecordingState argument. Can only be 'sleep', 'waking', 'REM'")
end

badChanDir = '';
scale = 1; 
IISflag = 1; %1 - check for spikes. 0 - do not check for spikedf,s

winSec = 2; % 2 human 1 rodent % pre and post event window in seconds
rippleWindowSec = 0.100; % 0.100; human 0.050 rodent %+/- window used to compute ripple ocsillation freq. in seconds
IIScheckSec = 0.500; %+/- in sec

fs = 1000; % recording data sampling rate

rippleband = [70 100]; % human [70 100] rodent [120 200]
species = 'human'; % 'human' or 'rodent'

%[theshold value, binary flag (0-inactive, 1-active)]
RejectParams.RBHFzscore = [0, 0]; %Reject zscore(Ripple) - zscore(HF) < RBHFzscore
RejectParams.sharpzscore = [7, 1]; %Reject zscore(100+hZ) > sharpzscore
RejectParams.LFPdiff = [50, 0]; %Reject zscore(LFPdiff) > LFPdiffzscore 
RejectParams.LFPdiffzscore = [4, 0]; %Reject stepwise jumps in UAB data
RejectParams.RBzscore = [3, 1]; %Reject RBAmp < 2 original [3, 1]
RejectParams.LFPthresh = [1000, 1];
RejectParams.minDuration = [0.015, 1]; % in seconds - consider 0.015 for rodent, 0.025 for human

RejectParams.bndParams.smthKrnl=100;
RejectParams.bndParams.srchWinSz=1/2; %of kernel size
RejectParams.bndParams.scoreThresh=0.75;

% NREM sleep files selected for analaysis
% sleepfiles_set = {[2,3,7],[2,5,9,10],[3,5,6,7],[1],...
%         [2,3,5],[4,10,12],[5,6,7,9,11],[2,3,6,8],[1,2,4,5,7],...
%         [3,7,8],[3,7,8,11],[1,2,5],[2,4,5,6,7,8,10,11],[3,8,9,12],...
%         [2,6,8],[2,6,7],[1,2,3,5,7,10],[1,3,4,5],[1,3,4]};
% 


sleepfiles_set_other = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};


    
    %%
fprintf('Running Ripple Selection ...\n')

for subj = runList %[1,3,6,7,10,11,13,16,17] %%  %[14, 15] %[11,13,14,15] %21 %15:17 %[3,4,5,8,10,11,13,14,15]  %loop through subjects
    subject = subj_list_filtered{subj};
    
    fprintf('Loading Subject %s ...\n', subject)

        
    if ~contains(subject,'CC') || ~strcmp(recordingState,'sleep')
        sleepfiles_set = sleepfiles_set_other;
    else
        sleepfiles_set = sleepfiles_set_CC;
    end

    if contains(subject,'CC')
        parcelFiles = dir(fullfile(parcelDir,['*',subject,'*.mat']));
        load(fullfile(parcelFiles(1).folder, parcelFiles(1).name)) 
    end
    
    badChanFile = fullfile(badChanDir,[subject,'_bad_chan.mat']);
    if isfile(badChanFile); load(badChanFile);
    else; bad_chan = cell(0);
    end
    

    for l = 1:length(locations)

        location = locations{l};
        
        fprintf('Loading %s Data...\n', location)
        
        tag = [recordingState,'_',location,'_',modifier];
        rippleStatsAllSleepSet = struct([]);
        shift = 0;
        saveRippleStats = true;
        for s = 1:length(sleepfiles_set{subj})
            sleep_ID = sleepfiles_set{subj}(s);
            
            fprintf('Loading segment %i ...\n', sleep_ID)               
            
            try
                f = sprintf('%s/%s.mat', lfpDirec, subject);
                load(f)
    %                     load(fripp)
            catch
               WarningString = ['could not load ',f];                   
               warning(WarningString)
               saveRippleStats = false;
       
                break
            end
            
            if exist('sfreq','var')
                fs = sfreq;
            end
            
            if exist('all_ripple_locs','var')
                good_ripple_locs_set = all_ripple_locs;
            end

            if exist('data','var')
                BroadbandData = data;
                clear data
            elseif exist('sleepdata','var')
                BroadbandData = sleepdata;
                clear sleepdata
            elseif exist('datNREM','var')
                BroadbandData = datNREM';
            elseif exist('datREM','var')
                BroadbandData = datREM';
            elseif exist('lfp', 'var')
                BroadbandData = lfp;
            end

            if exist('channel', 'var')
                channel = channel;
            end

             if exist('label', 'var')
                channel = label;
            end

            % nan out bad time points
            try 
                load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/%s/nan_edge_mask.mat', subject, recordingState),'nan_edge_mask') 
            catch 
                nan_edge_mask = false(1,size(BroadbandData,2));
            end

            if ~isa(BroadbandData, 'double'); BroadbandData = double(BroadbandData); end

            if size(BroadbandData,1) > size(BroadbandData,2); BroadbandData = BroadbandData'; end
            
            BroadbandData = scale * BroadbandData;

            if strcmp(recordingState,'sleep')
                polCheckFiles = dir(fullfile(preProcExportFolder,[subject,'_',num2str(sleep_ID),'_',tag,'_polCheck.mat']));

                if isempty(polCheckFiles)
                    polCheck = checkPolarity(BroadbandData,fs); %check polarity by comparing slow oscillations w HG
                    try
                        save(fullfile(preProcExportFolder,[subject,'_',num2str(sleep_ID),'_',tag,'_polCheck.mat']), 'polCheck', '-v7.3')
                    catch
                        warning('could not save polarity check')
                    end
                else
                    load(fullfile(preProcExportFolder,[subject,'_',num2str(sleep_ID),'_',tag,'_polCheck.mat']))
                end
            else
                if size(BroadbandData,1) == 96
                    polCheck = -1*ones(size(BroadbandData));
                else
                    polCheck = ones(size(BroadbandData));
                end

            end

            BroadbandData = polCheck .* BroadbandData;

            if IISflag
                spikeCheckFiles = dir(fullfile(preProcExportFolder,[subject,'_',num2str(sleep_ID),'_',tag,'_spikeCheck.mat']));
                if isempty(spikeCheckFiles)
                    [spikeMask,spikeBounds] = detectSpikes_IAV(BroadbandData', fs);
                    try
                        save(fullfile(preProcExportFolder,[subject,'_',num2str(sleep_ID),'_',tag,'_spikeCheck.mat']),'spikeBounds','spikeMask', '-v7.3')
                    catch
                        warning('could not save spike check')
                    end
                else
                    load(fullfile(preProcExportFolder,[subject,'_',num2str(sleep_ID),'_',tag,'_spikeCheck.mat']))
                end
            else 
                spikeMask = zeros(size(BroadbandData));
            end

            if isa(channel,'double'); channel = cellstr(string(channel)); end

            for ch = 1:size(BroadbandData,1)
%                 initialRipples = good_ripple_locs_set{ch}; % comment this out if you want to bypass Xi's method
                   
                [initialRipples, ~, ~]  = DetectHighFreqEvents(BroadbandData(ch,:),rippleband,RejectParams,fs,nan_edge_mask,3,1);
                fsNew = fs;
                win = round(winSec * fsNew); %win in samples.
                rippleWindow = round(rippleWindowSec * fsNew); %ripple window in samples
                IIScheck = round(IIScheckSec * fsNew); %spike check window in samples 
                
                checkDiff = [diff(initialRipples) 5];
                initialRipples(checkDiff<5) = [];
                [rippleSleepSet,rejectVec] = AnalyzeRipple(species, BroadbandData, rippleband, win, rippleWindow, IIScheck, spikeMask, nan_edge_mask, IISflag, ...
                                                          ch, initialRipples, RejectParams, fsNew, channel, shift,1);

                 
                rippleStatsSleepSet =  computeRippleStats(rippleSleepSet, RejectParams, fs);

                rippleStatsSleepSet.locs{1} = rippleSleepSet.goodRipples;
                rippleStatsSleepSet.rejectVec{1} = rejectVec;
                density = length(rippleStatsSleepSet.locs{1}) / length(BroadbandData) * fs * 60;
                rippleStatsSleepSet.density{ch} = density;
                
                rippleStatsSleepSet = patchRippleStats(rippleband, rippleStatsSleepSet, BroadbandData(ch,:), [], fs, 0, 0);

               
                if isempty(rippleStatsAllSleepSet)
                    fn = fieldnames(rippleStatsSleepSet);
                    for ch2 = 1:length(channel)
                        for f = 1:numel(fn)
                            rippleStatsAllSleepSet(1).(fn{f}){ch2} = [];
                        end
    
                    end
                end
    
    
                fn = fieldnames(rippleStatsSleepSet);
                if numel(rippleStatsSleepSet.locs{1}) > 5
                    for f = 1:numel(fn)
                        if isa(rippleStatsSleepSet.(fn{f}), 'cell')
                            fDat = rippleStatsSleepSet.(fn{f}){1};
                            DIM = size(fDat,1);
                            if DIM > 1 || size(fDat,2) == 4001
                                rippleStatsAllSleepSet.(fn{f}){ch} = [rippleStatsAllSleepSet.(fn{f}){ch}; fDat];
                            elseif DIM == 1
                                rippleStatsAllSleepSet.(fn{f}){ch} = [rippleStatsAllSleepSet.(fn{f}){ch}, fDat];
                            end
                        elseif isa(rippleStatsSleepSet.(fn{f}), 'struct')
                            fDat = rippleStatsSleepSet.(fn{f});
                           rippleStatsAllSleepSet.(fn{f}) = fDat;
                        end
    
                    end
                end
    
                rippleStatsAllSleepSet.rejectVec{ch} = rejectVec;
            end

            rippleStatsAllSleepSet.recordingLength(s) = size(BroadbandData,2);

            locs_sleep = nan(1, 5e5); %matrix to store ripple locs by sleep set
            for ch = 1:size(BroadbandData,1)
                       %save ripple indicies by sleep set
                    locs_sleep(1:length(rippleStatsAllSleepSet.locs{ch})) = rippleStatsAllSleepSet.locs{ch};
                    rippleStatsAllSleepSet.locs_sleep{ch}(s,1:size(locs_sleep,2)) = locs_sleep;
            end

            shift = shift + size(BroadbandData,2); %shift for ripple indexing


        end

        %% Exporting Ripple Data
        if saveRippleStats
            rippleStats = rippleStatsAllSleepSet;
            for ch = 1:size(BroadbandData,1)              

                density = length(rippleStatsAllSleepSet.locs{ch}) / sum(rippleStats.recordingLength) * fs * 60;
                rippleStats.density{ch} = density;
            end

            rippleStats.fs = fs;
            rippleStats.chanLabels = channel;
            rippleStats.RejectParams = RejectParams;
            rippleStats.RB = rippleband;
            rippleStats.IISflag = IISflag;
        %         save(sprintf('%s/%s_N%istats.mat', matExportFolder, subject, sleepStages(stage)), 'rippleStats', '-v7.3')  
        save(sprintf('%s/%s_ripple_stats_%s.mat', matExportFolder, subject, tag), 'rippleStats', '-v7.3')
        end
    end
    
    
end
    
    
    
   