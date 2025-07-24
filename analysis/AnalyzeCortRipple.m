clc
clear 
close all

% CC08_LA_N2_ULe3-2_CortRipple_MeanNC_TimeFreq_Comp
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
addpath(genpath('/home/iverzh/ripple/CortRipple'))
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/eeglab14_1_2b/'))
rmpath('/space/seh8/1/halgdev/projects/cdickey/eeglab14_1_2b/functions/octavefunc')
%% Inputs

% sEEG subjects from Cleveland Clinic
subj_list_full = {'CC04','CC08','CC15','CC17','CC18','CC20','CC23','CC24','CC25','CC26',...
                  'CC30','CC31','CC39','CC49','CC55','CC60','CC69','CC91','CC92'};

win = 2000; % pre and post event window       

spikeCheck = 100; %+/- for spike check 
rippleWindow = 10; %+/- window used to compute ripple ocsillation freq. 
fs = 1000; %Hz sample rate

masterFolder = '/home/iverzh/ripple/CortRipple/Figures_GoodRipples';
if ~isfolder(masterFolder); mkdir(masterFolder); end

matExportFolder = '/space/seh9/5/halgdev/projects/iverzh/ripples/matFiles';
if ~isfolder(matExportFolder); mkdir(matExportFolder); end

RejectParams.RBHFzscore = 1; %Reject zscore(Ripple) - zscore(HF) < RBHFzscore
RejectParams.RBLFzscore = 1; %Reject zscore(Ripple) - zscore(LF) < RBLFzscore
RejectParams.LFPdiffzscore = 7; %Reject zscore(LFPdiff) > LFPdiffzscore 
RejectParams.RBHFratio = 1; %Reject RBamp/HFamp < 1
RejectParams.RBAmp = 2; %Reject RBAmp < 2

%%
% NREM sleep files selected for analaysis
% sleepfiles_set = {[2,3,7],[2,5,9,10],[3,5,6,7],[1],...
%         [2,3,5],[4,10,12],[5,6,7,9,11],[2,3,6,8],[1,2,4,5,7],...
%         [3,7,8],[3,7,8,11],[1,2,5],[2,4,5,6,7,8,10,11],[3,8,9,12],...
%         [2,6,8],[2,6,7],[1,2,3,5,7,10],[1,3,4,5],[1,3,4]};

sleepfiles_set = {[2,3,7],[2,5,9],[3,5,6,7],[1],...
        [2,3,5],[4,10,12],[5,6,7,9,11],[2,3,6,8],[1,2,4,5,7],...
        [3,7,8],[3,7,8,11],[1,2,5],[2,4,5,6,7,8,10,11],[3,8,9,12],...
        [2,6,8],[2,6,7],[1,2,3,5,7,10],[1,3,4,5],[1,3,4]};
    
sleepStages = 2:3; %2 = N2, 4 = N3 


for subj = [2] %1:numel(subj_list_full) %loop through subjects
    subject = subj_list_full{subj}

    for s = 1%:length(sleepfiles_set{subj})
        sleep_ID = sleepfiles_set{subj}(s)
        
        
        load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/segment_indices/sleep_%01d_N23_seg_NC_polchecked_newRipp.mat', subject,  sleep_ID))
        load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/NC_ripple/%s/segment_indices/sleep_%01d_N23_seg_NC.mat', subject, sleep_ID));
        load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/N23_masks/%s_%02d_N23_mask.mat', subject, sleep_ID));

        ripple = struct([]);

        % notch filter and ripple bandpass data
        ripple(1).rippleband = NaN(size(sleepdata));
        for ch = 1:size(sleepdata,1) % loop through channels
            % notch filter data to remove line noise (60 Hz + harmonics)
            for nf = 60:60:240
                Wo = nf/(fs/2);
                BW = Wo/35;
                [b,a] = iirnotch(Wo, BW);
                sleepdata(ch,:) = filtfilt(b,a,sleepdata(ch,:));
            end

            % bandpass at ripple band (60-120 Hz)
            [b,a] = butter(3,[60 120]/(fs/2));
            ripple.rippleband(ch,:) = filtfilt(b,a,sleepdata(ch,:)); % zero phase

            % bandpass at low frequency band (25-50 Hz) 
            [b,a] = butter(3,[25 50]/(fs/2));
            ripple.lowfreqband(ch,:) = filtfilt(b,a,sleepdata(ch,:)); % zero phase

            % bandpass at high frequency band (200- Hz) 
            [b,a] = butter(3,200/(fs/2), 'high');
            ripple.highfreqband(ch,:) = filtfilt(b,a,sleepdata(ch,:)); % zero phase
            a = sprintf('notch and bandpass filter channel %i\n', ch);
            fprintf(a)
            
            ripple.RBzscore(ch,:) = zscore(abs(hilbert(ripple.rippleband(ch,:)))); %zscore of env of hilbert amplitudes
            ripple.HFzscore(ch,:) = zscore(abs(hilbert(ripple.highfreqband(ch,:))));
            ripple.LFzscore(ch,:) = zscore(abs(hilbert(ripple.lowfreqband(ch,:))));
            ripple.diffLFPzscore(ch,:)  = zscore(abs(diff(sleepdata(ch,:)))); %zscore of LFP numerical derivative 
        end

        ripple.LFP = sleepdata; 
        clear sleepdata;


        
        for stage = 1:length(sleepStages)
            stageMask = find(N23_mask_all{stage});
            NREM = find(NREM_mask);
            [~, stageMask, ~] = intersect(NREM, stageMask);
           
            
                for ch = 1:size(ripple.LFP,1)
                    
                   inStage = ismember(good_ripple_locs_set{ch}, stageMask);
                   goodRipples = good_ripple_locs_set{ch}(inStage);
                   
                   ripple.ind{ch} = goodRipples(goodRipples-win > 0 & goodRipples+win < size(ripple.LFP,2));

                   ripple.events.raw{ch}               = NaN(numel(ripple.ind{ch}), win*2 + 1);
                   ripple.events.rippleband{ch}        = NaN(numel(ripple.ind{ch}), win*2 + 1);
                   ripple.events.highFreqBand{ch}      = NaN(numel(ripple.ind{ch}), win*2 + 1);
                   ripple.events.lowFreqBand{ch}       = NaN(numel(ripple.ind{ch}), win*2 + 1);

                   ripple.events.RBzscore{ch}          = NaN(numel(ripple.ind{ch}), win*2 + 1);
                   ripple.events.HFzscore{ch}          = NaN(numel(ripple.ind{ch}), win*2 + 1);
                   ripple.events.LFzscore{ch}          = NaN(numel(ripple.ind{ch}), win*2 + 1);
                   ripple.events.diffLFPzscore{ch}     = NaN(numel(ripple.ind{ch}), win*2 + 1);  

                   ripple.events.highFreqAmp{ch}       = NaN(numel(ripple.ind{ch}),1);
                   ripple.events.lowFreqAmp{ch}        = NaN(numel(ripple.ind{ch}),1);
                   ripple.events.rippleAmp{ch}         = NaN(numel(ripple.ind{ch}),1);
                   ripple.events.highRippleRatio{ch}   = NaN(numel(ripple.ind{ch}),1);
                   ripple.events.lowRippleRatio{ch}    = NaN(numel(ripple.ind{ch}),1);
                   ripple.events.cov{ch}               = NaN(numel(ripple.ind{ch}),1);


                   fprintf(['extracting ripples and ripple stats from channel : ',chan_labels{ch}, '\n'])      


                   for rip_num = 1:numel(ripple.ind{ch}) % loop through ripple events

                        ripple.events.raw{ch}(rip_num,:)          = ripple.LFP(ch,ripple.ind{ch}(rip_num)-win:ripple.ind{ch}(rip_num)+win);
                        ripple.events.rippleband{ch}(rip_num,:)   = ripple.rippleband(ch,ripple.ind{ch}(rip_num)-win:ripple.ind{ch}(rip_num)+win);
                        ripple.events.highFreqBand{ch}(rip_num,:) = ripple.highfreqband(ch,ripple.ind{ch}(rip_num)-win:ripple.ind{ch}(rip_num)+win);
                        ripple.events.lowFreqBand{ch}(rip_num,:)  = ripple.lowfreqband(ch,ripple.ind{ch}(rip_num)-win:ripple.ind{ch}(rip_num)+win);

                        ripple.events.HFzscore{ch}(rip_num,:)     = ripple.HFzscore(ch,ripple.ind{ch}(rip_num)-win:ripple.ind{ch}(rip_num)+win);   
                        ripple.events.LFzscore{ch}(rip_num,:)     = ripple.LFzscore(ch,ripple.ind{ch}(rip_num)-win:ripple.ind{ch}(rip_num)+win);   
                        ripple.events.RBzscore{ch}(rip_num,:)     = ripple.RBzscore(ch,ripple.ind{ch}(rip_num)-win:ripple.ind{ch}(rip_num)+win);   
                        ripple.events.diffLFPzscore{ch}(rip_num,:)= ripple.diffLFPzscore(ch,ripple.ind{ch}(rip_num)-win:ripple.ind{ch}(rip_num)+win);   

                        % compute ripple stats
                        center = size(ripple.events.raw{ch},2)/2;
                        center = round(center);
                        rippleRange = center-rippleWindow:center+rippleWindow;

                        %High Frequency Hilbert Amp
                        signalhighF =ripple.events.highFreqBand{ch}(rip_num,rippleRange);
                        ripple.events.highFreqAmp{ch}(rip_num)       = max(abs(hilbert(signalhighF)));

                        % Low Frequency Ripple Amp
                        signallowF =ripple.events.lowFreqBand{ch}(rip_num,rippleRange);
                        ripple.events.lowFreqAmp{ch}(rip_num)        = max(abs(hilbert(signallowF)));

                        % Ripple Band Amp
                        signalRipple =ripple.events.rippleband{ch}(rip_num,rippleRange);
                        ripple.events.rippleAmp{ch}(rip_num)         = max(abs(hilbert(signalRipple)));

                        % Ratio RB amp : HF amp
                        ratio_HighRipple = max(abs(hilbert(signalRipple))) / max(abs(hilbert(signalhighF)));
                        ripple.events.highRippleRatio{ch}(rip_num)   = ratio_HighRipple;

                        % Ratio RB amp: LF amp
                        ratio_LowRipple = max(abs(hilbert(signalRipple))) / max(abs(hilbert(signallowF)));            
                        ripple.events.lowRippleRatio{ch}(rip_num)    = ratio_LowRipple;
                   end

                end


                cleanedEvents = rejectRipples(ripple, RejectParams, rippleWindow);
                plotStats =  computePlotStats(cleanedEvents, fs);

    %% Plotting

            for ch = 1:size(ripple.LFP,1)
               
                exportFolder = fullfile(masterFolder, subject, ['sleep_',num2str(sleep_ID)], ['N',num2str(sleepStages(stage))],chan_labels{ch});
                if ~isfolder(exportFolder); mkdir(exportFolder); end

                if numel(cleanedEvents.rippleAmp{ch}) > 40
                    ripple.meanLFP(ch,:) = mean(cleanedEvents.raw{ch});

                    for rip_num = 1:size(cleanedEvents.raw{ch},1)
                       rippleCOV = cov(cleanedEvents.raw{ch}(rip_num,:), ripple.meanLFP(ch,:));
                       cleanedEvents.cov{ch}(rip_num) = rippleCOV(1,2);     
                    end

%                     hsummary = chanOverviewPlot(cleanedEvents, plotStats, ch, ripple, win, fs, subject, chan_labels, sleepStages(stage));
%                     saveas(hsummary, sprintf('%s/%s_sleep%01d_%s_summary.png', exportFolder, subject, sleep_ID, chan_labels{ch}))
%                     
                    hTF = chanTFplot(cleanedEvents, plotStats, ch, ripple, win, fs, subject, chan_labels, sleepStages(stage),'TF');
                    saveas(hTF, sprintf('%s/../../%s_N%01d_%s_CortRipple_MeanNC_TimeFreq_Comp.fig', exportFolder, subject, sleepStages(stage), chan_labels{ch}))
                    
                    hTF = chanTFplot(cleanedEvents, plotStats, ch, ripple, win, fs, subject, chan_labels, sleepStages(stage),'LFP');
                    saveas(hTF, sprintf('%s/../../%s_N%01d_%s_CortRipple_MeanNC_LFP_Comp.fig', exportFolder, subject, sleepStages(stage), chan_labels{ch}))

                    [A,I] = sort(cleanedEvents.rippleAmp{ch}, 'descend');


                    randInd = randperm(size(cleanedEvents.raw{ch},1));
% 
%                     for f = 1:10 
%                         h3 = figure('Position', [100 100 3000 3000]);
%                         for i = 1:4
% 
%                             rip_num = I(4*(f-1) + i);
%                             rippleFullPlot(cleanedEvents, i, ch, rip_num)
% 
%                         end
%                         saveas(h3, sprintf('%s/%s_sleep%01d_%s_exampleRipples_%i.png', exportFolder, subject, sleep_ID, chan_labels{ch}, f))
%                         close
%                     end

%                     h4 = figure(4);
%                     histogram(cleanedEvents.lowRippleRatio{ch}, 100, 'EdgeColor', 'None', 'FaceColor', 'b', 'FaceAlpha', 0.5); hold on;
%                     histogram(cleanedEvents.highRippleRatio{ch}, 100, 'EdgeColor', 'None','FaceColor', 'r', 'FaceAlpha', 0.5); hold on;
%             %         hist2.FaceColor = 'r';
%                     l = legend('lowFreq', 'highFreq');
%                     saveas(h4, sprintf('%s/%s_sleep%01d_%s_RatioHist.png', exportFolder, subject, sleep_ID, chan_labels{ch}))



                    close all
                end
            end
            rippleStats = plotStats;
            for ch = 1:size(ripple.LFP,1)
                rippleStats.locs{ch} = cleanedEvents.goodRipples{ch}(plotStats.centeredInd{ch});
            end
%             save(sprintf('%s/%s_sleep%01d_N%istats.mat', matExportFolder, subject, sleep_ID, sleepStages(stage)), 'rippleStats', '-v7.3')
        
        end
    end
end