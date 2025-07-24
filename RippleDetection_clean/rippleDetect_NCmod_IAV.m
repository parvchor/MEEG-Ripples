% standalone ripple detector
% dependency: Fieldtrip 2014 and above; "newtimef" from EEGLAB
% new methods: hp_freq_ST,lp_freq_ST,order_ST,hasRMS,RMSfile,plotCtx,plotHipp
% note that this only accepts one hipp part at a time (but could have multiple channels for same hipp part in order to reject spikes)
function rippleDetect_NCmod_IAV(data,sfreq,savedir,sleepID,ctx_ind,chan_labels,stage,methods, N23_mask, state)

delete(gcp('nocreate'));
disp(['Do wavelet: ' num2str(methods.do_wavelet)])
tic
%% Set up default inputs
if isempty(methods)
%     methods.hp_freq_ST = 70;
%     methods.lp_freq_ST = 100;
    methods.order_ST = 6;
    methods.hasRMS = 0;
    methods.plotCtx = 0;
    methods.plotHipp = 0;
    
    methods.time_bin = 2;           % size of time bins, in seconds
    methods.wavthresh_Haar = 6;
    methods.wavthresh_RL1 = 6;
    methods.wavthresh_RL2 = 8;
    methods.wavthresh_RL3 = 8;
    methods.wavthresh_RL4 = 5;
    methods.wavthresh_RL5 = 4;
    methods.binsize = methods.time_bin;
    methods.bin_size = methods.time_bin;
    methods.use_relative_thresh = 1;
    methods.clean_prior = 1;
    methods.sfreq = sfreq;
    
    methods.Q_thresh = 0.9;
    methods.ZC_bin = 0.5;  % in seconds; for detrending and peak selection 
    methods.ZC_lowpass = 120;
    methods.ZC_chan = 1;  % the channel to define zero crossings and biphasic components from (and actually define ripples on)
    methods.IRI = 0.025; % in seconds; minimal distance between two ripple events
    methods.adjacent_merge = 1;
    methods.adjacency = 25;
    
    methods.separation_flag = 0;
    methods.separation_thresh_dot = 0.05;  % in 0.01*percentile
    methods.separation_thresh_pmax = 0.1;  % in 0.01*percentile
    methods.hipp_polarity = 0;
    
    % optional parameters
    methods.useAlterThresh = 0;  
    methods.ripple_chan = 1;  
    methods.biphasic_template = [];  % necessary if separation_flag == 1
    methods.HM_locs = [];            % necessary if separation_flag == 1 
    methods.HM_seg = [];             % necessary if separation_flag == 1
    
    methods.rmDotBase = 0;
    methods.looseEval = 0;
    methods.gap_thresh = 0.04;    % for SWR only
    methods.do_wavelet = 1;
    methods.do_pmax = 1;
    methods.HP_reject = 0;   % only useful when do_wavelet == 1
    methods.HP_thresh = 3;   % only useful when HP_reject == 1, noting # of std above mean
    
    methods.parpool_size = 12;  % parallel forloop for wavelet rejection
end

if ~isempty(ctx_ind)
%     data_ctx = data(1:ctx_ind,:);
%     chan_labels = chan_labels(1:ctx_ind);
    data = data(1:ctx_ind,:);
%     methods.ripple_chan=1:ctx_ind; % Charlie fix
    chan_labels = chan_labels(1:ctx_ind);
end
% if isfield(methods,'ripple_chan')
%     data_full = data(methods.ripple_chan,:);
% end
data_full = data; 

num_chans_rej = size(data,1);
num_chans = num_chans_rej;  % by default, if ANY channel has an artifact within a given time bin, anything detected in other channels during that time bin is discarded
%methods.parpool_size = num_chans;


% myCluster = parcluster('local');
% myCluster.NumWorkers = methods.parpool_size;
% saveProfile(myCluster);
%
% uncomment this if wish to skip files already made
% clear all_ripple_locs_concat
% try
%     load(sprintf('%s/sleep_%d_%s_polchecked_newRipp.mat',savedir,sleepID,stage),'all_ripple_locs_concat')
% catch
%     disp('OK')
% end
% if exist('all_ripple_locs_concat','var')
%     disp('AlreadyDone')
%     return
% end

%% Set up detection parameters (can skip if hasRMS == 1)
if methods.hasRMS ~= 1
    % Band pass filter for ZC check
    filt_set = zeros(num_chans,length(data));
    for iChan = 1:num_chans
        %filt_set(iChan,:) = ft_preproc_lowpassfilter(data(iChan,:),sfreq,methods.ZC_lowpass,methods.order_ST/2,'but');
        [b,a] = butter(3,methods.ZC_lowpass/(sfreq/2), 'low');
        filt_set(iChan,:) = filtfilt(b,a,data(iChan,:));
    end
    
    data_filtered_set = zeros(num_chans,length(data));
    for iChan = 1:num_chans
        %disp('high pass filtering')
        %data_filtered_ST = ft_preproc_highpassfilter(data(iChan,:),sfreq,methods.hp_freq_ST,methods.order_ST/2,'but');
        [b,a] = butter(3,methods.hp_freq_ST/(sfreq/2), 'high');
        data_filtered_ST = filtfilt(b,a,data(iChan,:));
        
        data_filtered_ST = abs(data_filtered_ST);
        %disp('low pass filtering');
        [b,a] = butter(3,methods.lp_freq_ST/(sfreq/2), 'low');
        data_filtered_ST = filtfilt(b,a,data_filtered_ST);
        %data_filtered_ST=ft_preproc_lowpassfilter(data_filtered_ST,sfreq,methods.lp_freq_ST,methods.order_ST/2,'but');
        data_filtered_set(iChan,:) = data_filtered_ST;
%     save(sprintf('%s/sleep_%d_%s.mat',savedir,sleepID,stage),'data_filtered_ST','-v7.3')
    end
    
    % get RMS with 20ms moving window
    %delete(gcp('nocreate'));
    rms_set = cell(1,num_chans);
    if num_chans > 1
        %myPool = parpool(8);
        for i = 1:num_chans % CD: deleted parfor loop 
            rms_set{i} = rms_mov(data_filtered_set(i,:),round(sfreq*0.02),round(sfreq*0.02)-1,0);
        end
        %delete(myPool);
    else
        rms_set{1} = rms_mov(data_filtered_set(1,:),round(sfreq*0.02),round(sfreq*0.02)-1,0);
    end
    
    % find threshold and select peaks
    rms_thresh = zeros(1,num_chans);
    if methods.useAlterThresh == 1 % assumes that methods.threshfile exists
        load(methods.threshfile,'rms_thresh')
    end
    
    rms_peakloc = cell(1,num_chans);
    peakeval = rms_peakloc;
    
% %     figure(100) % CD added
    for i = 1:num_chans
        if methods.useAlterThresh ~= 1
            rms_thresh(i) = quantile(rms_set{i},methods.Q_thresh);
        end
        rms_peakloc{i} = find(rms_set{i} >= rms_thresh(i));
        test = find(diff(rms_peakloc{i}) > 1);
        peakeval_curr = zeros(1,length(rms_peakloc{i}));
               
        for j = 1:length(test)
            if j == 1
%                 if rms_peakloc{i}(test(j)) - rms_peakloc{i}(1) >= round(sfreq*0.038) % Staresina limit for duration exceeding threshold
                if round((test(j)-1)/2) <= 0
                    peakeval_curr(1) = 1;
                else
                   peakeval_curr(round((test(j)-1)/2)) = 1;
                end
%                 end
            else
%                 if rms_peakloc{i}(test(j)) - rms_peakloc{i}(test(j-1)+1) >= round(sfreq*0.038) % Staresina limit for duration exceeding threshold
                   peakeval_curr(test(j-1)+1+round((test(j)-(test(j-1)+1))/2)) = 1;
%                 end
%                   title(num2str(j))
%                   if j > 5
%                       plot(-100:100, rms_set{i}(rms_peakloc{i}(test(j-1)+1+round((test(j)-(test(j-1)+1))/2))-100:rms_peakloc{i}(test(j-1)+1+round((test(j)-(test(j-1)+1))/2))+100));
%                       xlim([-100 100])
%                       vline(0)
%                       waitforbuttonpress
%                   end
            end
            if j == length(test)
%                 if rms_peakloc{i}(end) - rms_peakloc{i}(test(j)+1) >= round(sfreq*0.038)
                    peakeval_curr(test(j-1)+1+round((length(rms_peakloc{i})-(test(j)+1))/2)) = 1;
%                 end
            end
            
        end
        peakeval{i} = peakeval_curr;
        %disp('Herrrring')
    end
    
    save(sprintf('%s/%s_%s.mat',savedir,state, num2str(sleepID)),'rms_peakloc','rms_set','peakeval','rms_thresh', ...
                                                                 'data', 'chan_labels', 'N23_mask', 'sfreq','-v7.3','-nocompression')
    disp(['Sleep ',num2str(sleepID),' putative ripples detected (ST)!'])
else
    load(methods.RMSfile,'rms_peakloc','rms_set','peakeval','data')
    filt_set = zeros(num_chans,length(data));
    for iChan = 1:num_chans
        %filt_set(iChan,:) = ft_preproc_lowpassfilter(data(iChan,:),sfreq,methods.ZC_lowpass,methods.order_ST/2,'but');
        [b,a] = butter(3,methods.ZC_lowpass/(sfreq/2), 'low');
        filt_set(iChan,:) = filtfilt(b,a,data(iChan,:));
    end
end

toc 
%% Artifact rejection
disp('Beginning artifact rejection.')
all_ripple_locs = cell(1,num_chans);
for i = 1:num_chans
    all_ripple_locs{i} = rms_peakloc{i}(peakeval{i}>0)+round(0.5*sfreq*0.02); % position correction for moving rms
    
%     title(num2str(j))
%     pk_ind = find(peakeval{1});
%     for tl = all_ripple_locs{i}
%         if tl > 100
%           plot(-100:100, data(i,tl-100:tl+100));
%           xlim([-100 100])
%           vline(0)
%           waitforbuttonpress
%         end
%     end
end
% all_ripple_locs_concat = unique(cell2mat(all_ripple_locs));
    
putative_ripples_set = all_ripple_locs;
all_ripple_ranges_set = cell(1,num_chans);
put_ripple_set = cell(1,num_chans); ZC_outputs_set = cell(1,num_chans);
for iChan = 1:num_chans
    all_ripple_ranges_set{iChan} = zeros(length(putative_ripples_set{iChan}),2);
    all_ripple_ranges_set{iChan}(:,1) = putative_ripples_set{iChan} - round(0.5*sfreq)+1;
    all_ripple_ranges_set{iChan}(:,2) = putative_ripples_set{iChan} + round(0.5*sfreq);
end
    
for iChan = 1:num_chans
    putative_ripples = putative_ripples_set{iChan};
    all_ripple_ranges = all_ripple_ranges_set{iChan};
    ripple_qual = zeros(1,length(putative_ripples));  
    for iRipple = 1:length(all_ripple_ranges)
        [logi_idx,num_idx] = ismember(putative_ripples,[all_ripple_ranges(iRipple,1):all_ripple_ranges(iRipple,2)]);
        if length(find(num_idx > 0)) > 1
            good_locs = find(logi_idx);
            for j = 2:length(find(num_idx > 0))
                if j == 2
                    loc_comp = good_locs(j-1);
                end
                gap = num_idx(good_locs(j)) - num_idx(loc_comp);
%                 gap = num_idx(j) - num_idx(j-1);
                if ripple_qual(good_locs(j-1)) ~= -1                    
                    ripple_qual(good_locs(j-1)) = 1;
                end
                if gap >= (methods.IRI * sfreq) && ripple_qual(good_locs(j)) ~= -1 % if distance between consecutive high gamma peaks are at least (by default) 50 ms
                    ripple_qual(good_locs(j)) = 1;  % count both as separate ripples
                    loc_comp = good_locs(j);
                else
                    ripple_qual(good_locs(j)) = -1;
                end
            end
        else
            ripple_qual(logi_idx) = 1;
        end
    end

    all_ripple_locs_concat = putative_ripples(ripple_qual == 1);
    if methods.looseEval == 1  % check the eval scripts for specific parameters that go into zero-crossing-based frequency estimation
        [~,ZC_outputs.PKS_set,ZC_outputs.LOCS_set,ZC_outputs.evaluator_set,~] = CC_ripple_evaluator_loose(all_ripple_locs_concat,filt_set(iChan,:));
    else
        [~,ZC_outputs.PKS_set,ZC_outputs.LOCS_set,ZC_outputs.evaluator_set,~] = CC_ripple_evaluator(all_ripple_locs_concat,filt_set(iChan,:));
    end
    all_ripple_locs_concat = all_ripple_locs_concat(ZC_outputs.evaluator_set >= methods.hp_freq_ST);  % reject low freq ripple instances
    put_ripple_set{iChan} = all_ripple_locs_concat;
    ZC_outputs_set{iChan} = ZC_outputs;
end
save(sprintf('%s/%s_%d_polchecked_newRipp.mat',savedir,state,sleepID),'put_ripple_set','all_ripple_locs','ZC_outputs_set','chan_labels', 'N23_mask', 'sfreq','-v7.3')
disp('Artifact rejection complete')
toc 

% merge adjacent ripple locs; 
merged_idx_set = cell(1,num_chans);
if methods.adjacent_merge == 1
    disp('Beginning ripple merging.')
    for iChan = 1:num_chans
        cfg.artfctdef.visual.artifact = zeros(length(all_ripple_locs_concat),2);
        for i = 1:length(all_ripple_locs_concat)
            cfg.artfctdef.visual.artifact(i,1) = all_ripple_locs_concat(i)-methods.adjacency;
            cfg.artfctdef.visual.artifact(i,2) = all_ripple_locs_concat(i)+methods.adjacency;
        end
        [all_ripple_locs_concat,~,merged_idx] = CC_optimize_ripple_idx(data(iChan,:),data(iChan,:),data_filtered_set(iChan,:),cfg,ceil(0.1*sfreq)+1);  % magic number!

        while all_ripple_locs_concat(1) <= ceil(0.2*sfreq)
            all_ripple_locs_concat = all_ripple_locs_concat(2:end);
            %disp('pruning...')
        end
        while all_ripple_locs_concat(end)+ceil(0.5*sfreq) >= length(data)
            all_ripple_locs_concat(end) = [];
            %disp('...gninurp')
        end
        put_ripple_set{iChan} = all_ripple_locs_concat;
        merged_idx_set{iChan} = merged_idx;
    end
end

disp('Ripple merging complete.')
toc

good_ripple_locs_set = cell(1,num_chans);
bad_ripple_locs_set = cell(1,num_chans);
data_HP = [];
if methods.do_wavelet == 1
    disp('Beginning wavelet method.')
    for iChan = 1:num_chans
        all_ripple_locs_concat = put_ripple_set{iChan};
        tic;
        bad_bin_sum_all = zeros(1,length(all_ripple_locs_concat)); bad_bin_set = cell(1,length(num_chans_rej));
        delete(gcp('nocreate')); %%%%%%% NO PARFOR MOD
        %pool = parpool(methods.parpool_size); %%%%%%% NO PARFOR MOD
        for i = 1:num_chans_rej %%%%%%% NO PARFOR MOD - disabled for BW09
%     [~,bad_bin_set,~] = outlier_wavedecompC_dbtest2(data(i,:),all_ripple_locs_concat,methods);
%     bad_bin_sum = bad_bin_sum + bad_bin_set;
%     [~,bad_bin_set,~] = outlier_wavedecompC_dbtest2(data(i,:),all_ripple_locs_concat-sfreq,methods);
%     bad_bin_sum = bad_bin_sum + bad_bin_set;
%     [~,bad_bin_set,~] = outlier_wavedecompC_dbtest2(data(i,:),all_ripple_locs_concat+sfreq,methods);
%     bad_bin_sum = bad_bin_sum + bad_bin_set;
            bin_steps = [];
            for j = 1:length(all_ripple_locs_concat)
%         tic;
                bin_steps = [bin_steps,all_ripple_locs_concat(j)-round(0.1*sfreq):round(0.01*sfreq):all_ripple_locs_concat(j)+round(0.1*sfreq)];
%         [~,bad_bin_set,~] = outlier_wavedecompC_dbtest2(data(i,:),bin_steps,methods);
%         if sum(bad_bin_set) > 0
%             bad_bin_sum_all(j) = 1;
%         end
%         toc;
            end
            [~,bad_bin_set{i},~] = outlier_wavedecompC_dbtest2(data(i,:),bin_steps,methods);
        end
        %delete(pool) %%%%%%% NO PARFOR MOD

        bin_steps = [all_ripple_locs_concat(1)-round(0.1*sfreq):round(0.01*sfreq):all_ripple_locs_concat(1)+round(0.1*sfreq)];
        num_steps = length(bin_steps);
        for i = 1:num_chans_rej
            for j = 1:length(all_ripple_locs_concat)
                if sum(bad_bin_set{i}(num_steps*(j-1)+1:num_steps*j)) > 0
                    bad_bin_sum_all(j) = 1;
                end
            end
        end
    
        if methods.HP_reject == 1
            putative_bad_locs = all_ripple_locs_concat(bad_bin_sum_all > 0);
%             HM_locs = methods.HM_locs; HM_seg = methods.HM_seg;
%             HM_seg_HP = ft_preproc_highpassfilter(HM_seg,sfreq,200,4,'but');
%             HP_pool = zeros(1,length(HM_locs));
%         for iHM = 1:length(HM_locs)
%             HP_pool(iHM) = max(abs(HM_seg_HP(HM_locs(iHM)-round(0.25*sfreq):HM_locs(iHM)+round(0.25*sfreq))));
%         end
%         HP_thresh = mean(HP_pool) + methods.HP_thresh*std(HP_pool);
% %         HP_thresh = quantile(HP_pool,0.9);
% %         HP_thresh = mean(HP_pool);
            %data_HP(iChan,:) = ft_preproc_highpassfilter(data(iChan,:),sfreq,200,4,'but'); % charlie fix
            [b,a] = butter(3,200/(sfreq/2), 'high');
            data_HP(iChan,:) = filtfilt(b,a,data(iChan,:));
            for iPut = 1:length(putative_bad_locs)
%                 good_tally = 0;
%                 for iChan = 1:num_chans_rej
                    HP_eval = max(abs(data_HP(iChan,putative_bad_locs(iPut)-round(0.25*sfreq):putative_bad_locs(iPut)+round(0.25*sfreq))));
                if HP_eval < methods.HP_thresh
%                         good_tally = good_tally + 1;
% %                     break
%                     end
%                 end
%                 if good_tally == num_chans_rej
                    putative_bad_locs(iPut) = 0;
                end
            end
            putative_bad_locs(putative_bad_locs==0) = [];
            bad_ripple_locs_concat = putative_bad_locs;
            good_ripple_locs_concat = setdiff(all_ripple_locs_concat,bad_ripple_locs_concat);
        else
            good_ripple_locs_concat = all_ripple_locs_concat(bad_bin_sum_all == 0);
            bad_ripple_locs_concat = all_ripple_locs_concat(bad_bin_sum_all > 0);
        end
        good_ripple_locs_set{iChan} = good_ripple_locs_concat;
        bad_ripple_locs_set{iChan} = bad_ripple_locs_concat;
        toc;
    end
else
    disp('Skipping wavelet method.')
    good_ripple_locs_set = put_ripple_set;
    bad_ripple_locs_set = [];
end

toc
disp('Saving final results.')
save(sprintf('%s/%s_%d_polchecked_newRipp.mat',savedir,state,sleepID),'good_ripple_locs_set','bad_ripple_locs_set','merged_idx_set', 'chan_labels', 'N23_mask', 'methods', 'sfreq')
disp('Final results saved.')
% add. stuff here that separates SPW-Rs from NSPW-Rs; turned off for NC
% ripples since (as of May 2020) no need to shape-match like SWR
% the code below also hasn't been adjusted for different sets of ripples
% over multiple channels (for HC ripples were combined across electrode),
% and remain here only for reference
if methods.separation_flag == 1%     postRip_maxes_all = zeros(1,length(good_ripple_locs_concat));
%     for i = 1:length(good_ripple_locs_concat)
%         postRip_maxes_all(i) = max(data(methods.ZC_chan,good_ripple_locs_concat(i)+round(0.1*sfreq):good_ripple_locs_concat(i)+round(0.25*sfreq)) - data(good_ripple_locs_concat(i)));
%     end
%     good_ripple_locs_SPWR = good_ripple_locs_concat(postRip_maxes_all > methods.separation_thresh);  % biphasic
%     good_ripple_locs_NSPWR = good_ripple_locs_concat(postRip_maxes_all <= methods.separation_thresh);  % flat
    HM_template = methods.biphasic_template;
    HM_locs = methods.HM_locs;
    HM_evaluator = zeros(1,length(HM_locs));
    HM_seg = methods.HM_seg;
    for i = 1:length(HM_locs)
        temp_stretch = HM_locs(i)-round(0.01*sfreq):HM_locs(i)+round(0.01*sfreq);
        [~,temp_center] = min(HM_seg(temp_stretch));
        temp = HM_seg(temp_stretch(temp_center)-round(0.1*sfreq):temp_stretch(temp_center)+round(0.3*sfreq));
        HM_evaluator(i) = dot(temp,HM_template);
%         HM_evaluator(i) = pdist([temp;HM_template], 'cosine');
    end

    template_evaluator = zeros(1,length(good_ripple_locs_concat));
    for i = 1:length(good_ripple_locs_concat)
        temp_stretch = good_ripple_locs_concat(i)-round(0.01*sfreq):good_ripple_locs_concat(i)+round(0.01*sfreq);
        [~,temp_center] = min(data(methods.ZC_chan,temp_stretch));
        temp = data(methods.ZC_chan,temp_stretch(temp_center)-round(0.1*sfreq):temp_stretch(temp_center)+round(0.3*sfreq));
        if methods.rmDotBase == 1
            if good_ripple_locs_concat(i)-round(0.01*sfreq)-10*sfreq <= 0 
                baseline_stretch = 1:good_ripple_locs_concat(i)-round(0.01*sfreq)-1;
            else
                baseline_stretch = good_ripple_locs_concat(i)-round(0.01*sfreq)-10*sfreq:good_ripple_locs_concat(i)-round(0.01*sfreq);
            end
            baseline = mean(data(methods.ZC_chan,baseline_stretch));
            template_evaluator(i) = dot(temp-baseline,HM_template);
        else
            template_evaluator(i) = dot(temp,HM_template);
        end
    end
    temp_idx = template_evaluator > quantile(HM_evaluator,methods.separation_thresh_dot);
    good_ripple_locs_SPWR = good_ripple_locs_concat(temp_idx);
    flatness_evaluator = template_evaluator(temp_idx);
%     good_ripple_locs_NSPWR = good_ripple_locs_concat(~temp_idx);
    
    if methods.do_pmax == 1
        polarity = methods.hipp_polarity;
        postRip_maxes_all = zeros(1,length(good_ripple_locs_SPWR));
        for i = 1:length(good_ripple_locs_SPWR)
            if polarity == 0
                postRip_maxes_all(i) = max(data(methods.ZC_chan,good_ripple_locs_SPWR(i)+round(0.1*sfreq):good_ripple_locs_SPWR(i)+round(0.25*sfreq)) - data(methods.ZC_chan,good_ripple_locs_SPWR(i)));
            elseif polarity == 1
                postRip_maxes_all(i) = max(data(methods.ZC_chan,good_ripple_locs_SPWR(i)) - data(methods.ZC_chan,good_ripple_locs_SPWR(i)+round(0.1*sfreq):good_ripple_locs_SPWR(i)+round(0.25*sfreq)));
            end
        end

        postRip_maxes = zeros(1,length(HM_locs));
        for i = 1:length(HM_locs)
            if polarity == 0
                postRip_maxes(i) = max(HM_seg(HM_locs(i)+round(0.1*sfreq):HM_locs(i)+round(0.25*sfreq)) - HM_seg(HM_locs(i)));
            elseif polarity == 1
                postRip_maxes(i) = max(HM_seg(HM_locs(i)) - HM_seg(HM_locs(i)+round(0.1*sfreq):HM_locs(i)+round(0.25*sfreq)));
            end
        end
        evaluator_SPWR = flatness_evaluator(postRip_maxes_all > quantile(postRip_maxes,methods.separation_thresh_pmax));
%     evaluator_NSPWR = template_evaluator(~temp_idx);
        good_ripple_locs_SPWR = good_ripple_locs_SPWR(postRip_maxes_all > quantile(postRip_maxes,methods.separation_thresh_pmax));  % biphasic
    else
        evaluator_SPWR = flatness_evaluator;
    end
    [good_ripple_locs_NSPWR,idx] = setdiff(good_ripple_locs_concat,good_ripple_locs_SPWR);
    evaluator_NSPWR = template_evaluator(idx);

% re-scan for artifacts
    bad_ripple_locs_NSPWR = []; bad_ripple_locs_SPWR = [];
    if methods.do_wavelet == 1
        if ~isempty(good_ripple_locs_SPWR)
            tic;
            bad_bin_sum_all = zeros(1,length(good_ripple_locs_SPWR)); bad_bin_set = cell(1,num_chans_rej);
            for i = 1:num_chans_rej %%%%%%% NO PARFOR MOD
                bin_steps = [];
                for j = 1:length(good_ripple_locs_SPWR)
                    bin_steps = [bin_steps,good_ripple_locs_SPWR(j)-round(0.01*sfreq):1:good_ripple_locs_SPWR(j)+round(0.01*sfreq)];
                end
                [~,bad_bin_set{i},~] = outlier_wavedecompC_dbtest2(data_full(i,:),bin_steps,methods);
            end
            bin_steps = [good_ripple_locs_SPWR(1)-round(0.01*sfreq):1:good_ripple_locs_SPWR(1)+round(0.01*sfreq)];
            num_steps = length(bin_steps);
            for i = 1:num_chans_rej
                for j = 1:length(good_ripple_locs_SPWR)
                    if sum(bad_bin_set{i}(num_steps*(j-1)+1:num_steps*j)) > 0
                        bad_bin_sum_all(j) = 1;
                    end
                end
            end
            
            if methods.HP_reject == 1
                putative_bad_locs = good_ripple_locs_SPWR(bad_bin_sum_all > 0);
%                 HM_locs = methods.HM_locs; HM_seg = methods.HM_seg;
%                 HM_seg_HP = ft_preproc_lowpassfilter(HM_seg,sfreq,200,4,'but');
%                 HP_pool = zeros(1,length(HM_locs));
%                 for iHM = 1:length(HM_locs)
%                     HP_pool(iHM) = max(abs(HM_seg_HP(HM_locs(iHM)-round(0.1*sfreq):HM_locs(iHM)+round(0.1*sfreq))));
%                 end
%                 HP_thresh = mean(HP_pool) + 3*std(HP_pool);
%             data_HP = ft_preproc_lowpassfilter(data_full,sfreq,200,4,'but');
                for iPut = 1:length(putative_bad_locs)
                    good_tally = 0;
                    for iChan = 1:num_chans_rej
                        HP_eval = max(abs(data_HP(iChan,putative_bad_locs(iPut)-round(0.25*sfreq):putative_bad_locs(iPut)+round(0.25*sfreq))));
                        if HP_eval < methods.HP_thresh
                            good_tally = good_tally + 1;
%                             break
                        end
                    end
                    if good_tally == num_chans_rej
                        putative_bad_locs(iPut) = 0;
                    end
                end
                putative_bad_locs(putative_bad_locs==0) = [];
                bad_ripple_locs_SPWR = putative_bad_locs;
                [good_ripple_locs_SPWR,idx] = setdiff(good_ripple_locs_SPWR,bad_ripple_locs_SPWR);
                evaluator_SPWR = evaluator_SPWR(idx);
                
%                 for iPut = 1:length(good_ripple_locs_SPWR)
%                     for iChan = 1:num_chans_rej
%                         HP_eval = max(abs(data_HP(iChan,good_ripple_locs_SPWR(iPut)-round(0.25*sfreq):good_ripple_locs_SPWR(iPut)+round(0.25*sfreq))));
%                         if HP_eval > HP_thresh
%                             good_ripple_locs_SPWR(iPut) = 0;
%                             break
%                         end
%                     end
%                 end
%                 bad_ripple_locs_SPWR = sort([bad_ripple_locs_SPWR,good_ripple_locs_SPWR(good_ripple_locs_SPWR == 0)]);
%                 evaluator_SPWR(good_ripple_locs_SPWR == 0) = [];
%                 good_ripple_locs_SPWR(good_ripple_locs_SPWR == 0) = [];
            else
                bad_ripple_locs_SPWR = good_ripple_locs_SPWR(bad_bin_sum_all > 0);
                good_ripple_locs_SPWR = good_ripple_locs_SPWR(bad_bin_sum_all == 0); 
                evaluator_SPWR = evaluator_SPWR(bad_bin_sum_all == 0);
            end
            disp('SPWR rej time:')
            toc;
        end
    
        if ~isempty(good_ripple_locs_NSPWR)
            tic;
            bad_bin_sum_all = zeros(1,length(good_ripple_locs_NSPWR)); bad_bin_set = cell(1,num_chans_rej);
            for i = 1:num_chans_rej %%%%%%% NO PARFOR MOD
                bin_steps = [];
                for j = 1:length(good_ripple_locs_NSPWR)
                    bin_steps = [bin_steps,good_ripple_locs_NSPWR(j)-round(0.01*sfreq):1:good_ripple_locs_NSPWR(j)+round(0.01*sfreq)];
                end
                [~,bad_bin_set{i},~] = outlier_wavedecompC_dbtest2(data_full(i,:),bin_steps,methods);
            end
            for i = 1:num_chans_rej
                for j = 1:length(good_ripple_locs_NSPWR)
                    try
                        if sum(bad_bin_set{i}(num_steps*(j-1)+1:num_steps*j)) > 0
                            bad_bin_sum_all(j) = 1;
                        end
                    catch
                        disp([num2str(i),' ',num2str(j)]);
                    end
                end
            end
            
            if methods.HP_reject == 1
                putative_bad_locs = good_ripple_locs_NSPWR(bad_bin_sum_all > 0);
%                 HM_locs = methods.HM_locs; HM_seg = methods.HM_seg;
%                 HM_seg_HP = ft_preproc_lowpassfilter(HM_seg,sfreq,200,4,'but');
%                 HP_pool = zeros(1,length(HM_locs));
%                 for iHM = 1:length(HM_locs)
%                     HP_pool(iHM) = max(abs(HM_seg_HP(HM_locs(iHM)-round(0.1*sfreq):HM_locs(iHM)+round(0.1*sfreq))));
%                 end
%                 HP_thresh = mean(HP_pool) + 3*std(HP_pool);
%             data_HP = ft_preproc_lowpassfilter(data_full,sfreq,200,4,'but');
                for iPut = 1:length(putative_bad_locs)
                    good_tally = 0;
                    for iChan = 1:num_chans_rej
                        HP_eval = max(abs(data_HP(iChan,putative_bad_locs(iPut)-round(0.25*sfreq):putative_bad_locs(iPut)+round(0.25*sfreq))));
                        if HP_eval < methods.HP_thresh
                            good_tally = good_tally+1;
%                             break
                        end
                    end
                    if good_tally == num_chans_rej
                        putative_bad_locs(iPut) = 0;
                    end
                end
                putative_bad_locs(putative_bad_locs==0) = [];
                bad_ripple_locs_NSPWR = putative_bad_locs;
                [good_ripple_locs_NSPWR,idx] = setdiff(good_ripple_locs_NSPWR,bad_ripple_locs_NSPWR);
                evaluator_NSPWR = evaluator_NSPWR(idx);
            else
                bad_ripple_locs_NSPWR = good_ripple_locs_NSPWR(bad_bin_sum_all > 0);
                good_ripple_locs_NSPWR = good_ripple_locs_NSPWR(bad_bin_sum_all == 0); 
                evaluator_NSPWR = evaluator_NSPWR(bad_bin_sum_all == 0);
            end
            disp('NSPWR rej time:')
            toc;
        end
        delete(pool)
    end
    
    if methods.gap_thresh ~= 0.04
        removal_idx = zeros(1,length(good_ripple_locs_SPWR));
        for iSWR = 1:length(good_ripple_locs_SPWR)-1
            if removal_idx(iSWR) == 1
                continue
            end
            lead_SWR = good_ripple_locs_SPWR(iSWR);
            follow_SWR = good_ripple_locs_SPWR(iSWR+1);
            dist_SWR = follow_SWR - lead_SWR;
            if dist_SWR < methods.gap_thresh*sfreq
                lead_score = evaluator_SPWR(iSWR);
                follow_score = evaluator_SPWR(iSWR);
                if lead_score < follow_score
                    removal_idx(iSWR) = 1;
                else
                    removal_idx(iSWR+1) = 1;
                end
            end
        end
        good_ripple_locs_SPWR(logical(removal_idx)) = [];
        evaluator_SPWR(logical(removal_idx)) = [];
    end
    
    save(sprintf('%s/%s_%d_polchecked_newRipp.mat',savedir,state,sleepID),'good_ripple_locs_SPWR','good_ripple_locs_NSPWR','evaluator_SPWR','evaluator_NSPWR','bad_ripple_locs_SPWR','bad_ripple_locs_NSPWR', 'chan_labels', 'N23_mask', 'sfreq')
end
%% Plotting TF
% fig_dir = sprintf('%s/figs',savedir);
% if ~exist(fig_dir,'dir')
%     mkdir(fig_dir)
% end
% trial_dur = 0.4;  % in seconds
% window_size = round(trial_dur*sfreq);
% window_size_ctx = round(10*trial_dur*sfreq);
%     
%     identifiers.ctx_labels = chan_labels;
%     identifiers.fig_dir = fig_dir;
%     identifiers.subj = subj;
%     identifiers.sleepID = sleepID;
%     identifiers.stage = stage;
%     identifiers.hipp_parts_labels = hipp_parts_labels;
%     identifiers.ctx_labels = ctx_labels;
%     identifiers.hipp_parts = hipp_parts;
%     num_ctx_chans = size(data_ctx,1);
%     
%     if methods.separation_flag == 1
%         [trials_hipp_SPWR] = CC_ripples_TFplotter(num_chans,'short',window_size,good_ripple_locs_SPWR,data,identifiers);
%         [trials_hipp_SPWR_full] = CC_ripples_TFplotter(num_chans,'long',window_size_ctx,good_ripple_locs_SPWR,data,identifiers);
%         [trials_hipp_NSPWR] = CC_ripples_TFplotter(num_chans,'short',window_size,good_ripple_locs_NSPWR,data,identifiers);
%         [trials_hipp_NSPWR_full] = CC_ripples_TFplotter(num_chans,'long',window_size_ctx,good_ripple_locs_NSPWR,data,identifiers);
%         [trials_ctx_SPWR] = CC_ripples_TFplotter(num_ctx_chans,'long',window_size_ctx,good_ripple_locs_SPWR,data_ctx,identifiers);
%         [trials_ctx_NSPWR] = CC_ripples_TFplotter(num_ctx_chans,'long',window_size_ctx,good_ripple_locs_NSPWR,data_ctx,identifiers);
%         
%         save(sprintf('%s/sleep_%d_%s_trials_%s_newRipp_bisep.mat',savedir,sleepID,stage,hipp_parts{iHipp}),'trials_hipp_SPWR','trials_hipp_NSPWR','trials_hipp_SPWR_full','trials_hipp_NSPWR_full','-v7.3')
%         save(sprintf('%s/sleep_%d_%s_trials_%s_newRipp_bisep.mat',savedir,sleepID,stage,hipp_parts{iHipp}),'trials_ctx_SPWR','trials_ctx_NSPWR','-v7.3')
%     else
%         [trials_hipp_good] = CC_ripples_TFplotter(num_chans,'short',window_size,good_ripple_locs_auto_ST,data,identifiers);
%         [trials_hipp_good_full] = CC_ripples_TFplotter(num_chans,'long',window_size_ctx,good_ripple_locs_auto_ST,data,identifiers);
%         [trials_hipp_bad] = CC_ripples_TFplotter(num_chans,'short',window_size,bad_ripple_locs_auto_ST,data,identifiers);
%         [trials_hipp_bad_full] = CC_ripples_TFplotter(num_chans,'long',window_size_ctx,bad_ripple_locs_auto_ST,data,identifiers);
%         [trials_ctx_good] = CC_ripples_TFplotter(num_ctx_chans,'long',window_size_ctx,good_ripple_locs_auto_ST,data_ctx,identifiers);
%         [trials_ctx_bad] = CC_ripples_TFplotter(num_ctx_chans,'long',window_size_ctx,bad_ripple_locs_auto_ST,data_ctx,identifiers);
%         
%         save(sprintf('%s/sleep_%d_%s_trials_%s_newRipp.mat',savedir,sleepID,stage,hipp_parts{iHipp}),'trials_hipp_good','trials_hipp_bad','trials_hipp_good_full','trials_hipp_bad_full','-v7.3')
%         save(sprintf('%s/sleep_%d_%s_trials_%s_newRipp.mat',savedir,sleepID,stage,hipp_parts{iHipp}),'trials_ctx_good','trials_ctx_bad','-v7.3')
%     end
    
    disp(sprintf('Sleep %d ripples auto-evaluated!',sleepID))
    
end