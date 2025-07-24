function CortRipple_wrapper_MGH(subject, cond)

%%
delete(gcp('nocreate'));

%subject = 'UAB02'; % 02, 06, 08
%subjects = {'UAB01','UAB03','UAB05','UAB09','UAB14','UAB17','UAB29','UAB48'};
% subject = 'UAB26';
% cond = 'sleep';

% select epochs with minimal putative IIS, high frequency periods, or otherwise large artifacts
do_epoch_select = 1;

% exclude channels marked as bad?
remove_bad_chan = 0;

sfreq = 1000;

output_dir = 'NC_ripple';
savedir = sprintf('%s/%s/%s',output_dir,subject,cond);
dir_exist(savedir);

disp(subject)

dat.raw = [];

% load channel labels
% load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/labels/%s_labels.mat', subject));
% chan_labels = labels(:,1);

% load bipolar labels, require that channels included match those in labels above 
%load(sprintf('/space/seh10/3/halgdev/analysis/MULTI/segment_1kHz/%s/%s/%s_segment_%03d.mat',subject(1:3), subject, subject, 1), 'bpLabel');
load(sprintf('/space/seh10/3/halgdev/analysis/MULTI/segment_1kHz_singlePrecision/%s/%s/%s_segment_%03d.mat',subject(1:2), subject, subject, 1), 'bpLabel');
% gMsk = contains(bpLabel, chan_labels); % assumes order is maintained
       
clear n_points;
     
% load segment files
%flist = dir(sprintf('/space/seh10/3/halgdev/analysis/MULTI/segment_1kHz_singlePrecision/UAB/%s/', subject));
flist = dir(sprintf('/space/seh10/3/halgdev/analysis/MULTI/segment_1kHz_singlePrecision/MG/%s/', subject));

flist = {flist.name};
flist = flist(contains(flist, 'segment'));
last_seg = str2num(flist{end}(end-6:end-4));

% load(sprintf('/space/seh8/1/halgdev/projects/skajfez/SleepScoring/%s/%s/%s_hypnogram.mat',subject(1:3), subject, subject), 'hyp')
load(sprintf('/space/seh8/1/halgdev/projects/skajfez/SleepScoring/%s/%s/%s_N23.mat',subject(1:2), subject, subject), 'N23')
load(sprintf('/space/seh8/1/halgdev/projects/skajfez/ChannelSelection/%s/%s/%s_gMsk.mat',subject(1:2), subject, subject), 'gMsk')

binSec = 30; % pre-processed data are 30s binned
if strcmp(cond, 'sleep')
%     mask = repelem(hyp==2 | hyp==3, binSec*sfreq);  
    mask = repelem(N23, binSec*sfreq);
elseif strcmp(cond, 'waking')
    % waking periods selected by Sophie
%     load(sprintf('/space/seh8/1/halgdev/projects/skajfez/SleepScoring/UAB/%s/%s_wake.mat', subject,subject))
%     mask = repelem(wake, binSec*sfreq);
end

f = waitbar(0,'loading all SEEG segments...');
for seg = 1:last_seg 
    %load(sprintf('/space/seh10/3/halgdev/analysis/MULTI/segment_1kHz_singlePrecision/UAB/%s/%s_segment_%03d.mat', subject, subject, seg), 'data', 'bpOp', 'bpLabel', 'clockT');
    load(sprintf('/space/seh10/3/halgdev/analysis/MULTI/segment_1kHz_singlePrecision/MG/%s/%s_segment_%03d.mat', subject, subject, seg), 'data', 'bpOp', 'bpLabel', 'clockT');
    
    if ~exist('n_points', 'var')
        n_points_init = size(data,1); % almost but not totally robust
    end
    n_points = size(data,1);

    % select N23 epochs
    onset = (seg-1)*n_points_init+1; % i.e. n_points of last file
    offset = onset+n_points-1;
    if offset > numel(mask) % correct for laset file
        offset = numel(mask);
    end
    
    if size(mask,1) > size(mask,2)
        mask = mask';
    end

    dat_tmp = double(data)*bpOp;
    
    if strcmp(cond, 'sleep')
        % select between 9 PM and 6 AM
        clockT_perSecond = clockT(1:sfreq:end);
        timeVals = datestr(clockT_perSecond);
        timeVals = str2num(timeVals(:,end-7:end-6));
        timeVals = ismember(timeVals,[21:23 0:6])';
        
        timeVals = repelem(timeVals,sfreq);
    else
        timeVals = true(size(mask(onset:offset)));
    end
    
    if sum(mask(onset:offset))
        if numel(timeVals)>numel(mask(onset:offset))
            timeVals = timeVals(1:numel(mask(onset:offset)));
        end
        mask(onset:offset) = logical(mask(onset:offset) .* timeVals);
        
        N23_seg = mask(onset:offset);
        dat.raw = [dat.raw dat_tmp(N23_seg,gMsk)'];
    end
    clear data dat_tmp
    
    waitbar(seg/last_seg,f,'loading all SEEG segments...')

end
close all

% find mask edges
mask_edges = [0 diff(mask)];
mask_edges = abs(mask_edges(mask));

if size(dat.raw,2) ~= numel(mask_edges)
    error('mask and data sizes do not match')
end
        
chan_labels = bpLabel(gMsk);

% remove bad channels
if remove_bad_chan
    load(sprintf('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/bad_chan/%s_bad_chan.mat', subject));
    bad_chan_ind = ismember(chan_labels, bad_chan);
    chan_labels(bad_chan_ind) = [];
    dat.raw(bad_chan_ind,:) = [];
end
    

% check the raw data here before moving - e.g. plot(dat.raw(1,:))
% f2 = waitbar(0,'running epoch selection...');
if do_epoch_select
    
    % parameters for spikey regions to remove
    binSec = 30; % size of bins
    binSz0 = binSec;
    spikePad = 0.3; % padding around spikes
    spikeBinThresh = 0.5; % if spiking (with padding) this prop of time
    nBin = numel(N23); % number of 30 s bins
    %nBin = numel(hyp); % number of 30 s bins
    
    % very large artifacts to remove
    max_abs_thresh = 750; % microvolts per ms
    bigArt_win = 1; % +/- seconds around bounds of large artifact
    
    % find bins with many spikes
    load(sprintf('/space/seh10/3/halgdev/analysis/MULTI/preproc/%s/%s/%s_preproc.mat', subject(1:2), subject, subject), 'spike_regs')
    spike_regs = spike_regs(gMsk);
    spike_regs = cellfun(@(x) x*(sfreq/500), spike_regs, 'un', 0);
    spikeMsk0 = bounds2mask(cat(1,spike_regs{:}),nBin*binSec*sfreq,round(sfreq*spikePad));    
    spikeMsk = reshape(spikeMsk0(1:nBin*binSec*sfreq),binSec*sfreq,nBin);
    spikeMsk = sum(spikeMsk,1) > (binSz0*sfreq*spikeBinThresh)';
    spikeMsk = reshape(repelem(spikeMsk,1,binSec*sfreq)',1,[]);
    art_mask = spikeMsk(mask);
   
    % FIND LARGE ARTIFACTS
    RP = regionprops(max(abs(dat.raw))>max_abs_thresh, 'SubArrayIdx', 'area');
    RPidx = find([RP.Area]>0);
    art_center = NaN(numel(RPidx),1);
    for r_ind = 1:numel(RPidx)
        r = RPidx(r_ind);
        art_onset = RP(r).SubarrayIdx{2}(1);
        art_offset = RP(r).SubarrayIdx{2}(end);
        if art_onset < bigArt_win*sfreq
            art_mask(1:art_offset+bigArt_win*sfreq) = 1;
        elseif art_offset > size(dat.raw,2)-bigArt_win*sfreq
            art_mask(art_offset-bigArt_win*sfreq:end) = 1;
        else
            art_mask(art_onset-bigArt_win*sfreq:art_offset+bigArt_win*sfreq) = 1;
        end
    end

    % mark as artifactual if the moving RMS exceeds this value
    rms_thresh = 2;
    rms_window = 10; % minutes
    rms_overlap = 5; % minutes
    
    % moving RMS 10 minutes with 5 minute overlap
    dat.rms = NaN(size(dat.raw));
    for ch = 1:size(dat.raw,1)
        tmp = repelem(rms_mov(dat.raw(ch,:),rms_window*60*sfreq,rms_overlap*60*sfreq,1),rms_overlap*60*sfreq);
        if numel(tmp) > size(dat.raw,2)
            tmp(size(dat.raw,2)+1:end) = [];
        end
        dat.rms(ch,:) = tmp;
    end

    art_add = mean(dat.rms,1) > (prctile(mean(dat.rms,1),10)*rms_thresh); % 1.5 * 10th percentile
    art_mask = art_mask | art_add;
    
    % close small gaps
    RP = regionprops(not(art_mask), 'SubArrayIdx', 'area');
    RPidx = find([RP.Area]<(15*sfreq));
    for r_ind = 1:numel(RPidx)
        r = RPidx(r_ind);
        art_mask(RP(r).SubarrayIdx{2}(1):RP(r).SubarrayIdx{2}(end)) = 1;
    end
    
    mask_edges = mask_edges(~art_mask);
    mask_edges = find(mask_edges);
    art_mask_edges = [0 diff(art_mask)];
    art_mask_edges = art_mask_edges(~art_mask);
    art_mask_edges = find(abs(art_mask_edges));
    all_mask_edges = unique([mask_edges art_mask_edges]);
    
    nan_pad = 50; % samples
    nan_edge_mask = false(sum(~art_mask),1);
    for e = 1:numel(all_mask_edges)
        if all_mask_edges(e) < nan_pad
            nan_edge_mask(1:all_mask_edges(e)) = 1;
        elseif all_mask_edges(e) < (size(dat.raw,2)-nan_pad)
            nan_edge_mask(all_mask_edges(e)-nan_pad:all_mask_edges(e)+nan_pad) = 1;
        end
    end
    
    if numel(nan_edge_mask) ~= sum(~art_mask)
        error('nan_edge_mask does not match data size')
    end
    
    
    % remove bad epochs
    dat.raw(:,art_mask) = [];

end


load('/space/seh8/1/halgdev/projects/xjiang3/Ripples_newSubj/rippleDetectZC_HM/default_methods.mat')
methods.separation_flag = 0;
methods.adjacent_merge = 1;
methods.adjacency = 25; % in sampling points
% broader range for freq center investigation
methods.hp_freq_ST = 60;  % in Hz
methods.lp_freq_ST = 120;
methods.looseEval = 1;
methods.gap_thresh = 0.1;   % set as 100 ms for now; might want to turn this off for NC
methods.do_pmax = 0;   % flag off since no hand-marked template for cortical ripples (there's no "sharp-wave" equivalent as far as we know)
methods.do_wavelet = 0;  % next 7 params are to tweak outlier_wavedecompC_dbtest2.m
methods.time_bin = 2;
methods.wavthresh_Haar = 6;
methods.wavthresh_RL1 = 6;
methods.wavthresh_RL2 = 8;
methods.wavthresh_RL3 = 8;
methods.wavthresh_RL4 = 5;
methods.wavthresh_RL5 = 4;
methods.HP_reject = 1;
methods.HP_thresh = 3;   % reject based on high gamma power;
 
methods.parpool_size = 16;
% 
if exist('HM_thresh','var')
    methods.separation_thresh_pmax = HM_thresh;
    clear HM_thresh
end

if exist('HM_polarity','var')
    methods.hipp_polarity = HM_polarity;
    clear HM_polarity
end

% % notch filter data
% dat.raw = notch(dat.raw,sfreq);
% size(dat.raw)

fprintf('\n%f hours\n', size(dat.raw,2)/1000/60/60)

% identify HC, TH, and AMY channels
AMY_chan = contains(chan_labels, labels(contains(labels(:,10), 'Amygdala'),1));
HC_chan = contains(chan_labels, labels(contains(labels(:,10), 'Hippocampus'),1));
TH_chan = contains(chan_labels, labels(contains(labels(:,10), 'Thalamus'),1));
NC_chan = ~(AMY_chan | HC_chan | TH_chan);
clear dat.raw 

fprintf('\n%d NC\n', sum(NC_chan));
fprintf('\n%d TH\n', sum(TH_chan));
fprintf('\n%d HC\n', sum(HC_chan));
fprintf('\n%d AMY\n', sum(AMY_chan));

% save epoch selection information
if do_epoch_select
    save(sprintf('%s/epochSelection.mat',savedir),'spikeMsk', 'binSec', 'binSz0', 'spikePad', 'spikeBinThresh', 'nBin', 'max_abs_thresh', 'bigArt_win', 'art_mask')
    save(sprintf('%s/nan_edge_mask.mat',savedir),'nan_edge_mask', 'all_mask_edges')
end

% detect putative ripples
if any(NC_chan)
    dat.raw_NC = dat.raw(NC_chan,:);
    chan_labels_NC = chan_labels(NC_chan);
    save_dir_NC = [savedir '/NC'];
    dir_exist(save_dir_NC);
    rippleDetect_NCmod(dat.raw_NC,sfreq,save_dir_NC,1,[],chan_labels_NC,cond,methods, mask, cond);
end

if any(HC_chan)
    dat.raw_HC = dat.raw(HC_chan,:);
    chan_labels_HC = chan_labels(HC_chan);
    save_dir_HC = [savedir '/HC'];
    dir_exist(save_dir_HC);
    rippleDetect_NCmod(dat.raw_HC,sfreq,save_dir_HC,1,[],chan_labels_HC,cond,methods, mask, cond);
end

if any(TH_chan)
    dat.raw_TH = dat.raw(TH_chan,:);
    chan_labels_TH = chan_labels(TH_chan);
    save_dir_TH = [savedir '/TH'];
    dir_exist(save_dir_TH);
    rippleDetect_NCmod(dat.raw_TH,sfreq,save_dir_TH,1,[],chan_labels_TH,cond,methods, mask, cond);
end

if any(AMY_chan)
    dat.raw_AMY = dat.raw(AMY_chan,:);
    chan_labels_AMY = chan_labels(AMY_chan);
    save_dir_AMY = [savedir '/AMY'];
    dir_exist(save_dir_AMY);
    rippleDetect_NCmod(dat.raw_AMY,sfreq,save_dir_AMY,1,[],chan_labels_AMY,cond,methods, mask, cond);
end

rippleDetect_NCmod_IAV(dat.raw,sfreq,savedir,1,[],chan_labels,cond,methods, mask, cond);

end