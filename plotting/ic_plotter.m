addpath('/home/bqrosen/matlab/fieldtrip-20210825/') 
ft_defaults;

subj      = 'CB013';
bids_root = '/space/seh10/2/halgdev/projects/pchordiya/Cogigate/BIDS MEEG Data/COG_MEEG_EXP1_BIDS_SAMPLE'; % Ensure this path is correct

fif_file  = fullfile(bids_root, sprintf('sub-%s/ses-1/meg/sub-%s_ses-1_task-dur_run-01_meg.fif', subj, subj));

cfg            = [];
cfg.dataset    = fif_file;        
cfg.channel    = 'EEG';           
data_raw       = ft_preprocessing(cfg);

cfg            = [];
cfg.resamplefs = 500;    % new sampling rate in Hz
cfg.detrend    = 'no';   % no detrending before resampling
cfg.demean     = 'no';   % no demeaning here
data_resampled = ft_resampledata(cfg, data_raw);

cfg             = [];
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [70 100];    % bandpass between 70 and 100 Hz
cfg.bpfilttype  = 'but';     % windowed‐sinc FIR (or 'fir' / 'but' as needed)
cfg.bpfiltord   = 3;           % filter order, adjust if necessary
cfg.channel     = ft_channelselection('EEG', data_resampled.label); 
data_filt       = ft_preprocessing(cfg, data_resampled);

eeg_chan       = ft_channelselection('EEG', data_filt.label);
fprintf('Selected %d EEG channels for ICA.\n', numel(eeg_chan));

cfg              = [];
cfg.method       = 'runica';       
components       = ft_componentanalysis(cfg, data_filt);

cfg             = [];
cfg.elec        = data_resampled.elec;       % electrode structure with .elecpos & .label
layout          = ft_prepare_layout(cfg);

%% 7) Plot the scalp topographies for the first 10 ICs
cfg             = [];
cfg.component   = 1:10;              % plot IC 1 through IC 10
cfg.layout      = layout;            % use our custom layout
cfg.comment     = 'no';              % suppress overlay text (“IC 01”)
figure;
ft_topoplotIC(cfg, components);
sgtitle('ICA Component Topographies (IC 1–10)');

%% 8) Plot the time courses (waveforms) of the first 10 ICs
ics             = components.trial{1};               % [nICs × nTimePoints]
fs              = data_resampled.fsample;      % sampling frequency (e.g., 250 Hz)
nSamples        = size(ics, 2);
tvec            = (0:(nSamples-1)) / fs;       % time axis in seconds

figure;
for icIdx = 1:10
    subplot(5, 2, icIdx);
    plot(tvec, ics(icIdx, :), 'LineWidth', 1);
    xlim([tvec(1) tvec(end)]);
    xlabel('Time (s)');
    ylabel('Amplitude (a.u.)');
    title(sprintf('IC %d Time Course', icIdx));
    grid on;
end
sgtitle('Time Courses of Top 10 ICA Components');

%% 10) Compute PSD for each IC and plot in a 5×2 grid (IC 1–IC 10)
%   a) Construct a pseudo-“raw” data structure for IC channels so that 
%      ft_freqanalysis can run on it directly.
data_ic            = [];
data_ic.label      = components.label;     % {'EEG001', 'EEG002', …} → but now each is an IC
data_ic.fsample    = fs;                 % e.g., 250 Hz
data_ic.trial      = { ics };            % a single cell: [nIC × nTimePoints]
data_ic.time       = { tvec };           % a single cell: [1 × nTimePoints]
% We do NOT need elec or grad here for PSD.

%   b) Compute PSD (1–100 Hz in 1 Hz steps, Hanning taper)
cfg              = [];
cfg.method       = 'mtmfft';
cfg.output       = 'pow';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:140;
data_freq_ic     = ft_freqanalysis(cfg, data_ic);
% → data_freq_ic.powspctrm is [nIC × nFreq], data_freq_ic.freq is [1×nFreq]

%   c) Plot each IC’s PSD in a 5×2 subplot grid
figure('Name','IC PSDs (IC1–IC10)');
for icIdx = 1:10
    subplot(5,2,icIdx);
    plot(data_freq_ic.freq, squeeze(data_freq_ic.powspctrm(icIdx,:)), 'LineWidth', 1);
    xlim([60 140]);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(sprintf('IC %d PSD', icIdx));
    grid on;
end
sgtitle('Power Spectral Density of ICA Components 1–10');


ic_data           = [];
ic_data.label     = components.label;      % same as data_ic.label
ic_data.fsample   = fs;
ic_data.trial     = { ics };             % [nIC × nTimePoints]
ic_data.time      = { tvec };            % [1 × nTimePoints]
%   b) Open ft_databrowser in “vertical” mode to scroll through each IC’s time series:
cfg              = [];
cfg.viewmode     = 'vertical';           % stack ICs vertically
cfg.ylim         = 'maxmin';             % let FT choose y‐limits per channel
cfg.showlabels   = 'yes';                % show component name (e.g., “EEG001” → now “IC01”, etc.)
cfg.blocksize    = 10;                   % seconds of data per screen (adjust as you like)
figure('Name','IC Time Series Browser');
ft_databrowser(cfg, ic_data);