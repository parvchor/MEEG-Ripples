addpath('/home/bqrosen/matlab/fieldtrip-20210825/') 
ft_defaults;

subj      = 'CA140';
bids_root = '/space/seh10/2/halgdev/projects/pchordiya/Cogigate/BIDS MEEG Data/COG_MEEG_EXP1_BIDS_SAMPLE'; % Ensure this path is correct

fif_file  = fullfile(bids_root, sprintf('sub-%s/ses-1/meg/sub-%s_ses-1_task-dur_run-01_meg.fif', subj, subj));

cfg = [];
cfg.dataset = fif_file;         
cfg.channel = 'EEG';
data_raw = ft_preprocessing(cfg);

cfg = [];
cfg.resamplefs   = 500;    
cfg.detrend      = 'no';   
cfg.demean       = 'no';   
data_resampled = ft_resampledata(cfg, data_raw);

cfg = [];
cfg.method    = 'mtmfft';        
cfg.output    = 'pow';           
cfg.taper     = 'hanning';       
cfg.foi       = 1:1:100;         
data_freq = ft_freqanalysis(cfg, data_resampled);

avg_pow = mean(data_freq.powspctrm, 1);   
pow_db  = 10 * log10(avg_pow);

figure;
plot(data_freq.freq, pow_db, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title(sprintf('PSD (1–100 Hz) for sub-%s (all channels averaged)', subj));
xlim([1 100]);
grid on;