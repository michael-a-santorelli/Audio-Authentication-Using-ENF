%% Define global variables for ENF function
% Global variables for segmentation
global incr row col rows cols idx segmented_input 
% Global variables for windowing
global window windowed 
% Global variables for zero padding
global zeropadding zero_padded 
% Global variables for obtaining frequency data
global frequency_domain freq_per_div segmented_freq frequency_range X Y seg_freq_range
% Global variables for max and weighted energy
global max_en weight_en

%% Experiment 1: Obtain ENF for regular audio signal
% Close all windows to avoid too many open windows if any from prior codes
close all
% ENF function call shown below
% enf( x, Fs, Blocksize, Zeropad, Overlap, Window, Frequency )
% Load audio data
[data, fs] = audioread('recording.wav');
% Obtain ENF at 60 Hz
[max_60,weight_60] = enf_function(data, fs, fs*16, 0, 0.5, 'hanning', 60);
% Obtain ENF at 120 Hz
[max_120,weight_120] = enf_function(data, fs, fs*16, 0, 0.5, 'hanning', 120);
% Obtain ENF at 180 Hz
[max_180,weight_180] = enf_function(data, fs, fs*16, 0, 0.5, 'hanning', 180);
% Obtain ENF at 240 Hz
[max_240,weight_240] = enf_function(data, fs, fs*16, 0, 0.5, 'hanning', 240);

%% Experiment 2: Obtain ENF for ground truth, Comparison of Audio Signal to Ground Truth without Preprocessing
% This section of the script must be run after the previous section as it
% requires "weight_120" which is obtained in that section
% Close all windows to avoid too many open windows if any from prior codes
close all
% ENF function call shown below
% enf( x, Fs, Blocksize, Zeropad, Overlap, Window, Frequency )
% Load ground truth data
[data_gt, fs] = audioread('ground truth.wav');
% Obtain ENF for ground truth data at 120 Hz
[max_gt_120,weight_gt_120] = enf_function(data_gt, fs, fs*16, 0, 0.5, 'hanning', 120);
% Calculate mean of weighted energy for data at 120 Hz and ground truth at
% 120 Hz
m1 = mean(weight_120); m2 = mean(weight_gt_120);
% Calculate cross correlation between weighted energy of data at 120 Hz and
% ground truth at 120 Hz
corr = xcorr(weight_gt_120-m2, weight_120-m1);
% Plot cross correlation between weighted energy of data at 120 Hz and
% ground truth at 120 Hz
figure
plot(corr)
title('Cross Correlation between Weighted Energy of Audio Signal & Ground Truth')
% Determine the index which the highest cross correlation occurs; zero pad
% with delays
[~, delays] = find(ismember(corr, max(corr(:)))); delay_blocks = size(frequency_domain,2) - delays; delay_zeros = zeros(1,abs(delay_blocks));
% Pad weighted energy of audio data at 120 Hz to match delay of ground
% truth
new_weight_120 = horzcat(delay_zeros,weight_120);
% Scale for display purposes
for i=1:length(new_weight_120)
    if(new_weight_120(i) > 0)
        new_weight_120(i) = new_weight_120(i) - m1;
    end
end
for i=1:length(weight_gt_120)
    weight_gt_120(i) = weight_gt_120(i) - m2;
end
% Plot newly padded weighted energy of audio signal (normalized) and plot
% of ground truth (normalized)
figure
plot(new_weight_120)
hold on
plot(weight_gt_120)
hold off
legend('Audio Signal', 'Ground Truth')
title('Comparison of Weighted Energy between Audio Signal & Ground Truth')

%% Experiment 3: Comparison of Audio Signal and Ground Truth with Preprocessing (LFP, downsampling)
% Close all windows to avoid too many open windows if any from prior codes
close all
% ENF function call shown below
% enf( x, Fs, Blocksize, Zeropad, Overlap, Window, Frequency )
% Load audio signal data and ground truth data
[data,fs] = audioread('recording.wav');
[data_gt, fs] = audioread('ground truth.wav');
% Downsample sampling frequency by a factor of 100
fs = fs/100;
% G and SOS are the filter parameters obtained from designing the LPF
% required from the specifications using MATLAB's FDATOOL
G = [0.000163246154108335;0.000161899106682263;0.000160817321811931;0.000160119000399767;0.0126442830025082;1];
SOS = [1,2,1,1,-1.99049274831148,0.991145732927911;1,2,1,1,-1.97406793176450,0.974715528191229;1,2,1,1,-1.96087751419301,0.961520783480255;1,2,1,1,-1.95236274265371,0.953003218655308;1,1,0,1,-0.974711433994984,0];
% Apply filter to audio data
filtered = filtfilt(SOS,G,data);
decimated = [1:100:length(filtered)];
% Apply filter to ground truth data
filteredgt = filtfilt(SOS,G,data_gt);
decimatedgt = [1:100:length(filteredgt)];
% Obtain downsampled filtered data and ENF
new_data = filtered(decimated);
[max_new_data,weight_new_data] = enf_function(new_data, fs, fs*16, 16384-(fs*16), 0.5, 'hanning', 120);
% Obtain downsampled filtered data and ENF
new_datagt = filteredgt(decimatedgt);
[max_new_datagt,weight_new_datagt] = enf_function(new_datagt, fs, fs*16, 16384-(fs*16), 0.5, 'hanning', 120);
% Calculate mean of weighted energies
m1 = mean(weight_new_data); m2 = mean(weight_new_datagt);
% Obtain cross correlation between audio signal and ground truth after
% downsample and filter
corr = xcorr(weight_new_datagt-m2, weight_new_data-m1);
% Plot cross correlation
figure
plot(corr)
title('Cross Correlation between Weighted Energy of Audio Signal & Ground Truth')
% Obtain delay and create zeros equivalent to delay value for padding 
[~, delays] = find(ismember(corr, max(corr(:)))); delay_blocks = size(frequency_domain,2) - delays; delay_zeros = zeros(1,abs(delay_blocks));
% Apply delay padding to audio signal
new_weight_120 = horzcat(delay_zeros,weight_new_data);
% Scale audio signal weighted energy and ground truth weighted energy
for i=1:length(new_weight_120)
    if(new_weight_120(i) > 0)
        new_weight_120(i) = new_weight_120(i) - m1;
    end
end
for i=1:length(weight_new_datagt)
    weight_new_datagt(i) = weight_new_datagt(i) - m2;
end
% Plot comparison between audio signal and ground truth
figure
plot(new_weight_120)
hold on
plot(weight_new_datagt)
hold off
legend('Audio Signal', 'Ground Truth')
title('Comparison of Weighted Energy between Audio Signal & Ground Truth')

%% Experiment 4: LPF, downsampling, and zeropadding, second set
% Close all windows to avoid too many open windows if any from prior codes
close all
% ENF function call shown below
% enf( x, Fs, Blocksize, Zeropad, Overlap, Window, Frequency )
% Load audio signal data and ground truth data
[data,fs] = audioread('recording 2.wav');
[data_gt, fs] = audioread('ground truth 2.wav');
% Downsample sampling frequency by a factor of 100
fs = fs/100;
% G and SOS are the filter parameters obtained from designing the LPF
% required from the specifications using MATLAB's FDATOOL
G = [0.000163246154108335;0.000161899106682263;0.000160817321811931;0.000160119000399767;0.0126442830025082;1];
SOS = [1,2,1,1,-1.99049274831148,0.991145732927911;1,2,1,1,-1.97406793176450,0.974715528191229;1,2,1,1,-1.96087751419301,0.961520783480255;1,2,1,1,-1.95236274265371,0.953003218655308;1,1,0,1,-0.974711433994984,0];
% Apply filter to audio data
filtered = filtfilt(SOS,G,data);
decimated = [1:100:length(filtered)];
% Apply filter to ground truth data
filteredgt = filtfilt(SOS,G,data_gt);
decimatedgt = [1:100:length(filteredgt)];
% Obtain downsampled filtered data and ENF
new_data = filtered(decimated);
[max_new_data,weight_new_data] = enf_function(new_data, fs, fs*16, 16384-(fs*16), 0.5, 'hanning', 120);
% Obtain downsampled filtered data and ENF
new_datagt = filteredgt(decimatedgt);
[max_new_datagt,weight_new_datagt] = enf_function(new_datagt, fs, fs*16, 16384-(fs*16), 0.5, 'hanning', 120);
% Calculate mean of weighted energies
m1 = mean(weight_new_data); m2 = mean(weight_new_datagt);
% Obtain cross correlation between audio signal and ground truth after
% downsample and filter
corr = xcorr(weight_new_datagt-m2, weight_new_data-m1);
% Plot cross correlation
figure
plot(corr)
title('Cross Correlation between Weighted Energy of Audio Signal & Ground Truth')
% Obtain delay and create zeros equivalent to delay value for padding
[~, delays] = find(ismember(corr, max(corr(:)))); delay_blocks = size(frequency_domain,2) - delays; delay_zeros = zeros(1,abs(delay_blocks));
% Apply delay padding to audio signal
new_weight_120 = horzcat(delay_zeros,weight_new_data);
% Scale audio signal weighted energy and ground truth weighted energy
for i=1:length(new_weight_120)
    if(new_weight_120(i) > 0)
        new_weight_120(i) = new_weight_120(i) - m1;
    end
end
for i=1:length(weight_new_datagt)
    weight_new_datagt(i) = weight_new_datagt(i) - m2;
end
% Plot comparison between audio signal and ground truth
figure
plot(new_weight_120)
hold on
plot(weight_new_datagt)
hold off
legend('Audio Signal', 'Ground Truth')
title('Comparison of Weighted Energy between Audio Signal & Ground Truth')


