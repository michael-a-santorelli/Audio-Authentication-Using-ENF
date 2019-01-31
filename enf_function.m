function [ y1, y2 ] = enf( x, Fs, Blocksize, Zeropad, Overlap, Window, Frequency )

%enf: This function takes in 7 parameters: data, sampling frequency,
%     blocksize, zero padding length, overlap, type of window, and
%     frequency.
%
%     This function produces 2 outputs: y1 is the maximum energy vector of
%     the input data, and y2 is the weighted energy of the input data.
%
%     First, this function computes a segmentation of the input data using
%     the blocksize entered and the overlap entered. Then, it creates a
%     window using the specified window, and applies it to the segmented
%     input signal. Next, it applies zero padding to windowed input
%     using the specified zeropad length. Then, it transforms the zero
%     padded signal into the frequency domain using FFT. Next, it segments
%     the frequency domain signal to a range of +/- 1 Hz from the specified
%     frequency parameter, and calculates the ENF of the frequency domain 
%     signal and produces a surface plot of this. Lastly, it calculates and
%     plots the maximum and weighted energy of the segmented frequency
%     domain signal.

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

%% Input Segmentation
% Calculate the incrementing value for creating the columns matrix
incr = ceil(Blocksize*(1-Overlap));
% Define column matrix
col = [0:incr:length(x)-Blocksize]';
cols = repmat(col, [1 Blocksize]);
% Define row matrix
row = [1:1:Blocksize];
rows = repmat(row, [length(col) 1]);
% Create indexing matrix by combining the rows and columns matrix
idx = rows+cols;
% Segment input by indexing through using the above created matrix
segmented_input = x(idx);

%% Windowing
% Choose window to create based on entered "Window" parameter
if(Window == 'hamming')
    window = hamming(Blocksize, 'periodic');
elseif(Window == 'hanning')
    window = hanning(Blocksize, 'periodic');
elseif(Window == 'barlett')
    window = bartlett(Blocksize, 'periodic');
elseif(Window == 'none')
    window = ones(Blocksize);
end
% Change window from a column vector to a row vector, and repeat it using
% repmat
window = repmat(window, [1 length(col)]);
window = window';
% Apply window to segmented input, and make output a column matrix
windowed = (segmented_input.*window);

%% Zeropadding
% Create zeros matrix for zero padding
zeropadding = zeros(length(col), Zeropad);
% Apply zeropadding to windowed input
zero_padded = (horzcat(windowed, zeropadding))';

%% FFT
% Obtain frequency domain magnitude representation of zero padded input
frequency_domain = abs(fft(zero_padded));

%% Frequency of Interest
% Detemine frequency per division of the input data
freq_per_div = Fs/(Blocksize+length(zeropadding));
% Define frequency vector specifying the frequency at each row of the data
frequency_range = [0:freq_per_div:Fs-freq_per_div];
frequency_range = frequency_range';
% Segment the frequency domain data by extracting +/- 1 Hz from the desired
% frequency
segmented_freq = frequency_domain(find(frequency_range >= (Frequency-1-freq_per_div)& frequency_range <= (Frequency+1-freq_per_div)),:);
% Obtain segmented frequency range values in terms of Hz
seg_freq_range = (find(frequency_range >= Frequency-1 - freq_per_div & frequency_range <= Frequency+1 - freq_per_div))*freq_per_div;
seg_freq_range = repmat(seg_freq_range, [1, length(col)]);
% Define X and Y to be the time divison and frequency division of the
% segmented data in terms of its original time and frequency components
X = [1:1:length(col)];
Y = [Frequency-1:freq_per_div:Frequency+1]';
if(Zeropad > 0)
    Y = [Frequency-1:freq_per_div:Frequency+1-freq_per_div]';
end
% Plot the ENF of the input signal
figure
surf(X,Y,segmented_freq);
xlabel('Time Block');
ylabel('Frequency (Hz)')
zlabel('Magnitude')
if(Frequency == 60)
    title('Electric Network Frequency: 60 Hz')
elseif(Frequency == 120)
    title('Electric Network Frequency: 120 Hz')
elseif(Frequency == 180)
    title('Electric Network Frequency: 180 Hz')
elseif(Frequency == 240)
    title('Electric Network Frequency: 240 Hz')
end

%% Max and Weighted Energy
% Calculate maximum energy in each time division
[~, max_en] = max(segmented_freq);
max_en = Frequency-1 + max_en*freq_per_div;
% Plot the maximum energy spectrum of each time division
figure
plot(X,max_en)
xlabel('Time Block')
ylabel('Frequency')
if(Frequency == 60)
    title('Maximum Energy: 60 Hz')
elseif(Frequency == 120)
    title('Maximum Energy: 120 Hz')
elseif(Frequency == 180)
    title('Maximum Energy: 180 Hz')
elseif(Frequency == 240)
    title('Maximum Energy: 240 Hz')
end
% Calculate the weighted energy in each time division
weight_en = sum(segmented_freq.*seg_freq_range)./sum(segmented_freq);
% Plot the weighted energy spectrum of each time division
figure
plot(X,weight_en)
xlabel('Time Block')
ylabel('Frequency')
if(Frequency == 60)
    title('Weighted Energy: 60 Hz')
elseif(Frequency == 120)
    title('Weighted Energy: 120 Hz')
elseif(Frequency == 180)
    title('Weighted Energy: 180 Hz')
elseif(Frequency == 240)
    title('Weighted Energy: 240 Hz')
end

%% Assign outputs
% Assign outputs: y1 is the maximum energy spectrum and y2 is the weighted
% energy spectrum
y1 = max_en;
y2 = weight_en;


