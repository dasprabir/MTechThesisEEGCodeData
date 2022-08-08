%% Compute and plot TF-ITPC for one electrode

%load sampleEEGdata.mat
load T01.mat

% wavelet parameters
num_frex = 40;
min_freq =  2;
max_freq = 30;

channel2use = 'pz';

% set range for variable number of wavelet cycles
range_cycles = [ 3 10 ];

% parameters (notice using logarithmically spaced frequencies!)
frex  = logspace(log10(min_freq),log10(max_freq),num_frex);
nCycs = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex);
time  = -2:1/EEG.srate:2;
half_wave = (length(time)-1)/2;

% FFT parameters
nWave = length(time);
nData = EEG.pnts*EEG.trials;
nConv = nWave+nData-1;


% FFT of data (doesn't change on frequency iteration)
dataX = fft( reshape(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);

% initialize output time-frequency data
tf = zeros(num_frex,EEG.pnts);

% loop over frequencies
for fi=1:num_frex
    
    % create wavelet and get its FFT
    s = nCycs(fi)/(2*pi*frex(fi));
    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
    waveletX = fft(wavelet,nConv);
    
    % question: is this next line necessary?
    waveletX = waveletX./max(waveletX);
    
    % run convolution
    as = ifft(waveletX.*dataX,nConv);
    as = as(half_wave+1:end-half_wave);
    as = reshape(as,EEG.pnts,EEG.trials);
    
    % compute ITPC
    tf(fi,:) = abs(mean(exp(1i*angle(as)),2));
end

% plot results
figure(12), clf
contourf(EEG.times,frex,tf,40,'linecolor','none')
set(gca,'clim',[0 .45],'ydir','normal','xlim',[-300 1000])
title('ITPC')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
colormap jet