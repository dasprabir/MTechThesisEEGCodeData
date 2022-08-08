%% 1) All-to-all connectivity matrix
%  In this exercise, you are going to compute connectivity between
%    all pairs of channels. It will create a channel X channel matrix.

% First, import the v1 data and setup parameters for wavelet convolution 
%   at 8 hz and 20 hz.

clear

% load in data
%load v1_laminar
load T01.mat
%load N6_1.5s_pre.mat
%load N6_attrem_1s_128.mat

csd = EEG.data;
timevec = EEG.times;
% useful variables for later...
[nchans, npnts, ntrials] = size(EEG.data);

srate = EEG.srate;

% specify frequencies
frex = [8 20];
nCycles = [ 7 14 ];

% parameters for complex Morlet wavelets
wavtime  = -1:1/srate:1-1/srate; % why remove a sample point?!
half_wav = (length(wavtime)-1)/2;

% FFT parameters
nWave = length(wavtime);
nData = npnts*ntrials;
nConv = nWave+nData-1;

% and create wavelets
cmwX = zeros(length(frex),nConv);
for fi=1:length(frex)
    s       = nCycles(fi) / (2*pi*frex(fi));
    cmw      = exp(1i*2*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) ./ (2*s^2) );
    tempX     = fft(cmw,nConv);
    cmwX(fi,:) = tempX ./ max(tempX);
end

%% run convolution to extract phase values (don't need power)
% store the phase angle time series in a 
%    channels X frequency X time X trials matrix

allphases = zeros(nchans,length(frex),npnts,ntrials);

% spectrum of all channels using the fft matrix input 
% (check the matrix sizes and FFT inputs!)
dataX = fft( reshape(EEG.data,nchans,npnts*ntrials) ,nConv,2);

for fi=1:length(frex)
    
    % run convolution
    as = ifft( bsxfun(@times,dataX,cmwX(fi,:)) ,nConv,2 );
    as = as(:,half_wav+1:end-half_wav);
    as = reshape(as,size(csd));
    
    % phase values from all trials
    allphases(:,fi,:,:) = angle( as );
end


%% now compute connectivity
% Compute connectivity separately in two time windows:
%   .1 to .4, and .6 to .9 seconds.
% To do this, first compute synchronization over trials, then average the
%    synchronization values within those time windows.

% define time windows
tidx1 = dsearchn(timevec',[-1000 0]');
tidx2 = dsearchn(timevec',[254 1028]');


% initialize a channels X channels X frequency X time period matrix
connmat = zeros(nchans,nchans,length(frex),2);

% in a double for-loop, compute phase synchronization between each pair
% inside the for-loop 
for chani=1:nchans
    for chanj=1:nchans
        
        % compute eulerized phase angle differences
        phasediffs = exp(1i* squeeze( allphases(chani,:,:,:)-allphases(chanj,:,:,:) ) );
        
        % compute phase synchronization (ISPC) for all time points
        ispc = abs(mean(phasediffs,3));
        
        % get data averaged from the two time windows
        connmat(chani,chanj,:,1) = mean( ispc(:,tidx1(1):tidx1(2)) ,2);
        connmat(chani,chanj,:,2) = mean( ispc(:,tidx2(1):tidx2(2)) ,2);
    end
end

%% Make all-to-all connectivity plots
% In one figure, make six chan-by-chan matrices for 8 and 55 hz (upper and lower plots)
%   from averaged connectivity between .1-.4 s (left) and .6-.9 s (middle). 
% The right-most plot should show the difference of late-early connectivity.
% Use the same colorscaling for all 'raw' plots, 
%   and a different colorscaling for the difference plots.

% define color limits
clim  = [0 .8];
climD = [-.4 .4];


figure(1), clf
subplot(231)
imagesc(squeeze(connmat(:,:,1,1)))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ num2str(round(timevec(tidx1(1)),1)) '-' num2str(round(timevec(tidx1(2)),1)) 's, ' num2str(frex(1)) ' Hz' ])


subplot(232)
imagesc(squeeze(connmat(:,:,1,2)))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ num2str(round(timevec(tidx2(1)),1)) '-' num2str(round(timevec(tidx2(2)),1)) 's, ' num2str(frex(1)) ' Hz' ])


subplot(233)
imagesc(squeeze(connmat(:,:,1,2)-connmat(:,:,1,1)))
axis square
set(gca,'clim',climD,'xtick',1:nchans,'ytick',1:nchans)
title([ 'post-pre, ' num2str(frex(1)) ' Hz' ])






subplot(234)
imagesc(squeeze(connmat(:,:,2,1)))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ num2str(round(timevec(tidx1(1)),1)) '-' num2str(round(timevec(tidx1(2)),1)) 's, ' num2str(frex(2)) ' Hz' ])


subplot(235)
imagesc(squeeze(connmat(:,:,2,2)))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ num2str(round(timevec(tidx2(1)),1)) '-' num2str(round(timevec(tidx2(2)),1)) 's, ' num2str(frex(2)) ' Hz' ])

subplot(236)
imagesc(squeeze(connmat(:,:,2,2)-connmat(:,:,2,1)))
axis square
set(gca,'clim',climD,'xtick',1:nchans,'ytick',1:nchans)
title([ 'post-pre, ' num2str(frex(2)) ' Hz' ])

% QUESTIONS:
%   How do you interpret any of these things?!?!
