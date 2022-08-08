%% load EEG data and extract reaction times in ms

clear
%load sampleEEGdata.mat
load T01.mat

%%% note about the code in this cell:
%   this code extracts the reaction time from each trial
%   in the EEGLAB data format. You don't need to worry about
%   understanding this code if you do not use EEGLAB.

rts = zeros(size(EEG.epoch));

% loop over trials
for ei=1:EEG.trials
    
    % find the index corresponding to time=0, i.e., trial onset
    [~,zeroloc] = min(abs( cell2mat(EEG.epoch(ei).eventlatency) ));
    
    % reaction time is the event after the trial onset
    rts(ei) = EEG.epoch(ei).eventlatency{zeroloc+1};
end


% always good to inspect data, check for outliers, etc.
figure(17), clf
plot(rts,'ks-','markerfacecolor','w','markersize',12)
xlabel('Trial'), ylabel('Reaction time (ms)')

%% Create the design matrix
%  Our design matrix will have two regressors (two columns): intercept and RTs

X = [ ones(EEG.trials,1) rts' ];

%% Run wavelet convolution for time-frequency analysis
%  We didn't cover this in class, but this code extracts a time-frequency
%  map of power for each trial. These power values become the dependent
%  variables.

freqrange  = [2 25]; % extract only these frequencies (in Hz)
numfrex    = 30;     % number of frequencies between lowest and highest


% set up convolution parameters
wavtime = -2:1/EEG.srate:2;
frex    = linspace(freqrange(1),freqrange(2),numfrex);
nData   = EEG.pnts*EEG.trials;
nKern   = length(wavtime);
nConv   = nData + nKern - 1;
halfwav = (length(wavtime)-1)/2;
nCyc    = logspace(log10(4),log10(12),numfrex);

% initialize time-frequency output matrix
tf3d = zeros(numfrex,EEG.pnts,EEG.trials);

% compute Fourier coefficients of EEG data (doesn't change over frequency!)
eegX = fft( reshape(EEG.data(47,:,:),1,[]) ,nConv);

% loop over frequencies
for fi=1:numfrex
    
    %%% create the wavelet
    s    = nCyc(fi) / (2*pi*frex(fi));
    cmw  = exp(2*1i*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) / (2*s.^2) );
    cmwX = fft(cmw,nConv);
    cmwX = cmwX ./ max(cmwX);
    
    % second and third steps of convolution
    as = ifft( eegX .* cmwX );
    
    % cut wavelet back to size of data
    as = as(halfwav+1:end-halfwav);
    as = reshape(as,EEG.pnts,EEG.trials);
    
    % extract power from all trials
    tf3d(fi,:,:) = abs(as).^2;
    
end % end frequency loop

%% inspect the TF plots a bit

figure(18), clf

% show the raw power maps for three trials
for i=1:3
    subplot(2,3,i)
    imagesc(EEG.times,frex,squeeze(tf3d(:,:,i)))
    axis square, axis xy
    set(gca,'clim',[0 5],'xlim',[-200 1200])
    xlabel('Time (ms)'), ylabel('Frequency')
    title([ 'Trial ' num2str(i) ])
end


% now show the trial-average power map
subplot(212)
imagesc(EEG.times,frex,squeeze(mean(tf3d,3)))
axis square, axis xy
set(gca,'clim',[0 5],'xlim',[-200 1200])
xlabel('Time (ms)'), ylabel('Frequency')
title('All trials')

%% now for the regression model

% We're going to take a short-cut here, and reshape the 3D matrix to 2D.
% That doesn't change the values, and we don't alter the trial order.
% Note the size of the matrix below.
tf2d = reshape(tf3d,numfrex*EEG.pnts,EEG.trials)';


% Now we can fit the model on the 2D matrix
b = (X'*X)\X'*tf2d;

% reshape b into a time-by-frequency matrix
betamat = reshape(b(2,:),numfrex,EEG.pnts);

%% show the design and data matrices

figure(19), clf

ax1_h = axes;
set(ax1_h,'Position',[.05 .1 .1 .8])
imagesc(X)
set(ax1_h,'xtick',1:2,'xticklabel',{'Int';'RTs'},'ydir','norm')
ylabel('Trials')
title('Design matrix')


ax2_h = axes;
set(ax2_h,'Position',[.25 .1 .7 .8])
imagesc(tf2d)
set(ax2_h,'ydir','norm','clim',[0 20])
ylabel('Trials')
xlabel('Timefrequency')
title('Data matrix')

colormap gray


%%% QUESTION: Please interpret the matrices! 
%             What do they mean and what do they show?
% 
% 

%% show the results

figure(20), clf

% show time-frequency map of regressors
contourf(EEG.times,frex,betamat,40,'linecolor','none')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'xlim',[-200 1200],'clim',[-.012 .012]/4)
title('Regression against RT over trials')

%%% QUESTION: How do you interpret the results?
% 
% 