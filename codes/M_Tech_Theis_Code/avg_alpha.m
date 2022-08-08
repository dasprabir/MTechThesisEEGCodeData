%% 
% 
%   Quantify alpha (8-12)hz power over the scalp
% 
% 

clear
load N6_1.5s_pre.mat

% These data comprise 63 "epochs" of resting-state data. Each epoch is a
%   2-second interval cut from ~2 minutes of resting-state.
% The goal is to compute the power spectrum of each 2-second epoch
%   separately, then average together.
%   Then, extract the average power from 8-12 Hz (the "alpha band") and
%   make a topographical map of the distribution of this power.

% convert to double-precision
EEG.data = double(EEG.data);

chanpowr = (2*abs(fft(EEG.data,[],2)/EEG.pnts) ).^2;

% then average over trials
chanpowr = mean(chanpowr,3);

% vector of frequencies
hz = linspace(0,EEG.srate/2,floor(EEG.pnts/2)+1);


% do some plotting
% plot power spectrum of all channels
figure(15), clf
plot(hz,chanpowr(:,1:length(hz)),'linew',2)
xlabel('Frequency (Hz)'), ylabel('Power (\muV)')

set(gca,'xlim',[0 30],'ylim',[0 80])

%% now to extract alpha power

% boundaries in hz
alphabounds = [8 12];

% convert to indices
freqidx = dsearchn(hz',alphabounds');


% extract average power
alphapower = mean(chanpowr(:,freqidx(1):freqidx(2)),2);

% and plot
figure(16), clf
topoplotIndie(alphapower,EEG.chanlocs,'numcontour',0);
set(gca,'clim',[0 6]*4)
colormap hot