%% 
% 
%   Spectral analysis of resting-state EEG
% 
% 

% The goal of this cell is to plot a power spectrum of resting-state EEG data.

clear
%load EEGrestingState.mat
load N6_singlechandata.mat
srate = 256;
% create a time vector that starts from 0
npnts = length(singlechandata);
time  = (0:npnts-1)/srate;


% plot the time-domain signal 
figure(13), clf
plot(time,singlechandata)
xlabel('Time (s)'), ylabel('Voltage (\muV)') % alpha range in posterior channel in the back of the head during resting state
zoom on


% static spectral analysis
hz = linspace(0,srate/2,floor(npnts/2)+1);
ampl = 2*abs(fft(singlechandata)/npnts);
powr = ampl.^(2);

figure(14), clf, hold on
plot(hz,ampl(1:length(hz)),'k','linew',2)
plot(hz,powr(1:length(hz)),'r','linew',2)

xlabel('Frequency (Hz)')
ylabel('Amplitude or power')
legend({'Amplitude';'Power'})

% optional zooming in
%set(gca,'xlim',[0 30])


%%% QUESTION: What are the three most prominent features of the EEG spectrum?

% (1/f)(power is decreasing with increasing freq.)
% line noise spike
% Present of rythmic activity and oscillatory activity

%%% QUESTION: What do you notice about the difference between the amplitude
%             and power spectra?
% 
%%% QUESTION: Can you see the ~10 Hz oscillation in the raw time series data?
% 
