% Ploting erp on ANN data
clear all
close all
clc
%load T01.mat
%load N6_1s_epochs.mat
load N6_1.5s_pre.mat
EEG

erp = mean(EEG.data,3);
chan2plot = 'cz';

figure(1),clf
plot(EEG.times,erp(strcmpi({EEG.chanlocs.labels},chan2plot),:),'linew',2)
xlabel('Time (ms)'), ylabel('Activity (\muV)')
set(gca,'xlim',[-400 1200])

%Topoplot
% convert time in ms to time in indices
time2plot=300;
[~,tidx] = min(abs(EEG.times-time2plot));

figure(2), clf
topoplotIndie(erp(:,tidx),EEG.chanlocs);
title([ 'ERP from ' num2str(time2plot) ' ms' ])
colormap jet
colorbar
