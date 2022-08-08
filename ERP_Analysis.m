%% Pre-Processing - Sample Script

% Load Behavioural Data
EEG.etc.behavioural_data = readtable('C:\sample_data\beh_data\oddball_sample_behavioural_data.xlsx', 'Sheet', 'oddball_sample_1.TRC');

%% Plotting Waveforms

chan2use = {'FCz'};
chan2use_idx = ismember(lower({EEG.chanlocs.labels}),lower(chan2use));

standard_idx = strcmpi(EEG.etc.behavioural_data.Standard_Target,'standard');
target_idx = strcmpi(EEG.etc.behavioural_data.Standard_Target,'target');

standard_erp = mean(EEG.data(chan2use_idx,:,standard_idx),3);
target_erp = mean(EEG.data(chan2use_idx,:,target_idx),3);

figure;
subplot(2,1,1);
hold on;
plot(EEG.times,standard_erp,'color',[0 0 0],'linewidth',1,'linestyle','-','DisplayName','Standard');
plot(EEG.times,target_erp,'color',[1 0 0],'linewidth',1,'linestyle','-','DisplayName','Target');
legend('location','best','autoupdate','off');
plot([0 0],[-10 10],'k');
plot([min(EEG.times) max(EEG.times)],[0 0],'k');
text(0.01,1,chan2use{:},'units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom');
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');
title('Individual Average Waveform');

% Plot Scalpmap
load roma_inv

time2use = [250 450];
time2use_idx = dsearchn(EEG.times',time2use');

standard_scalp = mean(EEG.data(:,time2use_idx(1):time2use_idx(2),standard_idx),[2 3]);
target_scalp = mean(EEG.data(:,time2use_idx(1):time2use_idx(2),target_idx),[2 3]);
diff_scalp = target_scalp - standard_scalp;

crange = [-5 5];

subplot(2,2,3);
hold on;
topoplot(standard_scalp,EEG.chanlocs,'electrodes','on');
colormap(roma_inv);
caxis(crange);
title('Standard');
subplot(2,2,4);
hold on;
topoplot(target_scalp,EEG.chanlocs,'electrodes','on');
colormap(roma_inv);
caxis(crange);
cb = colorbar;
set(get(cb,'label'),'string','\muV'); 
title('Target');