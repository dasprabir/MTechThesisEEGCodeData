%% 2) Seeded synchronization with topographical maps.
%     The goal of this assignment is to explore "seeded" synchronization,
%     which means synchronization from one electrode to all other
%     electrodes.

% load in the sample EEG dataset
clear
%load sampleEEGdata
%load T01.mat
load N6_1.5s_pre.mat

% laplacian for later
EEG.lap = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);

%% Perform wavelet convolution on all channels/trials. Save the single-trial phase data for all channels.
%  Don't worry about the power data.


% soft-coded parameters for time-frequency analysis
freqrange  = [2 40];
numfrex    = 50;
numcycles  = linspace(3,9,numfrex);


% set up wavelet and convolution parameters
wavtime = -2:1/EEG.srate:2;
frex    = linspace(freqrange(1),freqrange(2),numfrex);
nData   = EEG.pnts*EEG.trials;
nKern   = length(wavtime);
nConv   = nData + nKern -1;
halfwav = (length(wavtime)-1)/2;


% create wavelets
waveletsX = zeros(numfrex,nConv);
for fi=1:numfrex
    
    % create time-domain wavelet
    s = numcycles(fi) / (2*pi*frex(fi));
    wavelet = exp(1i*2*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) / (2*s^2) );
    
    % compute fourier coefficients of wavelet. need to normalize?
    waveletsX(fi,:) = fft(wavelet,nConv);
end



% initialize matrix to store all phase values
allphases = zeros(EEG.nbchan,numfrex,EEG.pnts,EEG.trials);


% now do convolution!
% loop over channels
for chani=1:EEG.nbchan
    
    % compute Fourier coefficients of EEG data (doesn't change over frequency!)
    eegX = fft( reshape(EEG.data(chani,:,:),1,nData) ,nConv);
    
    % loop over frequencies
    for fi=1:numfrex
        
        % second and third steps of convolution
        as = ifft( eegX.*waveletsX(fi,:) );
        % reshape to time X trials
        as = reshape(as(halfwav+1:end-halfwav),EEG.pnts,EEG.trials);
        
        % save phase values from all time points and trials
        allphases(chani,fi,:,:) = angle( as );
        
    end % end frequency loop
end % end channel loop

%% 

% Compute FCz-to-all and Oz-to-all connectivity in the window of 0-600 ms
% from 8-12 Hz. Show topoplots.
% Then compute connectivity from -600 to 0. Subtract to show task-related
% connectivity. Interpret in context of volume conduction.


% time indices
tidx1 = dsearchn(EEG.times',[-600 0]');
tidx2 = dsearchn(EEG.times',[0 600]');

% frequency indices
fidx = dsearchn(frex',[8 12]');

% channel indices
FCzidx = strcmpi('fcz',{EEG.chanlocs.labels});
Ozidx  = strcmpi('oz',{EEG.chanlocs.labels});


% initialize seeded synchronization maps
[fcz_seeded_synch,oz_seeded_synch] = deal( zeros(2,EEG.nbchan) );


% loop over channels
for chani=1:EEG.nbchan
    
    % compute complex unit vectors defined by phase angle differences
    % from the 'seed' channel to the looping channel
    fcz_euler = exp(1i* (allphases(chani,:,:,:)-allphases(FCzidx,:,:,:)) );
    oz_euler  = exp(1i* (allphases(chani,:,:,:)-allphases(Ozidx,:,:,:)) );
    
    % compute ITPC seeded synchronization over all time-frequency points
    synch_fcz_all = abs(mean(fcz_euler,4));
    synch_oz_all  = abs(mean(oz_euler,4));
    
    % compute PLI seeded synchronization over all time-frequency points (used later)
%     synch_fcz_all = abs(mean(sign(imag(fcz_euler)),4));
%     synch_oz_all  = abs(mean(sign(imag(oz_euler)),4));
    
    
    
    % extract average synchronization from the time-frequency windows
    fcz_seeded_synch(1,chani) = mean(mean(synch_fcz_all(1,fidx(1):fidx(2),tidx1(1):tidx1(2)),2),3);
    oz_seeded_synch(1,chani)  = mean(mean(synch_oz_all(1,fidx(1):fidx(2),tidx1(1):tidx1(2)),2),3);
    
    fcz_seeded_synch(2,chani) = mean(mean(synch_fcz_all(1,fidx(1):fidx(2),tidx2(1):tidx2(2)),2),3);
    oz_seeded_synch(2,chani)  = mean(mean(synch_oz_all(1,fidx(1):fidx(2),tidx2(1):tidx2(2)),2),3);
    
end

%% now plotting

% use the same color limit
clim = [0 1]/2; % for the 'raw' synchronization
climdiff = [-.1 .1]; % for the difference in synchronization



figure(3), clf

subplot(321)
topoplotIndie(fcz_seeded_synch(1,:),EEG.chanlocs,'numcontour',0);
set(gca,'clim',[0 10])
title('FCz seed, pre')

subplot(323)
topoplotIndie(fcz_seeded_synch(2,:),EEG.chanlocs,'numcontour',0);
set(gca,'clim',clim)
title('FCz seed, post')

subplot(325)
topoplotIndie(diff(fcz_seeded_synch),EEG.chanlocs,'numcontour',0);
set(gca,'clim',climdiff)
title('FCz seed, post-pre')



subplot(322)
topoplotIndie(oz_seeded_synch(1,:),EEG.chanlocs,'numcontour',0);
set(gca,'clim',clim)
title('Oz seed, pre')

subplot(324)
topoplotIndie(oz_seeded_synch(2,:),EEG.chanlocs,'numcontour',0);
set(gca,'clim',clim)
title('Oz seed, post')

subplot(326)
topoplotIndie(diff(oz_seeded_synch),EEG.chanlocs,'numcontour',0);
set(gca,'clim',climdiff)
title('Oz seed, post-pre')

%% 3) Adapt the code above to use PLI instead of ISPC.
%     Plot in a different figure so you can compare them. 
%     What do you think of the results?


%% 4) Adapt the code above again using the Laplacian instead of voltage.
%     Comment!


%% 5) Time-frequency maps of inter-regional synchronization.
%     Make two time-frequency maps of FCz-Oz ISPC; one plotting the 'raw' connectivity
%     and one plotting the baseline-subtracted connectivity (use a baseline of -500 to -200 ms).

% baseline time window index
baseidx = dsearchn(EEG.times',[-500 -200]');



% extract the connectivity
phasediffs = exp(1i* (allphases(FCzidx,:,:,:)-allphases(Ozidx,:,:,:)) );
synch_raw  = squeeze(abs(mean(phasediffs,4)));
synch_base = bsxfun(@minus,synch_raw,mean(synch_raw(:,baseidx(1):baseidx(2)),2));



% plotting 
figure(5), clf
subplot(211)
contourf(EEG.times,frex,synch_raw,40,'linecolor','none')
set(gca,'clim',[0 .3],'xlim',[-300 1200])
xlabel('Time (ms)'), ylabel('Frequencies (Hz)')
title('FCz-Oz ISPC, "raw" phase synch.')
colorbar

subplot(212)
contourf(EEG.times,frex,synch_base,40,'linecolor','none')
set(gca,'clim',[-.3 .3],'xlim',[-300 1200])
xlabel('Time (ms)'), ylabel('Frequencies (Hz)')
title('FCz-Oz ISPC, baseline-subtracted')
colorbar


%% done.
