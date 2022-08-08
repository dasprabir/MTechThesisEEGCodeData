%%

%    Connectivity
%    Phase synchronization matrices

%
%%

% load dataset
%load v1_laminar.mat
load N6_1.5s_pre.mat
csd = double(EEG.data);
srate = EEG.srate;
% filter parameters (in hz)
centfreq = 10;
freqwidt = 2;

% time window parameters (in seconds)
timewin{1} = [0 500];
timewin{2} = [600 900];
% convert to indices
timevec = EEG.times;
tidx1 = dsearchn(timevec',timewin{1}');
tidx2 = dsearchn(timevec',timewin{2}');

%% filter and extract angles

% apply filter
filtdat = filterFGx(csd,srate,centfreq,freqwidt,1);

% extract angles via hilbert transform
angst = zeros(size(filtdat));
for triali=1:size(csd,3)
    angst(:,:,triali) = angle(hilbert(squeeze(filtdat(:,:,triali))').');
end

%% synchronization matrices in two time windows

nchans = size(csd,1);

% initialize
synchmat = zeros(2,nchans,nchans);

for chani=1:nchans
    for chanj=1:nchans
        
        %%% time window 1
        % extract angles
        tmpAi = angst(chani,tidx1(1):tidx1(2),:);
        tmpAj = angst(chanj,tidx1(1):tidx1(2),:);
        % synch on each trial
        trialsynch = abs(mean(exp(1i*( tmpAi-tmpAj )),2)); % averaging over time
        % average over trials
        synchmat(1,chani,chanj) = mean(trialsynch);
        
        
        
        %%% time window 2
        % extract angles
        tmpAi = angst(chani,tidx2(1):tidx2(2),:);
        tmpAj = angst(chanj,tidx2(1):tidx2(2),:);
        % synch on each trial
        trialsynch = abs(mean(exp(1i*( tmpAi-tmpAj )),2));
        % average over trials
        synchmat(2,chani,chanj) = mean(trialsynch);
    end
end


%% show the results


figure(1), clf
clim = [0 1.2];

% each map individually
for i=1:2
    subplot(1,3,i)
    imagesc(squeeze(synchmat(i,:,:)))
    set(gca,'clim',clim)
    axis square
    xlabel('Channel'), ylabel('Channel')
    title([ 'Synch: ' num2str(timewin{i}(1)) '-' num2str(timewin{i}(2)) 'ms, ' num2str(centfreq) ' Hz' ])
    colorbar
end

% show difference
subplot(133)
imagesc(squeeze(diff(synchmat)))
set(gca,'clim',[-1 1]*mean(clim)/2)
axis square
xlabel('Channel'), ylabel('Channel')
title('Difference: late - early')
colorbar
colormap jet
%% done.
