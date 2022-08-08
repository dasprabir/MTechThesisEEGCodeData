%% 
% 
%   Connectivity hubs
% 
% 


clear
%load sampleEEGdata.mat
load T01.mat
%load N6_1.5s_pre.mat
EEG.data = double(EEG.data);

% frequency in hz
frex = 10;

% time window for synchronization
tidx = dsearchn(EEG.times',[ 0 500 ]');


% time vector for wavelet
wtime = -1:1/EEG.srate:1;
fwhm = .2;


% convolution parameters
nData = EEG.pnts*EEG.trials;
nWave = length(wtime);
nConv = nData + nWave - 1;
halfW = floor(nWave/2);


% create wavelet
cmwX = fft( exp(1i*2*pi*frex*wtime) .* exp( -4*log(2)*wtime.^2 / fwhm^2 ) ,nConv );

dataX = fft( reshape(EEG.data,EEG.nbchan,[]) ,nConv,2 );


% convolution
as = ifft( bsxfun(@times,dataX,cmwX) ,[],2);
as = as(:,halfW:end-halfW-1);
as = reshape( as,size(EEG.data) );

% get angles
allphases = angle(as);

%% compute all-to-all PLI

pliall = zeros(EEG.nbchan);

for chani=1:EEG.nbchan
    for chanj=chani+1:EEG.nbchan
        
        % Euler-format phase differences
        cdd = exp( 1i*(allphases(chani,tidx(1):tidx(2),:)-allphases(chanj,tidx(1):tidx(2),:)) );
        cdd = squeeze(cdd);
        
        % compute PLI for this channel pair
        plitmp = mean( abs(mean(sign(imag(cdd)),1)) ,2);
        
        % enter into matrix!
        pliall(chani,chanj) = plitmp;
        pliall(chanj,chani) = plitmp;
    end
end
        
% let's see what it looks like
figure(10), clf
imagesc(pliall)
axis square
xlabel('Channels'), ylabel('Channels')
title([ 'All-to-all connectivity at ' num2str(frex) ' Hz' ])
set(gca,'clim',[0.1 0.7])
colormap jet
colorbar

%% now for hubness

% define a threshold
% gather unique data values into a vector (convenience)
distdata = nonzeros(triu(pliall));

% define a threshold
%thresh = mean(distdata) + std(distdata);
thresh = 0.48;

% threshold the matrix!
pliallThresh = pliall>thresh;


% plots!
figure(11), clf
subplot(311), hold on
histogram(distdata,50)
plot([1 1]*thresh,get(gca,'ylim'),'r--','linew',3)
xlabel('PLI (synch. strength)'), ylabel('Count')
legend({'Distribution';'Threshold'})


subplot(3,2,[3 5])
imagesc(pliallThresh)
axis square
xlabel('Channels'), ylabel('Channels')
title([ 'All-to-all connectivity at ' num2str(frex) ' Hz' ])


subplot(3,2,[4 6])
topoplotIndie(sum(pliallThresh)/EEG.nbchan,EEG.chanlocs,'numcontour',0);
set(gca,'clim',[.1 .4])
title('Topoplot of "hubness"')
colormap hot
colorbar

%%% QUESTION: Does the threshold affect the qualitative topographical distribution?
% 
% 
%%% QUESTION: Do the results look different for different frequencies or time windows?
% 
% 


%% done.
