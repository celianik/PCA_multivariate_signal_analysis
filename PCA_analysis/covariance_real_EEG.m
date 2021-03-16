load sampleEEGdata

% time window in which covariance will be computed
tidx = dsearchn(EEG.times',[0 800]');

% loop over trials and compute covariance for each trial
covmatT = zeros(EEG.nbchan);

for trial = 1:EEG.trials
  tmpdat = EEG.data(:,tidx(1):tidx(2),trial);
  % mean-center the data over time
  tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
  covmatT = covmatT + (tmpdat*tmpdat')/diff(tidx);
endfor
covmatT = covmatT/EEG.trials;

% cpvariance of trial anerage

erp = squeeze(mean(EEG.data(:,tidx(1):tidx(2),:),3));

tmpdat = bsxfun(@minus,erp,mean(erp,2));
covmatA = (tmpdat*tmpdat')/diff(tidx);

clim = [-1 1]*100;

figure(1), clf

subplot(1,2,1)
imagesc(covmatT), axis square
set(gca,'clim',clim)
xlabel('Channels'), ylabel('Channels')
title('Average of covariances')

subplot(1,2,2)
imagesc(covmatA), axis square
set(gca,'clim',clim/10)
xlabel('Channels'), ylabel('Channels')
title('Covariance of average')
