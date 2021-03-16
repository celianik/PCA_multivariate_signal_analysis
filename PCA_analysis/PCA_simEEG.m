%% simulate data
load emptyEEG

% index of dipole to simulate activity in
diploc = 109;   

figure(1), clf

subplot(1,2,1)
plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), 'bo','markerfacecolor','y')
hold on
plot3(lf.GridLoc(diploc,1), lf.GridLoc(diploc,2), lf.GridLoc(diploc,3), 'rs','markerfacecolor','k','markersize',10)
rotate3d on, axis square
title('Brain dipole locations')

subplot(1,2,2)
topoplotIndie(-lf.Gain(:,1,diploc), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
title('Signal dipole projection')

% time points and time vector
N = 1000;
EEG.times = (0:N-1)/EEG.srate;

% create random data in all brain dipoles
dipole_data = randn(N,length(lf.Gain));

% add signal to half of dataset
dipole_data(:,diploc) = 15*sin(2*pi*10*EEG.times);

% project all data from dipoles to scalp electrodes
EEG.data = ( dipole_data*squeeze(lf.Gain(:,1,:))' )';

%% compute PCA

% Following the 6-Step procedure

% Step 1: compute covariance matrix of the data
% First I need to mean-center the data over time
data2 = bsxfun(@minus,EEG.data,mean(EEG.data,2));
covmat = data2*data2'/size(EEG.data,2);

% Step 2: perform eigendecomposition
[evecs,evals] = eig(covmat);

% Step 3: sort eigenvectors by descending eigenvalue magnitude
[evals,idx] = sort(diag(evals),'descend');
% each column in evecs corresponds to an eigenvector, so we need to
% rearrange them by columns
evecs = evecs(:,idx);

% Step 4: -

% Step 5: convert eigenvalues to percent variance 
evals = 100*evals./sum(evals);

% Step 6: -

% principal compo8nent time series (Step 4 but for time series)
pcts = evecs(:,1)'*EEG.data;

%% plotting time

figure(2), clf

subplot(2,3,1)
plot(evals(1:20),'ko-','markerfacecolor','w','linew',2)
title('Eigenspectrum'), axis square
ylabel('Percent variance'), xlabel('Component number')

% topographical map of first eigenvector
subplot(2,3,2)
topoplotIndie(evecs(:,1),EEG.chanlocs,'numcontour',0,'shading','interp');
title('PC topomap')

% topographical map of dipole (ground truth)
subplot(2,3,3)
topoplotIndie(-lf.Gain(:,1,diploc), EEG.chanlocs,'numcontour',0,'shading','interp');
title('Ground truth topomap')

% plot time series
subplot(2,1,2)
plot(EEG.times,smooth(pcts,50), ...
     EEG.times,smooth(EEG.data(31,:),10), ...
     EEG.times,-dipole_data(:,diploc)*100 ,'linew',2 )
legend({'PC';'Chan';'Dipole'})
xlabel('Time (s)')


