%%
%     COURSE: Dimension reduction and source separation in neuroscience
%    SECTION: Dimension reduction with PCA
%      VIDEO: MATLAB: state-space representation via PCA
% Instructor: mikexcohen.com
%
%%

%% load data

load sampleEEGdata.mat
EEG.data = double( EEG.data );

%% compute PCA

erp = mean(EEG.data,3);

data = bsxfun(@minus,erp,mean(erp,2));
covd = data*data'/(EEG.pnts-1);

[evecs,evals] = eig( covd );

% sort according to eigenvalues
[evals,sidx] = sort(diag(evals),'descend');
evecs = evecs(:,sidx);
evals = 100*evals/sum(evals);

% principal component time series
pc_timeseries = evecs(:,1:2)'*erp;

%% plot PC results in time-voltage space

figure(1), clf

% eigenspectrum
subplot(231)
plot(evals(1:20),'ko-','markerfacecolor','w','linew',2)
title('Eigenspectrum'), axis square
ylabel('Percent variance'), xlabel('Component number')

% topographical map of first eigenvector
subplot(232)
topoplotIndie(evecs(:,1),EEG.chanlocs,'numcontour',0,'shading','interp','electrodes','numbers');
title('PC1 topomap')

subplot(233)
topoplotIndie(evecs(:,2),EEG.chanlocs,'numcontour',0,'shading','interp','electrodes','numbers');
title('PC2 topomap')


% plot time series
subplot(212)
plot(EEG.times,pc_timeseries,'linew',2 )
legend({'PC1';'PC2'})
xlabel('Time (s)')

%% plot state-space representation

% time windows for plotting in ms
timewin = [0 250; 250 600; 600 800];


% smooth the PC time series
pcts = pc_timeseries;
k = 5; % smoothing kernel
% running mean smoothing filter 
for ti=k+1:EEG.pnts-k-1
    pcts(1,ti) = mean(pc_timeseries(1,ti-k:ti+k));
    pcts(2,ti) = mean(pc_timeseries(2,ti-k:ti+k));
end



% now plot the data
figure(2), clf, hold on
leg = cell(1,size(timewin,1));
for twini=1:size(timewin,1)
    
    % convert to indices
    tidx = dsearchn(EEG.times',timewin(twini,:)');
    
    % and plot
    plot(pcts(1,tidx(1):tidx(2)),pcts(2,tidx(1):tidx(2)),'linew',3)
    
    % define legend
    leg{twini} = [ num2str(timewin(twini,1)) '-' num2str(timewin(twini,2)) ' ms' ];
end

xlabel('PC 1'), ylabel('PC 2')
legend(leg)

%%
