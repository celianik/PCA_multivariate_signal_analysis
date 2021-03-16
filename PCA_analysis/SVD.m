%% Simulate data

N = 10000; % time points
M = 20; % channels

% create time vector (radian units) -> 3 cycles
t = linspace(0,6*pi,N);

% set relationship between channels
chanrel = sin(linspace(0,2*pi,M))';

% generate data
data = bsxfun(@times,repmat(sin(t),M,1),chanrel) + randn(M,N);

%% PCA via eigendecomposition of cov matri
dataM = bsxfun(@minus,data,mean(data,2));
covMat = dataM*dataM'/(N-1);

% eigendecomposition of covariance
[evecs,evals] = eig(covMat);

% sort eigenvalues;
[evals,idx] = sort(diag(evals),'descend');
evecs = evecs(:,idx);

% convert eigenvalues to %
evals = 100*evals/sum(evals);

% time series for top principal component
eig_ts = evecs(:,1)'*data;

%% PCA via SVD
% [U,S,V] = svd(X) produces a diagonal matrix S, of the same 
% dimension as X and with nonnegative diagonal elements in
% decreasing order, and unitary matrices U and V so that
% X = U*S*V'
% [U,S,V] = svd(X,'econ') also produces the "economy size"
% decomposition
[U,S,V] = svd(data,'econ');

% convert singular values to % change
S = S.^2; % to make them comparable to the eigenvalues
s = 100*diag(S)/sum(diag(S));

%% Plotting time

figure(1), clf

% plot eigenvalue/singular-value spectrum
subplot(2,2,1), hold on
plot(evals,'bs-','markerfacecolor','w','markersize',10,'linew',2)
plot(s,'ro-','markerfacecolor','w','markersize',5)
xlabel('Component number'), ylabel('\lambda or \sigma')
legend({'eig';'svd'})
title('Eigenspectrum')


% show eigenvector/singular value
subplot(2,2,2), hold on
plot(evecs(:,1),'bs-','markerfacecolor','w','markersize',10,'linew',2)
plot(U(:,1),'ro-','markerfacecolor','w','markersize',5)
xlabel('Vector component'), ylabel('Weight')
title('Component weights')


% time series
subplot(2,1,2)
timevec = (0:N-1)/N;
plot(timevec,eig_ts, timevec,-V(:,1)*sqrt(S(1)),'.')
xlabel('Time (norm.)')
title('Component time series')
legend({'eig';'svd'})
zoom on
