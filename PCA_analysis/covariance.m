% simulation parameters
N = 1000; % time points 
ch = 20; % channels
trials = 50;

% time vector normalised in radian units 6*pi -> 3 cycles
t = linspace(0,6*pi,N);

data = cell(1,2); % array with dimensions 1x2

% impose relationship accross channels (imposing correlation/covariance)
chan_rel = sin(linspace(0,2*pi,ch))';

for trial = 1:trials
  % repmat (A, M, N) : form a block matrix M x N, with a copy of matrix A
  % as each element
  % bsxfun(@times,A,B) : multiplication of matrices A and B
  
  % simulation 1, phase-locked 
  % channels x time points x trials
  data{1}(:,:,trial) = bsxfun(@times,repmat(sin(t),ch,1),chan_rel) + randn(ch,N);
  
  % simulation 2, non phase-locked
  data{2}(:,:,trial) = bsxfun(@times,repmat(sin(t+rand*2*pi),ch,1),chan_rel) + randn(ch,N);
end

figure(1), clf

for i = 1:2
  subplot(1,2,i)
  
  % plot a random tial
  trial2plot = ceil(rand*trials);
  
  % squeeze(X): Remove singleton dimensions from X and return the result.
  plot(t,bsxfun(@plus,squeeze(data{i}(:,:,trial2plot)),(1:ch)'*3))
  
  set(gca,'ytick',[]), axis tight
  xlabel('Time'), ylabel('Channels')
  title(['Dataset ' num2str(i) ', trial' num2str(trial2plot)])
end

%% compute covariance matrix element-wise
covmat1 = zeros(ch);

% compute the trial average (ERP)
data_avg = mean(data{1}, 3); % channels x time points

for chani = 1:ch
  for chanj = 1:ch
    % mean-center data
    subi = data_avg(chani,:) - mean(data_avg(chani,:));
    subj = data_avg(chanj,:) - mean(data_avg(chanj,:));
    covmat1(chani, chanj) = sum(subi.*subj)/(N-1);
  end
end

%% compute covariance by matrix multiplication

% mean-centered over time
data_avg_mean = bsxfun(@minus,data_avg,mean(data_avg,2));
covmat2 = data_avg2*data_avg2'/(N-1);

% check size to see if we have the right coveriance matrix
size(covmat2)

%% compute covariance with Matlab function
covmat3 = cov(data_avg');

size(covmat3)

figure(2), clf

subplot(1,3,1), imagesc(covmat1)
axis square, title('Element-wise')
xlabel('Channels'), ylabel('Channels')

subplot(1,3,2), imagesc(covmat2)
axis square, title('Matrix multiplication')
xlabel('Channels'), ylabel('Channels')

subplot(1,3,3), imagesc(covmat3)
axis square, title('MATLAB cov function')
xlabel('Channels'), ylabel('Channels')
