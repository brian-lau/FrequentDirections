% Reproduce Figure 5 Desai et al. using Caltech Birds dataset
%
%     Desai, Ghashami, & Phillips (2016). Improved practical matrix sketching 
%       with guarantees. IEEE Transactions on Knowledge & Data Engineering,
%       28(7), 1678-1690

clear
% Check that Birds data exists somewhere
if ~exist('image_attribute_labels.txt','file')
   help('BirdsReader');
   error('Birds data must be downloaded first');
end

% Reader for Birds data
BR = BirdsReader('filename','image_attribute_labels.txt');

% Load entire data set into memory
BR.blockSize = inf;
A = BR();

% Desai et al. demean for this dataset
A = bsxfun(@minus,A,mean(A));

k = 20:10:100;
alpha = [1 1 0 0.2 0.2];
fast = [false true false false true];
id = {'FD' 'FastFD' 'iSVD' '0.2FD' 'Fast 0.2FD'};
symbol = ['s' 'h' '*' 'p' 'd'];
color = ['b' 'g' 'r' 'c' 'm'];

tic;
[U,S,V] = svd(A);
bruteRuntime = toc;
m = 10;
Am = U(:,1:m)*S(1:m,1:m)*V(:,1:m)'; % For projection error

%% This can take a little time
count = 1;
for kk = k
   for m = 1:numel(alpha)
      [kk m]
      if fast(m)
         % In Desai et al, fast variants did not double the rows of B, and
         % returned B. Their fast variants can yield sketch that sometimes 
         % only using half of its rows (i.e., the lower half of B may contain 
         % all zeroes or actual data samples). To replicate this, we halve
         % k and keep the 'full' sketch.
         k2 = kk/2;
         fullsize = true;
      else
         k2 = kk;
         fullsize = false;
      end
      sketcher = FrequentDirections(k2,'fast',fast(m),'alpha',alpha(m));
      tic;
      sketcher(A);
      runtime(count,m) = toc;
      coverr(count,m) = sketcher.coverr(A,fullsize);
      projerr(count,m) = sketcher.projerr(A,Am,10,fullsize);
   end
   count = count + 1;
end

if 1 % Include randomized sparse variant
   rng(1234);
   n = size(A,1);
   nbetak = 50;
   reps = 10;
   
   for r = 1:reps
      count = 1;
      for kk = k
         k2 = kk/2;
         
         sketcher = FrequentDirections(k2,'sparse',true);
         sketcher.beta = n/nbetak/k2;
         tic;
         sketcher(A);
         temp_runtime(count,r) = toc;
         temp_coverr(count,r) = sketcher.coverr(A,true);
         temp_projerr(count,r) = sketcher.projerr(A,Am,10,true);
         nSVD(count,r) = sketcher.nSVD;
         nSparseEmbed(count,r) = sketcher.nSparseEmbed;
         count = count + 1;
      end
   end
   
   id{end+1} = 'SpFD50';
   symbol(end+1) = 'v';
   color(end+1) = 'k';
   
   runtime = [runtime , mean(temp_runtime,2)];
   coverr = [coverr , mean(temp_coverr,2)];
   projerr = [projerr , mean(temp_projerr,2)];
end

%% Plot
figure;
subplot(131);
hold on;
for i = 1:numel(id)
   g1(i) = plot(k',coverr(:,i),[symbol(i) '-'],'color',color(i),...
      'Markerfacecolor',color(i));
end
legend(id);
axis([min(k) max(k) 0 0.07]);
xlabel('Sketch size');
ylabel('Covariance error');

subplot(132);
hold on;
for i = 1:numel(id)
   g1(i) = plot(k',projerr(:,i),[symbol(i) '-'],'color',color(i),...
      'Markerfacecolor',color(i));
end
axis tight
axis([min(k) max(k) 0.975 1.175])
xlabel('Sketch size');
ylabel('Projection error');

subplot(133);
hold on;
for i = 1:numel(id)
   g1(i) = plot(k',runtime(:,i),[symbol(i) '-'],'color',color(i),...
      'Markerfacecolor',color(i));
end
plot([k(1) k(end)],[bruteRuntime bruteRuntime],'k--');
axis tight
yy = get(gca,'ylim');
axis([min(k) max(k) 0 yy(2)]);
xlabel('Sketch size');
ylabel('Runtime (seconds)');
