% Figure 5 Desai et al.

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

k = 20:20:100;
alpha = [1 1 0 0.2 0.2];
fast = [false true false false true];
id = {'FD' 'FastFD' 'iSVD' '0.2FD' 'Fast 0.2FD'};
symbol = ['s' 'p' '*' 'h' 'd'];

count = 1;
for kk = k
   for m = 1:numel(alpha)
      [kk m]
      sketcher = FrequentDirections(kk,'fast',fast(m),'alpha',alpha(m));
      sketcher(A);
      coverr(count,m) = sketcher.coverr(A);
      %projerr(count,1) = sketcher.projerr(A);
   end
   count = count + 1;
end

subplot(211);
g1 = plot(k',coverr);
legend(id);
hold on;
for i = 1:numel(id)
   g2 = plot(k',coverr(:,i),symbol(i),...
      'Color',g1(i).Color,'Markerfacecolor',g1(i).Color);
end
axis([min(k) max(k) 0 0.07])

subplot(212);
g1 = plot(k',projerr);
legend(id);
hold on;
for i = 1:numel(id)
   g2 = plot(k',projerr(:,i),symbol(i),...
      'Color',g1(i).Color,'Markerfacecolor',g1(i).Color);
end
axis([min(k) max(k) 0.975 1.175])
