% Reproduce part of Figure 4 of Teng & Chu using Caltech Birds dataset
%
%     Teng & Chu (2017). Low-Rank approximation via sparse frequent directions.
%       arXiv preprint arXiv:1705.07140.

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

n = size(A,1);

k = 50:10:150;
sp = [true true true false];
nbetak = [5 10 50 1];
id = {'SpFD5' 'SpFD10' 'SpFD50' 'FastFD'};
symbol = ['^' 'd' '+' 's'];
color = ['g' 'g' 'g' 'k'];

tic;
[U,S,V] = svd(A);
bruteRuntime = toc;
m = 50;
Am = U(:,1:m)*S(1:m,1:m)*V(:,1:m)'; % For projection error

%% This can take a little time
count = 1;
for kk = k
   for m = 1:numel(nbetak)
      [kk m]

      sketcher = FrequentDirections(kk,'alpha',1,'sparse',sp(m));
      sketcher.beta = n/nbetak(m)/kk;
      
      tic;
      sketcher(A);
      runtime(count,m) = toc;
      
      coverr(count,m) = sketcher.coverr(A);
      
      Am_ = sketcher.approx(A,50);
      projerr(count,m) = norm(A-Am_,'fro')/norm(A-Am,'fro');
      
      nSVD(count,m) = sketcher.nSVD;
      nSparseEmbed(count,m) = sketcher.nSparseEmbed;
   end
   count = count + 1;
end

%% Plot
figure;
subplot(131);
hold on;
for i = 1:numel(id)
   g1(i) = plot(k',coverr(:,i),'-','color',color(i));
end
legend(id);
for i = 1:numel(id)
   g2 = plot(k',coverr(:,i),symbol(i),...
      'Color',g1(i).Color,'Markerfacecolor',g1(i).Color);
end
axis([min(k) max(k) 0 0.07]);
xlabel('Sketch size');
ylabel('Covariance error');

subplot(132);
hold on;
for i = 1:numel(id)
   g1(i) = plot(k',projerr(:,i),'-','color',color(i));
end
%legend(id);
for i = 1:numel(id)
   g2 = plot(k',projerr(:,i),symbol(i),...
      'Color',g1(i).Color,'Markerfacecolor',g1(i).Color);
end
axis tight
axis([min(k) max(k) 1 1.4])
xlabel('Sketch size');
ylabel('Projection error');

subplot(133);
hold on;
for i = 1:numel(id)
   g1(i) = plot(k',runtime(:,i),'-','color',color(i));
end
%legend(id);
for i = 1:numel(id)
   g2 = plot(k',runtime(:,i),symbol(i),...
      'Color',g1(i).Color,'Markerfacecolor',g1(i).Color);
end
plot([k(1) k(end)],[bruteRuntime bruteRuntime],'k--');
axis tight
yy = get(gca,'ylim');
axis([min(k) max(k) 0 yy(2)]);
xlabel('Sketch size');
ylabel('Runtime (seconds)');
