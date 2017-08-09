% Reproduce part of Figure 4 & 5 of Teng & Chu 
%
%     Teng & Chu (2017). Low-Rank approximation via sparse frequent directions.
%       arXiv preprint arXiv:1705.07140.
%
% TODO
%   o average over reps since sparse method is not deterministic
%   o runtimes don't seem to match Teng & Chu (their vanilla FD is slow?)
%
clear

if 0 % Birds data
   DR = BirdsReader();
   p = 50; % Approximating rank (k in table 4.1)
   k = 50:10:150;
   % Load entire data set into memory
   DR.blockSize = inf;
   A = DR();
else % MNIST data
   DR = DigitsReader();
   p = 100; % Approximating rank (k in table 4.1)
   k = 100:10:200;
   % Load entire data set into memory
   DR.blockSize = inf;
   A = DR();
   
   A = reshape(A,28*28,60000)';
end

n = size(A,1);

sp = [true true true false];
nbetak = [5 10 50 1];
id = {'SpFD5' 'SpFD10' 'SpFD50' 'FastFD'};
symbol = ['^' 'd' '+' 's'];
color = ['g' 'g' 'g' 'k'];

tic;
[U,S,V] = svd(A,'econ');
bruteRuntime = toc;
Am = U(:,1:p)*S(1:m,1:p)*V(:,1:p)'; % For projection error

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
      
      Am_ = sketcher.approx(A,p);
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
