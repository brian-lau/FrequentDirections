% Sketches created using Frequent Directions are mergeable, meaning that 
% sketches of partitions of a data stream can be merged to create a single
% sketch that inherits the error bounds (Ghashami et al, 2016).

clear

k = 100;            % sketch size
d = 300;            % data dimensionality

% Create a dataset to be partitioned
rng(1);
blk = 5000;
nSketch = 8;
data = randn(nSketch*blk,d);

% Apply FrequentDirections separately to each partition
for i = 1:nSketch
   sketch{i} = FrequentDirections(k);
   ind = ((i-1)*blk+1):i*blk;
   sketch{i}(data(ind,:));
end

% Passing individual sketches to merge yields a new sketch
s = merge(sketch{:});

s.coverr(data)

% Compare the covariance error to that of a single (unmerged) sketch
s1 = FrequentDirections(k);
s1(data);
s1.coverr(data)