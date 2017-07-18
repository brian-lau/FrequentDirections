% Check that Birds data exists somewhere
if ~exist('image_attribute_labels.txt','file')
   help('BirdsReader');
   error('Birds data must be downloaded first');
end

% Construct stream for Birds data
BR = BirdsReader('filename','image_attribute_labels.txt');

% Sketch object
k = 20;
alpha = .8;
sketcher = FrequentDirections(k,'alpha',alpha,'fast',false);

% Process streamed data samples
tic;
while ~BR.isDone()
   attributes = BR.step();
   sketcher(attributes);
end
toc

%
BR.release();
BR.blockSize = 12000;
A = BR();

sketcher.coverr(A)