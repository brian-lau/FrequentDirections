% Check that Birds data exists somewhere
if ~exist('image_attribute_labels.txt','file')
   help('BirdsReader');
   error('Birds data must be downloaded first');
end

% Construct stream for Birds data
BR = BirdsReader('filename','image_attribute_labels.txt');

% Initialize sketch object
k = 20;
alpha = .2;
sketcher = FrequentDirections(k,'alpha',alpha,'fast',false);

% Process streamed data samples
tic;
while ~BR.isDone()
   attributes = BR.step(); % load data sample
   sketcher(attributes);   % consume
end
toc