% Check that MNIST data exists somewhere
if ~(exist('train-images-idx3-ubyte','file') ...
      && exist('train-labels-idx1-ubyte','file'))
   help('DigitsReader');
   error('MNIST data must be downloaded first');
end

% Construct stream for MNIST digits data
DR = DigitsReader('imageFilename','train-images-idx3-ubyte',...
                  'labelFilename','train-labels-idx1-ubyte');

% Sketch object
k = 64;
alpha = 0.2;
sketcher = FrequentDirections(k,'alpha',alpha);

% Process streamed data samples
tic;
while ~DR.isDone()
   [image,label] = DR.step(); % load data sample
   sketcher(image(:)');       % consume
end
toc

% Plot singular vectors (row space of sketch)
[B,V] = get(sketcher);

figure;
for i = 1:36
   subplot(6,6,i);
   imshow(reshape(V(:,i),DR.nCols,DR.nRows),[]);
   title(strcat('SV ',int2str(i)));
end
