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
d = 28*28;
k = 64;
alpha = 0.2;
sketcher = FrequentDirections(d,k,'alpha',alpha);

% Process streamed data samples
tic;
while ~DR.isDone()
   [image,label] = DR.step();
   sketcher(image(:)');
end
toc

% Plot singular vectors
figure;
[~,~,V] = svd(get(sketcher),'econ');
for i = 1:36
   subplot(6,6,i);
   imshow(reshape(V(:,i),28,28),[]);
   title(strcat('SV ',int2str(i)));
end