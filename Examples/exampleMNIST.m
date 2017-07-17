DR = DigitsReader('imageFilename','train-images-idx3-ubyte',...
                  'labelFilename','train-labels-idx1-ubyte');

sketcher = FrequentDirections(28*28,64,'alpha',0.2);

tic;
while ~DR.isDone()
   [image,label] = DR.step();
   sketcher(image(:)');
end
toc

figure;
[U,~,~] = svd(get(sketcher)','econ');
for i = 1:36
   subplot(6,6,i);
   imshow(reshape(U(:,i)',28,28),[]);
   title(strcat('direction index = ',int2str(i)));
end
