d = 500;
k = 64;
monitor = true;
sketcher = FrequentDirections(k,'monitor',monitor);

%% Stream sample-by-sample
count = 0;
while count < 1000
   data = randn(1,d);
   sketcher(data);
   count = count + 1;
end

% Retrieve sketch
B = sketcher.get();

% Do something with sketch, e.g., approximate covariance matrix
covA = B'*B;

% To sketch a different matrix, release resources
sketcher.release();

%% Stream blocks of samples
blksz = 10;
count = 0;
while count < 1000
   data = randn(blksz,d);
   sketcher(data);
   count = count + blksz;
end
