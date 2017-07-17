d = 500;
sketcher = FrequentDirections(d,64,'alpha',1);

% Stream sample-by-sample
count = 0;
while count < 1000
   data = randn(1,d);
   sketcher(data);
   count = count + 1;
end

sketcher.release();

% Stream blocks of samples
blksz = 10;
count = 0;
while count < 1000
   data = randn(blksz,d);
   sketcher(data);
   count = count + blksz;
end
