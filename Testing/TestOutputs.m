classdef TestOutputs < matlab.unittest.TestCase
   methods (Test)
      function no_sketch(testCase)
         k = 16;
         fd = FrequentDirections(k);
         
         testCase.assertError(@() get(fd),...
            'FrequentDirections:NoOutput');
      end
            
      function bad_merge(testCase)
         fd1 = FrequentDirections(16);
         fd2 = FrequentDirections(8);
         
         testCase.assertError(@() merge(fd1,1),...
            'FrequentDirections:BadInput');

         testCase.assertError(@() merge(fd1,fd2),...
            'FrequentDirections:BadInput');

         fd2.alpha = .5;
         testCase.assertError(@() merge(fd1,fd2),...
            'FrequentDirections:BadInput');

         fd2.alpha = 1;
         fd2.fast = false;
         testCase.assertError(@() merge(fd1,fd2),...
            'FrequentDirections:BadInput');

         fd2.k = 16;
         fd1(rand(32,16));
         fd2(rand(32,32));

         testCase.assertError(@() merge(fd1,fd2),...
            'FrequentDirections:BadInput');
      end
   end
end

