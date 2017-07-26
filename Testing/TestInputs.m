classdef TestInputs < matlab.unittest.TestCase
   methods (Test)
      function good_d(testCase)
         k = 16;
         fd = FrequentDirections(k);
         B = reshape(1:k^2,k,k)';
         fd(B);
         
         testCase.assertEqual(get(fd),B);
         testCase.assertEqual(get(fd,true),[B ; zeros(size(B))]);
      end
      
      function bad_d(testCase)
         fd = FrequentDirections(16);
         
         testCase.assertError(@() fd([1 2 3]),...
            'FrequentDirections:BadDimension');
      end
      
      function bad_setup(testCase)
         fd = FrequentDirections(3);
         
         testCase.assertError(@() fd(randn(4,4,4)),...
            'FrequentDirections:BadDimension');
      end
      
      function bad_stream(testCase)
         fd = FrequentDirections(3);
         fd([1 1 1]);
         
         testCase.assertError(@() fd([1 1 1 1]),...
            'FrequentDirections:BadDimension');
      end
   end
end

