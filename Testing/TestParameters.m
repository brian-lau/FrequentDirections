classdef TestParameters < matlab.unittest.TestCase
   methods (Test)
      function testNoArgs(testCase)
         fd = FrequentDirections();
         testCase.assertClass(fd,'FrequentDirections');
      end
      
      function good_k(testCase)
         k = 16;
         fd = FrequentDirections(k);
         testCase.assertEqual(fd.k,k);
         testCase.assertEmpty(fd.d);
         
         testCase.assertEqual(fd.alpha,1);
         testCase.assertTrue(fd.fast);
         testCase.assertFalse(fd.monitor);
         testCase.assertEmpty(fd.figureAxis);
      end
      
      function good_k_overload(testCase)
         k = 16;
         fd = FrequentDirections('k',k);
         testCase.assertEqual(fd.k,k);
         testCase.assertEmpty(fd.d);
         
         testCase.assertEqual(fd.alpha,1);
         testCase.assertTrue(fd.fast);
         testCase.assertFalse(fd.monitor);
         testCase.assertEmpty(fd.figureAxis);
      end
      
      function bad_k(testCase)
         testCase.assertError(@() FrequentDirections(-1),...
            'FrequentDirections:BadInput');

         testCase.assertError(@() FrequentDirections(0),...
            'FrequentDirections:BadInput');

         testCase.assertError(@() FrequentDirections([1 2]),...
            'FrequentDirections:BadInput');
      end
      
      function good_alpha(testCase)
         k = 16;
         alpha = 0.5;
         fd = FrequentDirections(k,'alpha',alpha);

         testCase.assertEqual(fd.alpha,alpha);
      end
      
      function bad_alpha(testCase)
         testCase.assertError(@() FrequentDirections(16,'alpha',-1),...
            'FrequentDirections:BadInput');

         testCase.assertError(@() FrequentDirections(16,'alpha',1.1),...
            'FrequentDirections:BadInput');

         testCase.assertError(@() FrequentDirections(16,'alpha',[1 2]),...
            'FrequentDirections:BadInput');

      end
      
      function good_fast(testCase)
         k = 16;
         fd = FrequentDirections(k,'fast',true);

         testCase.assertEqual(fd.fast,true);
         
         fd = FrequentDirections(k,'fast',false);
         
         testCase.assertEqual(fd.fast,false);
      end
      
      function bad_fast(testCase)
         testCase.assertError(@() FrequentDirections(16,'fast',[1 2]),...
            'FrequentDirections:BadInput');
      end

      function good_sparse(testCase)
         k = 16;
         fd = FrequentDirections(k,'sparse',true);

         testCase.assertEqual(fd.sparse,true);
         
         fd = FrequentDirections(k,'sparse',false);
         
         testCase.assertEqual(fd.sparse,false);
      end
      
      function bad_sparse(testCase)
         testCase.assertError(@() FrequentDirections(16,'fast',0,'sparse',1),...
            'FrequentDirections:BadInput');
      end
      
      function good_monitor(testCase)
         k = 16;
         fd = FrequentDirections(k,'monitor',true);

         testCase.assertEqual(fd.monitor,true);
         
         fd = FrequentDirections(k,'monitor',false);
         
         testCase.assertEqual(fd.monitor,false);
      end
      
      function bad_monitor(testCase)
         testCase.assertError(@() FrequentDirections(16,'monitor',[1 2]),...
            'FrequentDirections:BadInput');
      end
      
      function good_figureAxis(testCase)
         k = 16;
         figure;
         ax = subplot(1,1,1);
         fd = FrequentDirections(k,'figureAxis',ax);
         close;
         testCase.assertInstanceOf(fd.figureAxis,'matlab.graphics.axis.Axes');
      end
      
      function bad_figureAxis(testCase)
         testCase.assertError(@() FrequentDirections(16,'figureAxis',1),...
            'FrequentDirections:BadInput');
      end
   end
end

