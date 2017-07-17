% FREQUENTDIRECTIONS          Streaming low-rank matrix approximation
% 
%     sketcher = FrequentDirections(d,k,varargin)
%
%     Given an [n x d] matrix A, builds a [k x d] sketch B, where typically
%     k << n, using the Frequent Directions algorithm (Liberty 2013). This 
%     object works for matrices that are stored completely in-memory as well 
%     as data streams (see examples).
%
%     Implementing a number of FD variants (Desai et al 2016),
%        Classic FD: alpha = 1, fast = false
%           As defined in Liberty (2013)
%        Fast FD: alpha = 1, fast = true
%           Fast variant of FD that ensures that at most half of the rows
%           (k*alpha/2) of B are zeroed at each iteration (Liberty 2013). 
%           We follow the convention of Desai et al (2016), ensuring at most 
%           half of the rows of the sketch are all zeros after each rank 
%           reduction step.
%           Reduces the runtime from O(ndk^2) to O(ndk) at the expense of a 
%           sketch sometimes only using half of its rows (i.e., the lower 
%           half of B may contain all zeroes or actual data samples). To 
%           make identical to Liberty's version, double k.
%        Incremental SVD (iSVD): alpha = 0, fast = false
%        Parameterized FD: alpha = scalar in (0,1), fast = false
%        Fast Parameterized FD: alpha = scalar in (0,1), fast = true
%           alpha = 0.2, fast = true produces 'Fast 0.2FD' in Desai et al.
%
%     INPUTS
%     d - data dimensionality (rows are samples, columns are dimensions)
%     k - scalar in [1,d], sketch size. Note that this is frequently
%         referred to as l (ell) in references and other implementations
%
%     OPTIONAL (as name/value pairs, order irrelevant)
%     fast       - boolean, true indicates fast algorithm (default = TRUE)
%     alpha      - scalar in [0,1], controls fraction of sketch rows zeroed
%                  on each rank reduction (default = 1)
%     monitor    - boolean, true plots singular values at each rank reduction
%                  (default = FALSE)
%     figureAxis - axis handle for use when monitor = TRUE
%
%     PROPERTIES
%     d - data dimensionality determined by data given to object
%         (rows are samples, columns are dimensions)
%
%     METHODS
%     step    - given a [n x d] matrix, runs FD until all samples consumed
%     get     - returns sketch, B [k x d]
%     coverr  - given [n x d] matrix A, returns covariance error of sketch
%               ||A'A - B'B||_2 / ||A||_F^2
%     release - release resources to change parameters
%     reset   - delete current sketch and reset counters
%
%     EXAMPLE
%     k = 64;          % sketch size
%     monitor = false; % set true if you want to watch evolution
%
%     % Initialize object
%     sketcher = FrequentDirections(k,'monitor',monitor);
%     
%     d = 512;         % data dimensionality
%     count = 0;
%     while count < 1000
%        data = randn(1,d);   % Random sample
%        sketcher(data);      % consume sample
%        count = count + 1;
%     end
%
%     REFERENCE
%     Desai, Ghashami, & Phillips (2016). Improved practical matrix sketching 
%       with guarantees. IEEE Transactions on Knowledge & Data Engineering,
%       28(7), 1678-1690
%     Liberty (2013). Simple and deterministic matrix sketching. In 
%       Proceedings of the 19th ACM SIGKDD international conference on 
%       Knowledge discovery and data mining, 581-588

%     $ Copyright (C) 2017 Brian Lau, brian.lau@upmc.fr $
%     The full license and most recent version of the code can be found at:
%     https://github.com/brian-lau/
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.

% TODO
% o  Trigger final SVD for fast = true?

classdef FrequentDirections < matlab.System

   properties(Dependent)
      d
   end
   
   properties(Nontunable)
      k                 % sketch size
      alpha = 1         % [0,1] skrinkage control, 0 = iSVD, 1 = original FD
      fast = true       % boolean
   end
   
   properties
      monitor = false   % true plots singular values to axis
      figureAxis        % axis handle for plotting singular values
   end
   
   properties(Access = private)
      d_
      B_                % sketch
      reduceRankHandle_ % handle to rank reduction algorithm
      n_                % counter tracking # of data samples consumed
      count_            % counter tracking # of SVD calls
   end
   
   methods
      function self = FrequentDirections(varargin)
         setProperties(self,nargin,varargin{:},'k');
         %setProperties(self,nargin,varargin{:},'d','k');
      end
      
      function set.d_(self,d)
         if ~isempty(d)
            d = fix(d);
            assert(d>0,'d must be an integer > 0');
            assert(self.k<=d,'Sketch size k must <= data dimensionality');
            self.d_ = d;
         end
      end
      
      function set.k(self,k)
         k = fix(k);
         assert(k>0,'k must be an integer > 0');
         self.k = k;
      end
      
      function set.alpha(self,alpha)
         assert((alpha>=0)&&(alpha<=1),'alpha must be in [0,1]');
         self.alpha = alpha;
      end
      
      function set.fast(self,fast)
         if fast
            self.fast = true;
         else
            self.fast = false;
         end
      end
      
      function d = get.d(self)
         d = self.d_;
      end
      
      function B = get(self)
         B = self.B_;
      end
      
      function err = coverr(self,A)
         B = get(self);
         err = norm(A'*A - B'*B)/norm(A,'fro')^2;
      end
      
      function plot(self,S,Sprime)
         s = diag(S);
         sprime = diag(Sprime);
         ind = 1:numel(s);
         
         if isempty(self.figureAxis) || ~self.figureAxis.isvalid
            ax = subplot(1,1,1);
            self.figureAxis = ax;
         else
            ax = self.figureAxis;
         end
         
         if isempty(ax.Children)
            hold on;
            plot(ind,s,'ro');
            plot(ind,sprime,'bs');
            ax.XLim = [ind(1) ind(end)];
            ax.YLabel.String = 'Singular value';
         else
            ax.Children(1).YData = s;
            ax.Children(2).YData = sprime;
         end
         
         ax.Title.String = {sprintf('d=%g, k=%g, fast=%g, alpha=%1.2f',...
            self.d,self.k,self.fast,self.alpha) ...
            sprintf('#data=%g, #SVD=%g',self.n_,self.count_)};
         
         drawnow;
      end
   end
   
   methods(Access = protected)
      function setupImpl(self,A)
         if self.fast
            self.reduceRankHandle_ = @self.reduceRankFast;
         else
            self.reduceRankHandle_ = @self.reduceRankOriginal;
         end
         
         [~,p] = size(A);
         self.d_ = p;
         self.B_ = zeros(self.k,p);
      end
      
      function stepImpl(self,A)
         if isempty(A)
            return;
         end
         
         [n,p] = size(A);
         assert(p==self.d_,...
            'Input dimensionality does not match past samples!');

         k = self.k;                                          %#ok<*PROPLC>
         alpha = self.alpha;
         reduceRank = self.reduceRankHandle_;
         B = self.B_;

         %% Generic Frequent Directions algorithm
         count = 0; % keep track of SVD calls
         ind = [];  % keep track of zero rows of B
         i = 1;     % keep track of data samples appended
         while i <= n
            if isempty(ind)
               % Index all-zero rows of B
               ind = find(~any(B,2));
            end
            
            if ~isempty(ind)
               % Insert next data sample into first non-zero row of B
               B(ind(1),:) = A(i,:);
               ind(1) = [];
               i = i + 1;
            else
               % Update sketch
               [~,S,V] = svd(B,'econ');
               Sprime = reduceRank(S,k,alpha);
               B = Sprime*V';
               count = count + 1;
               
               if self.monitor
                  plot(self,S,Sprime);
               end
            end
         end
         
         self.B_ = B;
         self.count_ = self.count_ + count;
         self.n_ = self.n_ + i - 1;
      end
      
      function releaseImpl(self)
         self.B_ = [];
         self.d_ = [];
         if self.monitor
            close(self.figureAxis.Parent);
         end
      end
      
      function resetImpl(self)
         self.n_ = 0;
         self.count_ = 0;
      end
   end
   
   methods(Static)
      function Sprime = reduceRankOriginal(S,k,alpha)
         s = diag(S);
         sprime = zeros(size(s));
         
         skip = min(floor(k*(1-alpha)) + 1,k);
         if skip > 1
            sprime(1:skip) = s(1:skip);
         end
         
         if skip <= k
            dirac = s(k)^2;
            sprime(skip:end) = sqrt( s(skip:end).^2 - dirac );
         end
         
         Sprime = diag(sprime);
      end
      
      function Sprime = reduceRankFast(S,k,alpha)
         % Fast variant of FD that ensures that at most half of the rows
         % (k*alpha/2) of B are zeroed at each iteration
         s = diag(S);
         sprime = zeros(size(s));
         
         skip = floor(k*(1-alpha)/2) + 1;
         if skip > 1
            sprime(1:skip) = s(1:skip);
         end
         
         dirac_ind = k - floor(k*alpha/2) + 1;
         if (skip < k) && (dirac_ind <= k)
            dirac = s(dirac_ind)^2;
            sprime(skip:end) = sqrt( max(s(skip:end).^2 - dirac,0) );
         end
         
         Sprime = diag(sprime);
      end
   end
end
