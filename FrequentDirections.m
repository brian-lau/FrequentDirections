% FREQUENTDIRECTIONS          Streaming low-rank matrix approximation
% 
%     sketcher = FrequentDirections(k,varargin)
%
%     Given an [n x d] matrix A, builds a [k x d] sketch B, where typically
%     k << n, using the Frequent Directions algorithm (Liberty 2013). This 
%     object works for matrices that are stored completely in-memory as well 
%     as data streams (see examples).
%
%     Data dimensionality is determined at run-time from data provided 
%     (rows are samples, columns are dimensions).
%
%     Implements a number of FD variants (Desai et al 2016),
%        Classic FD: alpha = 1, fast = false
%           As defined in Liberty (2013)
%        Fast FD: alpha = 1, fast = true
%           Fast variant of FD that ensures that at most half of the rows
%           of B are zeroed at each iteration (Liberty 2013). Reduces the 
%           runtime from O(ndk^2) to O(ndk), at the expense of double
%           storage size of sketch while algorithm is running.
%        Iterative SVD (iSVD): alpha = 0, fast = false
%        Parameterized FD: alpha = scalar in (0,1), fast = false
%        Fast Parameterized FD: alpha = scalar in (0,1), fast = true
%           alpha = 0.2, fast = true produces 'Fast 0.2FD' in Desai et al.
%        
%     INPUTS
%     k - scalar in [1,d], sketch size. Note that this is commonly referred 
%         to as l (ell) in references and other implementations
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
%     d - data dimensionality determined at run-time from data provided
%         (rows are samples, columns are dimensions)
%
%     METHODS
%     step    - Given a [n x d] matrix, runs FD until all samples consumed.
%               After the first call, object parameters are locked, and
%               subsequent steps must have the same number of columns (d),
%               and each step is used to build on the current sketch.
%               obj.step(A) is equivalent to obj(A)
%     get     - returns sketch, B [k x d]
%               Setting the input true (i.e. obj.get(true) as opposed to
%               obj.get() or get(obj)) will return a [2k x d] matrix when
%               fast = true.
%     coverr  - given [n x d] matrix A, returns covariance error of sketch
%               ||A'A - B'B||_2 / ||A||_F^2
%     projerr - given [n x d] matrix A, returns projection error of sketch
%               ||A - proj(A,B)||_F^2 / ||A - A_m||_F^2
%     release - delete current sketch & release resources to change parameters
%     reset   - reset counters
%
%     EXAMPLE
%     k = 16;            % sketch size
%     monitor = false;   % set true to watch evolution of singular values
%
%     % Initialize object
%     sketcher = FrequentDirections(k,'monitor',monitor);
%     
%     d = 64;            % data dimensionality
%
%     % Sketch matrix entirely in-memory
%     data = randn(1000,d);
%     sketcher(data);
%     get(sketcher)
%
%     % Sketch streaming data
%     release(sketcher); % release object to build new sketch
%     
%     d = 512;           % different data dimensionality
%     sketcher.k = 32    % change sketch size
%
%     count = 0;
%     while count < 1000
%        data = randn(1,d);   % random sample
%        sketcher(data);      % consume sample
%        count = count + 1;
%     end
%
%     REFERENCE
%     Desai, Ghashami, & Phillips (2016). Improved practical matrix sketching 
%       with guarantees. IEEE Transactions on Knowledge & Data Engineering,
%       28(7), 1678-1690
%     Ghashami et al (2016). Frequent directions: Simple and deterministic 
%       matrix sketching. SIAM Journal on Computing, 45(5), 1762-1792.
%     Liberty (2013). Simple and deterministic matrix sketching. In 
%       Proceedings of the 19th ACM SIGKDD international conference on 
%       Knowledge discovery and data mining, 581-588

%     $ Copyright (C) 2017 Brian Lau, brian.lau@upmc.fr $
%     The full license and most recent version of the code can be found at:
%     https://github.com/brian-lau/FrequentDirections
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
      d                 % data dimensionality (# columns of input matrix)
   end
   
   properties(Nontunable)
      k                 % sketch size
      alpha = 1         % [0,1] skrinkage control parameter, 0 = iSVD, 1 = original FD
      fast = true       % boolean
   end
   
   properties
      monitor = false   % true plots singular values to axis
      figureAxis        % axis handle for plotting singular values
   end
   
   properties(Access = private)
      d_                % data dimensionality
      k2_               % Temporary sketch size (doubled for fast=true)
      B_                % sketch
      reduceRankHandle_ % handle to rank reduction algorithm
      n_                % counter tracking # of data samples consumed
      count_            % counter tracking # of SVD calls
   end
   
   properties(SetAccess = immutable)
      version = '0.3.0' % Version string
   end
   
   methods
      function self = FrequentDirections(varargin)
         setProperties(self,nargin,varargin{:},'k');
      end
      
      function set.d_(self,d)
         if ~isempty(d)
            d = fix(d);
            assert(d>0,'FrequentDirections:BadDimension',...
               'd must be an integer > 0');
            assert(self.k<=d,'FrequentDirections:BadDimension',...
               'Sketch size k must be <= data dimensionality');
            self.d_ = d;
         end
      end
      
      function set.k(self,k)
         k = fix(k);
         assert(isscalar(k)&&(k>0),'FrequentDirections:BadInput',...
            'k must be an scalar integer > 0');
         self.k = k;
      end
      
      function set.alpha(self,alpha)
         assert(isscalar(alpha)&&(alpha>=0)&&(alpha<=1),...
            'FrequentDirections:BadInput',...
            'alpha must be in scalar in [0,1]');
         self.alpha = alpha;
      end
      
      function set.fast(self,fast)
         assert(isscalar(fast),'FrequentDirections:BadInput',...
            'fast must be a scalar boolean');
         self.fast = logical(fast);
      end
      
      function set.monitor(self,monitor)
         assert(isscalar(monitor),'FrequentDirections:BadInput',...
            'monitor must be a scalar boolean');
         self.monitor = logical(monitor);
      end
      
      function set.figureAxis(self,h)
         assert(isa(h,'matlab.graphics.axis.Axes'),...
            'FrequentDirections:BadInput',...
            'Input must be of type matlab.graphics.axis.Axes');
         
         self.figureAxis = h;
      end
      
      function d = get.d(self)
         d = self.d_;
      end
      
      % GET         Return matrix sketch
      %
      % INPUT
      % fullsize - boolean, only relevant when fast = true
      %            true indicates returning sketch with 2*k rows, the bottom
      %            half of which may be all zeroes or actual data samples
      %            default = FALSE
      %
      % OUTPUT
      % B        - [k x d] sketch
      %            [2k x d] sketch if fullsize = true && fast = true
      function B = get(self,fullsize)
         if nargin < 2
            fullsize = false;
         end
         
         assert(~isempty(self.B_),'FrequentDirections:NoOutput',...
            'No sketch to get yet!');
         
         if self.fast && ~fullsize
            B = self.B_(1:self.k,:);
         else
            B = self.B_;
         end
      end
      
      % COVERR      Covariance error
      function err = coverr(self,A,fullsize)
         if nargin < 3
            fullsize = false;
         end
         B = get(self,fullsize);
         err = norm(A'*A - B'*B)/norm(A,'fro')^2;
      end
      
      % PROJERR     Projection error
      function err = projerr(self,A,Am,m,fullsize)
         if nargin < 5
            fullsize = false;
         end
         
         if nargin < 4
            m = 10;
         end
         
         assert(m<=self.k,'m must be less than k');
         
         if (nargin < 3) || isempty(Am)
            % Rank m approximation of A
            [U,S,V] = svd(A);
            Am = U(:,1:m)*S(1:m,1:m)*V(:,1:m)';
         end
         
         B = get(self,fullsize);
         Bm = B(1:m,:);
         
         Am_ = A*Bm'*pinv(Bm*Bm')*Bm;
         err = norm(A-Am_,'fro')^2 / norm(A-Am,'fro')^2;
      end
      
      % MERGE       Merge separate sketches
      % 
      % Sketches created using Frequent Directions are mergeable, meaning 
      % that sketches of data stream partitions can be merged to create a 
      % single sketch that inherits the error bounds (Ghashami et al, 2016).
      %
      % INPUTS
      % Individual FrequentDirections objects to be merged
      %
      % OUTPUT
      % obj - a new FrequentDirections object containing merged sketch
      %
      % EXAMPLE
      % s1 = FrequentDirections(16);
      % s1(randn(1000,16));
      % s2 = FrequentDirections(16);
      % s2(randn(1000,16));
      % s = merge(s1,s2);
      %
      % SEE ALSO
      % exampleMerge
      function obj = merge(varargin)
         tf = all(cellfun(@(x) isa(x,'FrequentDirections'),varargin));
         assert(tf,'FrequentDirections:BadInput',...
            'Inputs must all be FrequentDirections objects.');
         
         k = cellfun(@(x) x.k,varargin);                        %#ok<*PROP>
         assert(numel(unique(k))==1,'FrequentDirections:BadInput',...
            'Merging sketches requires the same k');
         
         d = cellfun(@(x) x.d,varargin);
         assert(numel(unique(d))==1,'FrequentDirections:BadInput',...
            'Merging sketches requires the same d');
         
         alpha = cellfun(@(x) x.alpha,varargin);
         assert(numel(unique(alpha))==1,'FrequentDirections:BadInput',...
            'Merging sketches requires the same alpha');

         fast = cellfun(@(x) x.fast,varargin);
         assert(numel(unique(fast))==1,'FrequentDirections:BadInput',...
            'Merging sketches requires the same fast setting');
         
         obj = FrequentDirections(k(1),'alpha',alpha(1),'fast',fast(1));
         
         B = cellfun(@(x) get(x),varargin,'uni',false);
         B = cat(1,B{:});
         
         % FD on concatenated sketches
         obj.step(B);
         
         % Force a rank reduction if none performed
         if obj.count_ == 0
            obj.step(zeros(1,obj.d));
         end
         
         % Update count & n
         obj.n_ = sum(cellfun(@(x) x.n_,varargin));
         obj.count_ = obj.count_ + sum(cellfun(@(x) x.count_,varargin));
      end

   end
   
   methods(Access = protected)
      function setupImpl(self,A)
         assert(ismatrix(A),'FrequentDirections:BadDimension',...
            'Input must be 2D matrix.');
         [~,d] = size(A);
         self.d_ = d;

         if self.fast
            self.reduceRankHandle_ = @self.reduceRankFast;
            self.k2_ = self.k*2;
         else
            self.reduceRankHandle_ = @self.reduceRankOriginal;
            self.k2_ = self.k;
         end
         
         self.B_ = zeros(self.k2_,d);
      end
      
      function stepImpl(self,A)
         if isempty(A)
            return;
         end
         
         [n,p] = size(A);
         assert(p==self.d_,'FrequentDirections:BadDimension',...
            'Input dimensionality does not match past samples!');

         k = self.k2_;                                        %#ok<*PROPLC>
         alpha = self.alpha;
         reduceRank = self.reduceRankHandle_;
         B = self.B_;
         monitor = self.monitor;

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
               
               if monitor
                  plot(self,S,Sprime,i-1,count);
               end
            end
         end
         
         self.B_ = B;
         self.count_ = self.count_ + count;
         self.n_ = self.n_ + i - 1;
         
         if monitor
            plot(self,S,Sprime,self.n_,self.count_);
         end
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
      
      % PLOT        Plot singular values
      function plot(self,S,Sprime,n,count)
         s = diag(S);
         sprime = diag(Sprime);
         ind = 1:numel(s);
         
         if isempty(self.figureAxis) || ~self.figureAxis.isvalid
            figure;
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
            sprintf('#data=%g, #SVD=%g',n,count)};
         drawnow;
      end
   end
   
   methods(Static)
      % REDUCERANKORIGINAL    Original rank reduction
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
      
      % REDUCERANKFAST        Fast rank reduction
      %
      %     Fast variant of FD that ensures that at most half of the rows
      %     (k*alpha/2) of B are zeroed at each iteration.
      function Sprime = reduceRankFast(S,k,alpha)
         s = diag(S);
         sprime = zeros(size(s));
         
         skip = floor(k*(1-alpha)) + 1;
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
