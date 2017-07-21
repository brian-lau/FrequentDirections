%    BR = BirdsReader() returns a System object that streams data from
%    the Caltech-UCSD Birds-200-2011.
%
%    This object returns the binary attributes (d = 312) for each sample
%    (n = 11,788). The correctly formatted data file is available here:
%
%    http://www.vision.caltech.edu/visipedia/CUB-200-2011.html
%
%    where the attribute data is the text file 'image_attribute_labels.txt'. 
%    For example, to stream the entire data set:
%
%    BR = BirdsReader('filename','image_attribute_labels.txt');
%    
%    while ~BR.isDone()
%       attributes = BR.step();
%       % do some processing on the current attributes
%    end

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


classdef BirdsReader < matlab.System & matlab.system.mixin.FiniteSource
   
   properties(Nontunable)
      filename = 'image_attribute_labels.txt'
      certainty = false
   end
   
   properties(SetAccess = immutable)
      formatSpec = '%f %f %f %f %f %*[^\n]';
      nAttributes = 312;
   end
   
   properties
      blockSize = 1      % # of images to read on each iteration
   end
   
   properties(SetAccess = private,GetAccess = public)
      FID = -1           % file id for image data
      currentCount       % index of current image
   end
   
   methods
      function self = BirdsReader(varargin)
         setProperties(self,nargin,varargin{:});
      end
      
      function set.blockSize(self,blockSize)
         blockSize = fix(blockSize);
         assert(blockSize>0,'blockSize must be >= 1');
         self.blockSize = blockSize;
      end
      
      function goToImage(self,i)
         if self.imageFID == -1
            self.setup();
         end
         assert(i<self.nImages,'index must be <= total number of images');
         nRows = self.nRows;                                  %#ok<*PROPLC>
         nCols = self.nCols;
         spf = nRows*nCols;

         fseek(self.imageFID,4*4 + (i-1)*spf,'bof');
         fseek(self.labelFID,2*4 + (i-1),'bof');
      end
   end
   
   methods(Access = protected)
      function setupImpl(self)
         getWorkingFID(self);
      end
      
      function resetImpl(self)
         goToStartOfData(self);
         self.currentCount = 0;
      end
      
      function [attributes,id] = stepImpl(self)
         nAttributes = self.nAttributes;                        %#ok<*PROP>
         n = nAttributes*self.blockSize;
         if self.certainty
            col = 4;
         else
            col = 3;
         end

         data = textscan(self.FID,self.formatSpec,n,...
            'Delimiter','Whitespace',...
            'CollectOutput',true,...
            'MultipleDelimsAsOne',true);
         
         id = unique(data{1}(:,1));
         if self.blockSize == 1
            attributes = data{1}(:,col)';
         else
            x = data{1}(:,col);
            attributes = reshape(x,nAttributes,numel(x)/nAttributes)';
         end

         if ~isempty(id)
            self.currentCount = id(end);
         end
      end
      
      function releaseImpl(self)
         fclose(self.FID);
         self.FID = -1;
      end
      
      function tf = isDoneImpl(self)
         tf = logical(feof(self.FID));
      end
   end
   
   methods(Access = private)
      function getWorkingFID(self)
         if(self.FID < 0)
            [self.FID, err] = fopen(self.filename,'rt');
            if ~isempty(err)
               error(message('BirdsReader:fileError',err));
            end
         end
      end
      
      function goToStartOfData(self)
         fseek(self.FID,0,'bof');
      end
   end
end


