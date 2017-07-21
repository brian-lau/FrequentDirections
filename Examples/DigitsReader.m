%    DR = DigitsReader() returns a System object that streams images from
%    the MNIST database of handwritten digits.
%
%    The binary data files must be somewhere in the Matlab path.
%    Correctly formatted data is available from Yann Lecun's website:
%
%    http://yann.lecun.com/exdb/mnist/
%
%    where the image data and labels are provided separately. For example,
%    to stream the training data set (saved in 'train-images-idx3-ubyte'
%    and 'train-labels-idx1-ubyte'):
%
%    DR = DigitsReader('imageFilename','train-images-idx3-ubyte',...
%                      'labelFilename','train-labels-idx1-ubyte');
%    
%    while ~DR.isDone()
%       [image,label] = DR.step();
%       % do some processing on the current image
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

classdef DigitsReader < matlab.System & matlab.system.mixin.FiniteSource
   
   properties(Nontunable)
      imageFilename = 'train-images-idx3-ubyte'
      labelFilename = 'train-labels-idx1-ubyte'
   end
   
   properties
      blockSize = 1      % # of images to read on each iteration
   end
   
   properties(SetAccess = private,GetAccess = public)
      imageFID = -1      % file id for image data
      labelFID = -1      % file id for corresponding label data
      
      nRows              % # of rows for each image
      nCols              % # of columns for each image
      nImages            % # of images
      currentImageCount  % index of current image
   end
   
   methods
      function self = DigitsReader(varargin)
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
         
         self.nImages = fread(self.imageFID,1,'int32',0,'ieee-be');
         self.nRows = fread(self.imageFID,1,'int32',0,'ieee-be');
         self.nCols = fread(self.imageFID,1,'int32',0,'ieee-be');
         
         nLabels = fread(self.labelFID,1,'int32',0,'ieee-be');
         assert(self.nImages==nLabels,'mismatch');
         self.currentImageCount = 0;
      end
      
      function resetImpl(self)
         goToStartOfData(self);
         self.currentImageCount = 0;
      end
      
      function [images,labels] = stepImpl(self)
         n = self.blockSize;
         nRows = self.nRows;                                    %#ok<*PROP>
         nCols = self.nCols;
         spf = nRows*nCols;
         imageFID = self.imageFID;
         labelFID = self.labelFID;
         
         [images,count] = fread(imageFID,n*spf,'unsigned char');
         if count > 0
            images = reshape(images,nCols,nRows,count/spf);
            images = permute(images,[2 1 3]);
         end
         
         labels = fread(labelFID,n,'unsigned char');
         if ~isempty(labels)
            self.currentImageCount = self.currentImageCount + count/spf;
         end
      end
      
      function releaseImpl(self)
         fclose(self.imageFID);
         self.imageFID = -1;

         fclose(self.labelFID);
         self.labelFID = -1;
      end
      
      function tf = isDoneImpl(self)
         tf = logical(feof(self.imageFID));
      end
   end
   
   methods(Access = private)
      function getWorkingFID(self)
         if(self.imageFID < 0)
            [self.imageFID, err] = fopen(self.imageFilename,'rb');
            if ~isempty(err)
               error(message('DigitsReader:fileError',err));
            else
               magic = fread(self.imageFID,1,'int32',0,'ieee-be');
               assert(magic == 2051,['Bad magic # in ' self.imageFilename]);
            end
         end
         
         if(self.labelFID < 0)
            [self.labelFID, err] = fopen(self.labelFilename,'rb');
            if ~isempty(err)
               error(message('DigitsReader:fileError',err));
            else
               magic = fread(self.labelFID,1,'int32',0,'ieee-be');
               assert(magic == 2049,['Bad magic # in ' self.labelFilename]);
            end
         end
      end
      
      function goToStartOfData(self)
         % Skip header
         fseek(self.imageFID,4*4,'bof');
         fseek(self.labelFID,2*4,'bof');
      end
   end
end


