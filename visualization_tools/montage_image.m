% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [bigImage, h] = montage_image(varargin)
%MONTAGE_IMAGE Returns the montage image. Original montage function is
%slightly modified to return the image itself. The image is also displayed
%when the function is called with no output variable or with exactly two
%output variables. 
%
%  help montage for details.

warning('off', 'images:removing:function');

[I,cmap,mSize,indices,displayRange] = parse_inputs_MONTAGEIMAGE(varargin{:});

if isempty(indices) || isempty(I)
   if (nargout == 0 || nargout == 2)
      h = imshow([]);
   end
   return;
end

% Function Scope
nFrames = numel(indices);
nRows = size(I,1);
nCols = size(I,2);

montageSize = calculateMontageSize(mSize);

bigImage = createMontageImage;

if isempty(cmap)
    if isempty(displayRange)
        num = numel(I(:,:,:,indices));
        displayRange(1) = min(reshape(I(:,:,:,indices),[1 num]));
        displayRange(2) = max(reshape(I(:,:,:,indices),[1 num]));
    end
    if (nargout == 0 || nargout == 2)
       h = imshow(bigImage, displayRange);
    end
else
   if (nargout == 0 || nargout == 2)
      h = imshow(bigImage,cmap);
   end
end

% if nargout > 0
%     h = hh;
% end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function montageSize = calculateMontageSize(mSize)

        if isempty(mSize) || all(isnan(mSize))
            %Calculate montageSize for the user
            
            % Estimate nMontageColumns and nMontageRows given the desired
            % ratio of Columns to Rows to be one (square montage).
            aspectRatio = 1;
            montageCols = sqrt(aspectRatio * nRows * nFrames / nCols);

            % Make sure montage rows and columns are integers. The order in
            % the adjustment matters because the montage image is created
            % horizontally across columns.
            montageCols = ceil(montageCols);
            montageRows = ceil(nFrames / montageCols);
            montageSize = [montageRows montageCols];
        
        elseif any(isnan(mSize))
            montageSize = mSize;
            nanIdx = isnan(mSize);
            montageSize(nanIdx) = ceil(nFrames / mSize(~nanIdx));

        elseif prod(mSize) < nFrames
            eid = sprintf('Images:%s:sizeTooSmall', mfilename);
            error(eid, ...
                  'SIZE must be big enough to include all frames in I.');

        else
            montageSize = mSize;
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function b = createMontageImage

        nMontageRows = montageSize(1);
        nMontageCols = montageSize(2);
        nBands = size(I, 3);

        sizeOfBigImage = [nMontageRows*nRows nMontageCols*nCols nBands 1];
        if islogical(I)
            b = false(sizeOfBigImage);
        else
            b = zeros(sizeOfBigImage, class(I));
        end
        
        rows = 1 : nRows;
        cols = 1 : nCols;
        k = 1;

        for i = 0 : nMontageRows-1
            for j = 0 : nMontageCols-1,
                if k <= nFrames
                    b(rows + i * nRows, cols + j * nCols, :) = ...
                        I(:,:,:,indices(k));
                else
                    return;
                end
                k = k + 1;
            end
        end

    end

end %MONTAGE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,cmap,montageSize,idxs,displayRange] = parse_inputs_MONTAGEIMAGE(varargin)

iptchecknargin(1, 8, nargin, mfilename);

% initialize variables
cmap = [];
montageSize = [];

charStart = find(cellfun('isclass', varargin, 'char'));

if iscell(varargin{1})
    %MONTAGE(FILENAMES.,..)
    [I,cmap] = getImagesFromFiles(varargin{1});
    
else
    %MONTAGE(I,...) or MONTAGE(X,MAP,...)
    I = varargin{1};
    iptcheckinput(varargin{1}, ...
        {'uint8' 'double' 'uint16' 'logical' 'single' 'int16'}, {}, ...
        mfilename, 'I, BW, or RGB', 1);
end

nframes = size(I,4);
displayRange = getrangefromclass(I);
idxs = 1:nframes;

if isempty(charStart)
    %MONTAGE(FILENAMES), MONTAGE(I) or MONTAGE(X,MAP)
    if nargin == 2
        %MONTAGE(X,MAP)
        cmap = validateColormapSyntax(I,varargin{2});
    end
    return;
end

charStart = charStart(1);

if charStart == 3
    %MONTAGE(X,MAP,Param1,Value1,...)
    cmap = validateColormapSyntax(I,varargin{2});
end

paramStrings = {'Size', 'Indices', 'DisplayRange'};
    
for k = charStart:2:nargin

    param = lower(varargin{k});
    inputStr = iptcheckstrs(param, paramStrings, mfilename, 'PARAM', k);
    valueIdx = k + 1;
    if valueIdx > nargin
        eid = sprintf('Images:%s:missingParameterValue', mfilename);
        error(eid, ...
            'Parameter ''%s'' must be followed by a value.', ...
            inputStr);
    end

    switch (inputStr)
        case 'Size'
            montageSize = varargin{valueIdx};
            iptcheckinput(montageSize,{'numeric'},...
                {'vector','positive'}, ...
                mfilename, 'Size', valueIdx);
            if numel(montageSize) ~= 2
                eid = sprintf('Images:%s:invalidSize',mfilename);
                error(eid, 'Size must be a 2-element vector.');
            end

            montageSize = double(montageSize);

        case 'Indices'
            idxs = varargin{valueIdx};
            iptcheckinput(idxs, {'numeric'},...
                {'vector','integer','nonnan'}, ...
                mfilename, 'Indices', valueIdx);

            invalidIdxs = ~isempty(idxs) && ...
                any(idxs < 1) || ...
                any(idxs > nframes);

            if invalidIdxs
                eid = sprintf('Images:%s:invalidIndices',mfilename);
                error(eid, ...
                    'An index in INDICES cannot be less than 1 %s', ...
                    'or greater than the number of frames in I.');
            end

            idxs = double(idxs);

        case 'DisplayRange'
            displayRange = varargin{valueIdx};
            % displayRange = checkDisplayRange(displayRange, mfilename);

    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cmap = validateColormapSyntax(I,cmap)

if isa(I,'int16')
    eid = sprintf('Images:%s:invalidIndexedImage',mfilename);
    error(eid, 'An indexed image can be uint8, uint16, %s', ...
        'double, single, or logical.');
end

iptcheckinput(cmap,{'double'},{},mfilename,'MAP',1);

if size(cmap,1) == 1 && prod(cmap) == numel(I)
    % MONTAGE(D,[M N P]) OBSOLETE
    eid = sprintf('Images:%s:obsoleteSyntax',mfilename);
    error(eid, ...
        'MONTAGE(D,[M N P]) is an obsolete syntax.\n%s', ...
        'Use multidimensional arrays to represent multiframe images.');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I, map] = getImagesFromFiles(fileNames)

if isempty(fileNames)
    eid = sprintf('Images:%s:invalidType', mfilename);
    msg = 'FILENAMES must be a cell array of strings.';
    error(eid, msg);
end
    
nframes = length(fileNames);
 
[img, map] = getImageFromFile(fileNames{1}, mfilename);
classImg = class(img);
sizeImg = size(img);

if length(sizeImg) > 2 && sizeImg(3) == 3
    nbands = 3;
else
    nbands = 1;
end
    
sizeImageArray = [sizeImg(1) sizeImg(2) nbands nframes]; 
if islogical(img)
    I = false(sizeImageArray);
else
    I = zeros(sizeImageArray, classImg);
end

I(:,:,:,1) = img;

for k = 2 : nframes
    [img,tempmap] = getImageFromFile(fileNames{k}, mfilename);
    
    if ~isequal(size(img),sizeImg)
        eid = sprintf('Images:%s:imagesNotSameSize', mfilename);
        error(eid, ...
            'FILENAMES must contain images that are the same size.');
    end
    
    if ~strcmp(class(img),classImg)
        eid = sprintf('Images:%s:imagesNotSameClass', mfilename);
        error(eid, ...
            'FILENAMES must contain images that have the same class.');
    end

    if isempty(map) && ~isempty(tempmap)
        map = tempmap;
    end
    I(:,:,:,k) = img;
end


warning('on', 'images:removing:function');

end
