% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function s = IPEDtxtReader(txt_fname)
% Reads IPED-text file and loads diffusion params. The IPED-text is a
% simple plain text file which describes different parameters for DWIs
% acquired. The file follows a table type format, where each row of the file
% describes parameter for one set of DWIs in either of the following
% formats:
%  
%   4D-nifti-filename.nii.gz  PED  bvec-filename  bval-filename
%
%   4D-nifti-filename.nii.gz  PED  bmat-filename
%
%   4D-nifti-filename.nii.gz  PED  bvec-filename  bval-filename bmat-filename
%
% In addition, any text coming after a # (hash) is ignored for rest of the
% line.
%
% Filenames can be specified either with their full-path or without any path (for the later case,
% files are assumed to be in the same directory as the IPED-text file. Relative paths (eg:
% subdirectory/file.nii.gz) are NOT supported. When the  filenames/paths contain spaces then they
% must be enclosed in double quotes (eg: "subject23 trial2.nii.gz"). Different columns of the
% textfile can also be delimited by comma, like that in standard csv files
% (in case of csv format, filenames with spaces can be directly specified and
% must NOT use double quotes). 
% 
% Valid options for PED (phase encoding direction) are x, x-, y, y-, z and
% z-. x direction indicates towards the right side of the subject, while x-
% increases towards the left side of the subject. Similarly, y and y- are
% along the anterior-posterior direction of the subject, and z & z- are
% along the inferior-superior direction.
%


fid = fopen(txt_fname, 'r');
if fid<0
   error('Could not open/find the IPED-text file: %s', txt_fname);
end
p_dir = fileparts(txt_fname);

% get first row
row1 = stripComments(fgetl(fid));
while ischar(row1) && isempty(row1)
   row1 = stripComments(fgetl(fid));
end
if ~ischar(row1)
   error('Input IPED-text file is empty! \nPlease add IPED data information to the file: %s', txt_fname);
end

% check for csv file
if isempty(strfind(row1, ','))
   csv_input = false;
else
   csv_input = true;
end


% get number of cols etc.
temp = line2cols(row1, csv_input);
nCols = length(temp);
if nCols<3
   error('IPED-text file must have atleast three columns: %s', txt_fname);
   
elseif nCols>5
   error('IPED-text file can have atmost five columns: %s', txt_fname);
   
elseif nCols==3
   fprintf('\nThree columns found in IPED-text file. \nLast column will be used as b-Matrices file.')
   
elseif nCols==4
   fprintf('\nFour columns found in IPED-text file. \nLast two column will be used as bvec-file followed by bval-file.')
   
elseif nCols==5
   fprintf('\nFive columns found in IPED-text file. \nLast column will be used as b-Matrices file.')
end


% read parameters line by line
row_str = row1;
s = [];
while ischar(row_str)
   c = line2cols(row_str, csv_input);
   if ~isempty(c)
      if length(c)==nCols
         temp = cols2PEDstruct(c, nCols, p_dir);
         if isempty(s)
            s = temp;
         else
            s = cat(1, s, temp);
         end
      else
         error('Number of columns in IPED-txt file is not consistent: %s', txt_fname)
      end
   end
   row_str = fgetl(fid);
end

fclose(fid);

end

function s = cols2PEDstruct(c, nCols, p_dir)

% add p_dir, when applicable
for k = [1 3:nCols]
   if isempty(fileparts(c{k}))
      c{k} = fullfile(p_dir, c{k});
   end
end

if nCols==3
   s = checkDiffusionEncodingScheme(c{3});
   s.bvec_file = '';
   s.bval_file = '';
   s.bmat_file = c{3};
   
elseif nCols==4
   s = checkDiffusionEncodingScheme(c{3}, c{4});
   s.bvec_file = c{3};
   s.bval_file = c{4};
   s.bmat_file = '';
   
elseif nCols==5
   s = checkDiffusionEncodingScheme(c{5}); % ignore bvec-bval when bmat is present
   s.bvec_file = c{3};
   s.bval_file = c{4};
   s.bmat_file = c{5};
   
end

s.dwi_file = c{1};
if ismember(lower(c{2}), {'x' 'x-' 'y' 'y-' 'z' 'z-'})
   s.PED = lower(c{2});   
else
   error('\nUnidentified PED-direction: %s \nValid options for PED-direction are x, x-, y, y-, z and z-', c{2});
end

end

function str = stripComments(str)
str = regexprep(str, '#.*', '');
str = strtrim(str);
end

function C = line2cols(str, csv_input)
% separate a line into columns

str = stripComments(str);
if isempty(str)
   C = {}; 
   return;
end

if csv_input
   C = textscan(str, '%s', 'Delimiter', ',');
   C = C{1};
   
else
   [matchstart, matchend] = regexp(str, '"[^"]*"');
   
   if isempty(matchstart)
      C = textscan(str, '%s'); % delimit using white space
      C = C{1};
      
   else % double quotes found
      str_len = length(str);      
      
      % add fake-end points for simpler for-loop in next part
      matchstart = [0, matchstart, str_len+1];
      matchend = [0, matchend, str_len+1];
      last_match = length(matchstart);
      
      C = {};
      for k = 2:length(matchstart)
         % string before double quotes - delimit by space
         ind1 = max(1, matchend(k-1)+1);
         ind2 = min(matchstart(k)-1, str_len);
         temp = line2cols(str(ind1:ind2), false);
         if ~isempty(temp)
            C = cat(1, C, temp);
         end
         
         % string in double quotes - concat as it is
         if k<last_match % skip last (fake) match
            ind1 = matchstart(k)+1;
            ind2 = matchend(k)-1;
            C = cat(1, C, str(ind1:ind2));
         end
      end
   end
   
end


end
