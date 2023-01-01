% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [ fileName ] = extract_header(niiData, fileName)
%Extract the NIFTI header to text file or command prompt.
%This version is intended to work with non-standard NIFTI type files, like
%.eig format. (Should work without any problem with standard NIFTI files)
%
%   niiData - nii structure or .nii/.nii.gz file.
%   fileName - (Optional) Name (string) of output text file. If the string
%   is empty [] or equals 'no_file' then the text is not written to file
%   but is displayed on standard output (the screen / command prompt).
%
%   Requires NIFTI toolbox.
%

%fprintf('\n==== Writing header text dump ====\n');

if ischar(niiData)
	data.hdr = load_untouch_header_only_non_standard_gz(niiData);
	data.fileprefix = remove_extension(niiData);
else
	data = niiData;
end

if ~exist('fileName', 'var')
	fileName = [data.fileprefix '.hdr.txt'];
	fileID = fopen(fileName, 'w');
elseif isempty(filename) || strcmpi(fileName, 'no_file')
	fileID = 1;
else
	fileID = fopen(fileName, 'w');
end

[~, fname] = file_base_name(data.fileprefix);

if isfield(data, 'original')
	fprintf('''Original'' field found. Will be followed by main header.\n')
	fprintf(fileID, '\n  Note: Check ''original'' field followed by main header section \n');
end


if isfield(data, 'hdr')
	dump_field(data, 'hdr', fileID, fname);
else
	fprintf('No header field found.\n')
end


if isfield(data, 'original')
	fprintf(fileID, '\n  - - - - ''original'' field - - - - \n');
	dump_field(data, 'original', fileID, fname);
end

if fileID~=1
	fclose(fileID);
end

end


function dump_field (S, f, fid, pStr)
%   writes the content of field to text file
%
%   S - structure
%   f - name of field in char
%   fid - file identifier
%   pStr - parent string
%

sField = getfield(S, f);

if isstruct(sField)
	pStr = [pStr '.' f];
	subFields = fieldnames(sField);
	subFieldNums = size(subFields);
	
	fprintf(fid, '\n%s\n',pStr);
	for k = 1:subFieldNums(1)
		dump_field(sField, char(subFields(k)), fid, pStr);
	end
	fprintf(fid, '\n\n',f);
	
else
	fprintf(fid, '\n    %-14s    ',f);
	switch class(sField)
		case 'char'
			fprintf(fid, ' %s',sField);
		case {'double', 'single'}
			fprintf(fid, '%- 1.5g    ',sField);
		otherwise
			fprintf(fid, 'Warning: Unknow class:  %s ====\n',class(sField));
	end
	
end

end

