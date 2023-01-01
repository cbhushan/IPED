% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function bstfnames = bdpWriteBST(opts)
% Writes BST files based on options and inputs

blank_struct = struct(...
   'path',	[],...
   'volume', [],...
   'label', [],...
   'labeltext', [],...
   'overlay', [], ...
   'overlay2', [],...
   'mask', [], ...
   'surface', [], ...
   'curveset', [], ...
   'tractset', [], ...
   'SHC', [], ...
   'connectivity_matrix', []);
bstfnames = {};

bdpPrintSectionHeader('Writing BST files');

fprintf(bdp_linewrap(['BST files can be opened in BrainSuite for conveniently loading ' ...
   'most data/files required for Tractography and Connectivity analysis.']));
fprintf('\nBDP is checking for required files...');

% find DWI filename
if (opts.fieldmap_distortion_correction || opts.registration_distortion_correction) ...
      && exist([opts.file_base_name opts.dwi_fname_suffix '.RAS' opts.dwi_corrected_suffix '.nii.gz'], 'file')
   dwi_file = [opts.file_base_name opts.dwi_fname_suffix '.RAS' opts.dwi_corrected_suffix '.nii.gz'];
   
elseif ~(opts.fieldmap_distortion_correction || opts.registration_distortion_correction) ...
      && exist([opts.file_base_name opts.dwi_fname_suffix '.RAS.nii.gz'], 'file')
   dwi_file = [opts.file_base_name opts.dwi_fname_suffix '.RAS.nii.gz'];
   
else
   err_msg = ['BDP can not find the diffusion file written by it. Please make sure that you ran '...
      'BDP on fileprefix: ' escape_filename(opts.file_base_name) '\nIf you did run BDP previously, '...
      'please make sure that all flags, except --generate-only-stats, are exactly same as previous '...
      'BDP run. You can find the command used in BDP summary file (<fileprefix>.BDPSummary.txt).'];
   error('BDP:ExpectedFileNotFound', bdp_linewrap(err_msg));
end
dwi_filebase = fileBaseName(dwi_file);

% find T1 filename
struct_fileprefix = [];
if exist([opts.bfc_file_base '.nii.gz'], 'file')==2
   [p, n, e] = fileparts([opts.bfc_file_base '.nii.gz']);
   struct_file = [n e];
   if ~isempty(opts.diffusion_coord_output_folder)
      [p, n, e] = fileparts([opts.bfc_file_base opts.Diffusion_coord_suffix '.nii.gz']);
      struct_file_D = [n e];
   end
   
elseif exist([opts.bfc_file_base '.nii'], 'file')==2
   [p, n, e] = fileparts([opts.bfc_file_base '.nii']);
   struct_file = [n e];
   if ~isempty(opts.diffusion_coord_output_folder)
      [p, n, e] = fileparts([opts.bfc_file_base opts.Diffusion_coord_suffix '.nii.gz']);
      struct_file_D = [n e];
   end
   
else
   [p, n, e] = fileparts(opts.bfc_file);
   struct_file = [n e];
   if ~isempty(opts.diffusion_coord_output_folder)
      [p, n, e] = fileparts([opts.bfc_file_base '.bfc' opts.Diffusion_coord_suffix '.nii.gz']);
      struct_file_D = [n e];
   end
end

if ~isequal(opts.bfc_file_base, opts.file_base_name) % output-subdir
   struct_fileprefix = '../';
   struct_file = [struct_fileprefix struct_file];
end

% Look for masks and SVReg outputs
if exist([opts.bfc_file_base '.mask.nii.gz'], 'file')==2
   [p, n, e] = fileparts([opts.bfc_file_base '.mask.nii.gz']);
   T1_mask = [struct_fileprefix n e];
   
   if ~isempty(opts.diffusion_coord_output_folder)
      [p, n, e] = fileparts([opts.bfc_file_base opts.Diffusion_coord_suffix '.mask.nii.gz']);
      T1_mask_D = [n e];
   end
else
   T1_mask = [];
   T1_mask_D = [];
   msg = {'\n\n', ['BDP could not find the default brain mask (.mask.nii.gz) file. Use of brain mask reduces '...
      'the run-time for tractography significantly and avoids spurious fiber tracks. Brain mask can '...
      'be generated and hand-edited in BrainSuite.'], '\n'};
   fprintf(bdp_linewrap(msg));
end

if exist([opts.bfc_file_base '.right.pial.cortex.svreg.dfs'], 'file')==2
   [p, n, e] = fileparts([opts.bfc_file_base '.right.pial.cortex.svreg.dfs']);
   R_surf = [struct_fileprefix n e];  
else
   R_surf = [];   
end

if exist([opts.bfc_file_base '.svreg.corr.label.nii.gz'], 'file')==2
   [p, n, e] = fileparts([opts.bfc_file_base '.svreg.corr.label.nii.gz']);
   svreg_label = [struct_fileprefix n e];
   if ~isempty(opts.diffusion_coord_output_folder)
      [p, n, e] = fileparts([opts.bfc_file_base '.svreg.corr' opts.Diffusion_coord_suffix '.label.nii.gz']);
      svreg_label_D = [n e];
   end
   
elseif exist([opts.bfc_file_base '.svreg.label.nii.gz'], 'file')==2
   [p, n, e] = fileparts([opts.bfc_file_base '.svreg.label.nii.gz']);
   svreg_label = [struct_fileprefix n e];
   if ~isempty(opts.diffusion_coord_output_folder)
      [p, n, e] = fileparts([opts.bfc_file_base '.svreg' opts.Diffusion_coord_suffix '.label.nii.gz']);
      svreg_label_D = [n e];
   end
   
else
   svreg_label = [];
   svreg_label_D = [];
   xml_label = [];
   xml_label_D = [];
   msg = {'\n\n', ['BDP could not find the default label (.label.nii.gz) file written by SVReg. Tractography '...
      'can be run without label file. However, label file is required for connectivity analysis. Label file '...
      'can be generated by running Surface/Volume registration in BrainSuite.'], '\n'};
   fprintf(bdp_linewrap(msg));
end

if ~isempty(svreg_label) % makes sense only when labels are found
   if exist(fullfile(fileparts(opts.bfc_file_base), 'brainsuite_labeldescription.xml'), 'file')==2
      xml_label = [struct_fileprefix 'brainsuite_labeldescription.xml'];
      if ~isempty(opts.diffusion_coord_output_folder)
         xml_label_D = [struct_fileprefix '../brainsuite_labeldescription.xml'];
      end
      
   elseif exist(fullfile(fileparts(opts.bfc_file_base), 'brainsuite_labeltext.xml'), 'file')==2
      xml_label = [struct_fileprefix 'brainsuite_labeltext.xml'];
      if ~isempty(opts.diffusion_coord_output_folder)
         xml_label_D = [struct_fileprefix '../brainsuite_labeltext.xml'];
      end
      
   else
      xml_label = [];
      xml_label_D = [];
      msg = {'\n\n', ['BDP could not find the default label description (.xml) file written by SVReg. '...
         'It will use default label description (and colors) for ROIs.'], '\n'};
      fprintf(bdp_linewrap(msg));
   end
end

bcnt = 1;
% Tensor estimation
if opts.estimate_tensor
   bst_fname = [opts.file_base_name '.tensor' opts.mprage_coord_suffix '.bst'];
   bsts = blank_struct;
   bsts.volume = [dwi_filebase opts.mprage_coord_suffix '.eig.nii.gz'];
   bsts.overlay = struct_file;
   bsts.mask = T1_mask;
   bsts.label = svreg_label;
   bsts.labeltext = xml_label;
   bsts.surface = R_surf;
   writeBST(bsts, bst_fname)
   bstfnames{bcnt} = bst_fname;
   bcnt = bcnt +1;
   
   if opts.diffusion_coord_outputs
      bst_fname = fullfile(opts.diffusion_coord_output_folder, ...
         [fileBaseName(opts.bfc_file_base) '.tensor' opts.Diffusion_coord_suffix '.bst']);
      bsts = blank_struct;
      bsts.volume = [dwi_filebase '.eig.nii.gz'];
      bsts.overlay = struct_file_D;
      bsts.mask = T1_mask_D;
      bsts.label = svreg_label_D;
      bsts.labeltext = xml_label_D;      
      writeBST(bsts, bst_fname)
      bstfnames{bcnt} = bst_fname;
      bcnt = bcnt +1;
   end
end

% FRACT
if opts.estimate_odf_FRACT
   bst_fname = [opts.file_base_name '.FRACT' opts.mprage_coord_suffix '.bst'];
   bsts = blank_struct;
   bsts.SHC = fullfile('FRACT', [dwi_filebase '.SH.FRACT' opts.mprage_coord_suffix '.odf']);
   bsts.volume = struct_file;
   bsts.overlay = [dwi_filebase '.FA.color' opts.mprage_coord_suffix '.nii.gz'];
   bsts.mask = T1_mask;
   bsts.label = svreg_label;
   bsts.labeltext = xml_label;
   bsts.surface = R_surf;
   writeBST(bsts, bst_fname)
   bstfnames{bcnt} = bst_fname;
   bcnt = bcnt +1;
   
   if opts.diffusion_coord_outputs
      bst_fname = fullfile(opts.diffusion_coord_output_folder, ...
         [fileBaseName(opts.bfc_file_base) '.FRACT' opts.Diffusion_coord_suffix '.bst']);
      bsts = blank_struct;
      bsts.SHC = fullfile('FRACT', [dwi_filebase '.SH.FRACT.odf']);
      bsts.volume = struct_file_D;
      bsts.overlay = [dwi_filebase '.eig.nii.gz'];
      bsts.mask = T1_mask_D;
      bsts.label = svreg_label_D;
      bsts.labeltext = xml_label_D;      
      writeBST(bsts, bst_fname)
      bstfnames{bcnt} = bst_fname;
      bcnt = bcnt +1;
   end
end

% FRT
if opts.estimate_odf_FRT
   bst_fname = [opts.file_base_name '.FRT' opts.mprage_coord_suffix '.bst'];
   bsts = blank_struct;
   bsts.SHC = fullfile('FRT', [dwi_filebase '.SH.FRT' opts.mprage_coord_suffix '.odf']);
   bsts.volume = struct_file;
   bsts.overlay = [dwi_filebase '.FA.color' opts.mprage_coord_suffix '.nii.gz'];
   bsts.mask = T1_mask;
   bsts.label = svreg_label;
   bsts.labeltext = xml_label;
   bsts.surface = R_surf;
   writeBST(bsts, bst_fname)
   bstfnames{bcnt} = bst_fname;
   bcnt = bcnt +1;
   
   if opts.diffusion_coord_outputs
      bst_fname = fullfile(opts.diffusion_coord_output_folder, ...
         [fileBaseName(opts.bfc_file_base) '.FRT' opts.Diffusion_coord_suffix '.bst']);
      bsts = blank_struct;
      bsts.SHC = fullfile('FRT', [dwi_filebase '.SH.FRT.odf']);
      bsts.volume = struct_file_D;
      bsts.overlay = [dwi_filebase '.eig.nii.gz'];
      bsts.mask = T1_mask_D;
      bsts.label = svreg_label_D;
      bsts.labeltext = xml_label_D;      
      writeBST(bsts, bst_fname)
      bstfnames{bcnt} = bst_fname;
      bcnt = bcnt +1;
   end
end

% 3DSHORE
if opts.estimate_odf_3DSHORE
   bst_fname = [opts.file_base_name '.3DSHORE' opts.mprage_coord_suffix '.bst'];
   bsts = blank_struct;
   bsts.SHC = fullfile('3DSHORE', [dwi_filebase '.SH.3DSHORE' opts.mprage_coord_suffix '.odf']);
   bsts.volume = struct_file;
   %bsts.overlay = [dwi_filebase '.FA.color' opts.mprage_coord_suffix '.nii.gz'];
   bsts.mask = T1_mask;
   bsts.label = svreg_label;
   bsts.labeltext = xml_label;
   bsts.surface = R_surf;
   writeBST(bsts, bst_fname)
   bstfnames{bcnt} = bst_fname;
   bcnt = bcnt +1;
   
   if opts.diffusion_coord_outputs
      bst_fname = fullfile(opts.diffusion_coord_output_folder, ...
         [fileBaseName(opts.bfc_file_base) '.3DSHORE' opts.Diffusion_coord_suffix '.bst']);
      bsts = blank_struct;
      bsts.SHC = fullfile('3DSHORE', [dwi_filebase '.SH.3DSHORE.odf']);
      bsts.volume = struct_file_D;
      %bsts.overlay = [dwi_filebase '.eig.nii.gz'];
      bsts.mask = T1_mask_D;
      bsts.label = svreg_label_D;
      bsts.labeltext = xml_label_D;      
      writeBST(bsts, bst_fname)
      bstfnames{bcnt} = bst_fname;
      bcnt = bcnt +1;
   end
end   

% GQI
if opts.estimate_odf_GQI
   bst_fname = [opts.file_base_name '.GQI' opts.mprage_coord_suffix '.bst'];
   bsts = blank_struct;
   bsts.SHC = fullfile('GQI', [dwi_filebase '.SH.GQI' opts.mprage_coord_suffix '.odf']);
   bsts.volume = struct_file;
   %bsts.overlay = [dwi_filebase '.FA.color' opts.mprage_coord_suffix '.nii.gz'];
   bsts.mask = T1_mask;
   bsts.label = svreg_label;
   bsts.labeltext = xml_label;
   bsts.surface = R_surf;
   writeBST(bsts, bst_fname)
   bstfnames{bcnt} = bst_fname;
   bcnt = bcnt +1;
   
   if opts.diffusion_coord_outputs
      bst_fname = fullfile(opts.diffusion_coord_output_folder, ...
         [fileBaseName(opts.bfc_file_base) '.GQI' opts.Diffusion_coord_suffix '.bst']);
      bsts = blank_struct;
      bsts.SHC = fullfile('GQI', [dwi_filebase '.SH.GQI.odf']);
      bsts.volume = struct_file_D;
      %bsts.overlay = [dwi_filebase '.eig.nii.gz'];
      bsts.mask = T1_mask_D;
      bsts.label = svreg_label_D;
      bsts.labeltext = xml_label_D;      
      writeBST(bsts, bst_fname)
      bstfnames{bcnt} = bst_fname;
      bcnt = bcnt +1;
   end   
end

if ~isempty(bstfnames)
   fprintf('\nBDP saved following BST files:')
   for bcnt = 1:length(bstfnames)
      fprintf('\n\t%s', bstfnames{bcnt});
   end
   fprintf('\n');
end
end

function writeBST(file_struct, filename)
fid = fopen(filename, 'w');
if fid<0
   error('BDP:FileOpenFailed', 'BDP could not open file for writing BST file: %s', escape_filename(filename));
end

if ~isempty(file_struct.path)
   fprintf(fid, 'path %s\n', file_struct.path);
end

if ~isempty(file_struct.SHC)
   fprintf(fid, 'SHC %s\n', file_struct.SHC);
end

if ~isempty(file_struct.volume)
   fprintf(fid, 'volume %s\n', file_struct.volume);
end

if ~isempty(file_struct.label)
   fprintf(fid, 'label %s\n', file_struct.label);
end

if ~isempty(file_struct.labeltext)
   fprintf(fid, 'labeltext %s\n', file_struct.labeltext);
end

if ~isempty(file_struct.overlay)
   fprintf(fid, 'overlay %s\n', file_struct.overlay);
end

if ~isempty(file_struct.overlay2)
   fprintf(fid, 'overlay2 %s\n', file_struct.overlay2);
end

if ~isempty(file_struct.mask)
   fprintf(fid, 'mask %s\n', file_struct.mask);
end

if ~isempty(file_struct.surface)
   fprintf(fid, 'surface %s\n', file_struct.surface);
end

if ~isempty(file_struct.curveset)
   fprintf(fid, 'curveset %s\n', file_struct.curveset);
end

if ~isempty(file_struct.tractset)
   fprintf(fid, 'tractset %s\n', file_struct.tractset);
end

if ~isempty(file_struct.connectivity_matrix)
   fprintf(fid, 'connectivity_matrix %s\n', file_struct.connectivity_matrix);
end

fclose(fid);

end




