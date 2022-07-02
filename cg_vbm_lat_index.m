function cg_vbm_lat_index(mask_raw)
% Optimized processing of asymmetry voxel based morphometry
%
% It is recommended to use your own gray matter template:
% First, run the script cg_create_template with all raw images from your sample
% to create a customized template.
%
% Run this script by selecting all raw (not normalized) images of your sample
% to be analyzed. This function creates the following files:
%
% output files:
% _______________________________________________________________________
%	w*.img		warped raw image
%	LI_Gw*.img:	gray matter of warped raw image
%	LI_Ww*.img:	white matter of warped raw image
%	LI_Cw*.img:	CSF of warped raw image
%	avgT1.img	average of all warped T1 images
% 	vbm_volumes.txt	absolute volumes of gray/white matter and csf
%
% If HMRF is applied all segmented images are indicated by _HMRF.
%
% use of customized prior images
% _______________________________________________________________________
% Depending on your dataset it is sometimes useful to use cutomized prior
% images resulting from cg_create_template.
% The use of own priors is recommended if the gray/white matter
% distribution of your sample differs from a control population
% (e.g. children, Alzheimer disease...).
%
% application of HMRF
% _______________________________________________________________________
% We can optionally use prior information by applying a Hidden Markov Random
% Field (HMRF) model. Using this model we introduce spatial constraints based
% on neighbouring voxels of a 3x3x3 cube. The center voxel has 26 neighbours
% and we can calculate MRF energy by counting the number of neighbours. The
% idea is to remove isolated voxels of one tissue class which are unlikely to
% be member of this tissue type. This procedure also closes holes in a cluster
% of connected voxels of one tissue type. In the resulting segmentation the
% noise level will be lowered.
%
% calculation of absolute gray/white matter and csf volumes
% _______________________________________________________________________
% A file "vbm_volumes.txt" is saved in the actual directory with the
% absolute volumes of gray/white matter and CSF (in ml). If this file already
% exists the new values will be appended. This file can be opened with
% Excel (or gnumeric) as text-file with tabs as delimiters. You can
% also use cg_read_vbm_volumes to read raw volumes from this file.
%
% mask raw lesions
% _______________________________________________________________________
% There is a hidden function to use masks for lesions or other regions
% which should not contribute to the estimation of normalization
% parameters. If you want to use this option call cg_vbm_lat_index(1).
%
% _______________________________________________________________________
% @(#)cg_vbm_lat_index.m	1.02 Christian Gaser 2006/07/31

global defaults

dseg				= defaults.segment;
dseg.write.wrt_cor		= 0;		% don't write bias corrected image
dseg.estimate.affreg.weight	= fullfile(spm('Dir'),'apriori','brainmask.mnc');	% always use brainmask
dseg.write.cleanup   		= 1;		% clean up partitions

dnrm				= defaults.normalise;
dnrm.write.vox			= [1 1 1];
dnrm.write.bb(1,3)		= -70;		% enlarge z-axis to include more cerebellum
dnrm.estimate.graphics		= 0;		% Crashes otherwise

if nargin == 0
    mask_raw = 0;
end

V   = spm_vol(spm_get(Inf,'*.IMAGE','Select raw (not normalized) images'));
VG0 = spm_vol(spm_get([0 1],'T1*.IMAGE','Select T1 template or press done for default'));
if isempty(VG0), VG0 = spm_vol(fullfile(spm('Dir'),'templates','T1.mnc')); end

PG1 = spm_get([0 1],'gray*.IMAGE','Select gray matter template for normalization or press done for default');
if isempty(PG1), PG1 = fullfile(spm('Dir'),'apriori','gray.mnc'); end

max_seg    = spm_input('Calculate',1,'m','gray matter|gray+white matter|all',[1 2 3],1);
own_prior  = spm_input('Use own prior images?','+1','yes|no',[1 0],2);
dseg.write.mrf = spm_input(['HMRF weighting?'],'+1','m',...
    ['large HMRF weighting (0.45)|medium HMRF weighting (0.3)|small HMRF weighting (0.15)|no HMRF (0)'],...
    [0.45 0.3 0.15 0], 2);

if own_prior
	dseg.estimate.priors = spm_get(3,'*.IMAGE','Select gray/white and CSF prior images');
end

if mask_raw
	Vmask = spm_vol(spm_get(size(V,1),'*.IMAGE','Select mask images in same order as raw images'));
end

% change default parameters for (non)linear normalization
dnrm.estimate.cutoff = get_cutoff(dnrm.estimate.cutoff);
if dnrm.estimate.cutoff ~= Inf
	dnrm.estimate.reg    = get_reg(dnrm.estimate.reg);
	dnrm.estimate.nits   = get_nits(dnrm.estimate.nits);
end

% make GM template symmetrical
PG1 = cg_make_symmetrical(PG1);
VG1 = spm_vol(PG1);

% make priors symmetrical
gm_sym  = cg_make_symmetrical(dseg.estimate.priors(1,:));
wm_sym  = cg_make_symmetrical(dseg.estimate.priors(2,:));
csf_sym = cg_make_symmetrical(dseg.estimate.priors(3,:));

dseg.estimate.priors = str2mat(gm_sym, wm_sym, csf_sym);

% create graphics window
fg=spm_figure('FindWin','Graphics');
if isempty(fg), fg=spm_figure('Create','Graphics'); end

% save gray/white matter and csf volumes in log-file
logname = 'vbm_volumes.txt';
if exist(logname,'file');
	fid = fopen(logname,'a');
else
	fid = fopen(logname,'w');
	fprintf(fid,'Volume [ml]\tgray\twhite\tcsf');
end

n = (0);
for i=1:length(V)
	[pth,nam,ext] = fileparts(V(i).fname);

	% first segmentation of raw image in original space
	dseg0 = dseg;
	dseg0.write.mrf = 0;   % don't use MRF in raw space because of huge memory demands
	VT  = spm_segment(V(i),VG0,dseg0);
	
	% calculate volumes of GM, WM and CSF
	factor = abs(det(VT(1).mat(1:3,1:3)))/1000;
	gray  = factor*sum(VT(1).dat(:))/255;
	white = factor*sum(VT(2).dat(:))/255;
	csf   = factor*sum(VT(3).dat(:))/255;
	
	% save volumes in log-file
	fprintf(fid,'\n%s\t%5.3f\t%5.3f\t%5.3f\t%5.3f',nam,gray,white,csf);

	VT = VT(1);

	% Estimate spatial normalisation parameters
	% by matching GM with a GM template
	matname = fullfile(pth,[nam '_sn.mat']);
	if mask_raw
		prm = spm_normalise(VG1,VT,matname,'',Vmask,dnrm.estimate);
	else
		prm = spm_normalise(VG1,VT,matname,'','',dnrm.estimate);	
	end
	clear VT
	
	% Create a spatially normalised version of the original image
	VN = spm_write_sn(V(i),prm,dnrm.write);

	% Segment the spatially normalised image
	VT = spm_segment(VN,eye(4),dseg);
	VN.fname = fullfile(pth,['w' nam ext]);
	spm_write_vol(VN,VN.dat);
	
	% make average of warped images
	n = n + double(VN.dat);

	% filenames for segmented images
	if dseg.write.mrf
	   nam2 = [nam '_HMRF'];
	else
	   nam2 = nam;
	end
	VT(1).fname = fullfile(pth,['LI_Gw' nam2 ext]);
	VT(2).fname = fullfile(pth,['LI_Ww' nam2 ext]);
	VT(3).fname = fullfile(pth,['LI_Cw' nam2 ext]);
		
	% CSF is better estimated using gray and white matter values
	t1 = ((double(VT(1).dat) + double(VT(2).dat) + double(VT(3).dat)) > 0.5*255);
	t2 = (1 - (double(VT(1).dat) + double(VT(2).dat))/255);
	VT(3).dat = uint8(round((t1.*t2)*255));

    % calculate LI=2*(R-L)/(R+L)
	% and finally save segmented images
    for j=1:max_seg
        VT(j).dat = double(VT(j).dat);
        VT(j).dat = 2*(VT(j).dat-VT(j).dat(end:-1:1,:,:))./(VT(j).dat+VT(j).dat(end:-1:1,:,:)+eps);
        VT(j).dim(4) = 4;
		spm_write_vol(VT(j),VT(j).dat);
		VT(j).dat = [];
    end
    
	clear VT
end;

fclose(fid);

% save average image
n = n/length(V);
VN.fname   = 'avgT1.img';
VN.descrip = 'Average of warped T1 images';
spm_write_vol(VN,n);

% remove temporary symmetrical priors and template
remove_images(PG1);
remove_images(gm_sym);
remove_images(wm_sym);
remove_images(csf_sym);

return

%_______________________________________________________________________
function remove_images(P)
% remove temporary images

[pth, nam, ext, ver] = fileparts(P);
spm_unlink(fullfile(pth,[nam '.img']),fullfile(pth,[nam '.hdr']));

return
%_______________________________________________________________________

%_______________________________________________________________________
function nits = get_nits(nits)
% Get number of nonlinear iterations
% FORMAT nits = get_nits(nits)

if prod(nits) > 0,
        tmp2 = [1 3 5 8 12 16];
        tmp = find(tmp2 == nits);
        if isempty(tmp) tmp = length(tmp2); end;
        nits = spm_input(['# Nonlinear Iterations?'],'+1','m',...
        	['1  nonlinear iteration |3  nonlinear iterations'...
        	'|5  nonlinear iterations|8  nonlinear iterations'...
        	'|12 nonlinear iterations|16 nonlinear iterations'],tmp2, tmp);
else, nits = 0; end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function cutoff = get_cutoff(cutoff)
% Get cutoff frequency of DCT bases
% FORMAT cutoff = get_cutoff(cutoff)

tmp2 = [15 20 25 30 35 40 45 50 60 70 80 Inf];
tmp = find(tmp2 == cutoff);
if isempty(tmp) tmp = length(tmp2); end;
cutoff = spm_input('Cutoff spatial normalization','+1','m',...
        ['15mm cutoff|20mm cutoff|25mm cutoff|30mm cutoff|'...
	 '35mm cutoff|40mm cutoff|45mm cutoff|50mm cutoff|'...
	 '60mm cutoff|70mm cutoff|80mm cutoff|Affine only'], tmp2, tmp);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function reg = get_reg(reg)
% Get amount of regularisation
% FORMAT reg = get_reg(reg)

tmp2 = [100 10 1 0.1 0.01];
tmp = find(tmp2 == reg);
if isempty(tmp) tmp = length(tmp2); end;
reg = spm_input('Nonlinear Regularization','+1','m',...
        ['Extremely Heavy regularization (100)|Heavy regularization (10)|'...
         'Medium regularization (1)|Light regularization (0.1)|'...
         'Very Light regularization (0.01)'], tmp2, tmp);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function PO = cg_make_symmetrical(P)
% Make a symmetrical copy of input images (for symmetrical templates and priors)
% New files are saved as temporary file

PO = [tempname '.img'];
spm_imcalc_ui(P,PO,'(i1+flipud(i1))/2',{0 0 2 1});

return

