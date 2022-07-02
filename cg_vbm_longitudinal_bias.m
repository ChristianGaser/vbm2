function cg_vbm_longitudinal(bb)
% Optimized processing of longitudinal data with voxel based morphometry
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
%	Gw*.img:	gray matter of warped raw image
%	Ww*.img:	white matter of warped raw image
%	Cw*.img:	CSF of warped raw image
%	avgT1.img	average of all warped T1 images
% 	vbm_volumes.txt	absolute volumes of gray/white matter and csf
%
% If HMRF is applied all segmented images are indicated by _HMRF.
%
% use of customized prior images
% _______________________________________________________________________
% Depending on your dataset it is also useful to use cutomized prior
% images resulting from cg_create_template.
% The use of own priors is recommended if the gray/white matter
% distribution of your sample differs from a control population
% (e.g. children, Alzheimer disease...) and if you have enough
% subjects to represent a population.
%
% application of HMRF
% _______________________________________________________________________
% We can optionally use prior information by applying a Hidden Markov Random
% Field (HMRF) model. Using this model we introduce spatial constraints based
% on neighbouring voxels of a 3x3x3 cube. The center voxel has 26 neighbours
% and we can calculate MRF energy by counting the number of neighbours. The
% idea is to remove isolated voxels of one tissue class which are unlikely to
% be member of this tissue type. This procedure also closes holes in a clsuter
% of connected voxels of one tissue type. In the resulting segmentation the
% noise level will be lowered.
%
% calculation of absolute gray/white matter and csf volumes
% _______________________________________________________________________
% A file "vbm_volumes.txt" is saved in the actual directory with the
% absolute volumes of gray/white matter and CSF (in ml). If this file
% exists the new values will be appended. This file can be opened with
% Excel (or gnumeric) as text-file with tabs as delimiters.
%
% additional bias correction between time points
% _______________________________________________________________________
% Because longitudinal scans are often acquired with a large gap between the
% time points this may cause different distributions of intensity nonuni-
% formities. Although each scan will be corrected independently using bias 
% correction subtle but systematical differences in the bias field remain and may
% influence results of segmentation. To minimize these intensity differences
% between the time points we approximate the difference bias field by using 
% the intracraniel parts of the difference image which is smoothed with a
% very large gaussian kernel of around 30 mm. However, more exact methods might
% exist for this purpose.
%f
%_______________________________________________________________________
% @(#)cg_vbm_longitudinal.m	1.06 Christian Gaser 2006/11/02

global defaults

dseg				= defaults.segment;
dseg.write.wrt_cor		= 0;                            % don't write bias corrected image
dseg.estimate.affreg.weight	= fullfile(spm('Dir'),'apriori','brainmask.mnc');	% always use brainmask
dseg.write.cleanup   		= 1;                         % clean up partitions

dnrm				= defaults.normalise;
dnrm.write.vox			= [1 1 1];
dnrm.write.bb(1,3)            = -70;                    % enlarge z-axis to include more cerebellum
dnrm.estimate.graphics		= 0;		% Crashes otherwise

% hidden function to change bounding box
if nargin == 1
	dnrm.write.bb =		[-84 -118 -70; 83 81 89];		% dimension 168 x 200 x 160
end

T1_template			= fullfile(spm('Dir'),'templates','T1.mnc');	% T1 template for first normalization

linfun = inline('fprintf([''%-60s%s''],x,[char(sprintf(''\b'')*ones(1,60))])','x');

SPMid = spm('FnBanner',mfilename,'1.05');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','VBM longitudinal');
spm_help('!ContextHelp',mfilename);

% select images for each subject
don = 0;
for i = 1:1000,
	PF = spm_get(Inf,'IMAGE',...
		['Select all time points for subj ' num2str(i)]);
	if isempty(PF), don = 1; break; end;
	VF{i} = spm_vol(PF);
end;

VG0 = spm_vol(T1_template);
GM = spm_get([0 1],'gray*.IMAGE','Select GM template or press done for default');
if isempty(GM)
    GM = fullfile(spm('Dir'),'apriori','gray.mnc');
end
VG1 = spm_vol(GM);

Vm	= spm_vol(fullfile(spm('Dir'),'apriori','brainmask.mnc'));

max_seg    = spm_input('Calculate',1,'m','gray matter|gray+white matter|all',[1 2 3],1);
own_prior  = spm_input('Use own prior images?','+1','yes|no',[1 0],2);
dseg.write.mrf = spm_input(['HMRF weighting?'],'+1','m',...
    ['large HMRF weighting (0.45)|medium HMRF weighting (0.3)|small HMRF weighting (0.15)|no HMRF (0)'],...
    [0.45 0.3 0.15 0], 2);

% smoothing of difference image to remove inhomogeneities
FWHM_diff = spm_input('FWHM for add. bias corr. between TPs (0 for no correction)','+1','e',50);

if own_prior
	dseg.estimate.priors = spm_get(3,'*.IMAGE','Select gray/white and CSF prior images');
end

% change default parameters for (non)linear normalization
dnrm.estimate.cutoff = get_cutoff(dnrm.estimate.cutoff);
if dnrm.estimate.cutoff ~= Inf
	dnrm.estimate.reg    = get_reg(dnrm.estimate.reg);
	dnrm.estimate.nits   = get_nits(dnrm.estimate.nits);
end

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

n = 0;
for i = 1:length(VF),
	V = VF{i};
	
	[pth1,nam1,ext1] = fileparts(V(1).fname);

	% first segmentation of first raw image in original space
	linfun(['Segmentation of raw image: ' nam1]);
	dseg0 = dseg;
	dseg0.write.mrf = 0;   % don't use MRF in raw space because of huge memory demands
	VT = spm_segment(V(1),VG0,dseg0);

	% calculate volumes of GM, WM and CSF for first image
	factor = abs(det(VT(1).mat(1:3,1:3)))/1000;
	gray  = factor*sum(VT(1).dat(:))/255;
	white = factor*sum(VT(2).dat(:))/255;
	csf   = factor*sum(VT(3).dat(:))/255;
	
	% save volumes in log-file
	fprintf(fid,'\n%s\t%5.3f\t%5.3f\t%5.3f\t%5.3f',nam1,gray,white,csf);
	
	% compute binary mask for realignment using GM+WM+CSF of first image
	mask = (double(VT(1).dat)+double(VT(2).dat)+double(VT(3).dat))>128;
	VW = VT(1);
	VW.fname = fullfile(pth1,['brainmask_' nam1 ext1]);
	spm_write_vol(VW,mask);
	
	% flags for realignment: use weighting mask and high quality
	flagsC = struct('quality',1,'fwhm',4,'rtm',0,'PW',VW);
	% realign T1 images
	linfun(['Realign: ' nam1]);
	V = spm_realign(V,flagsC);
	
	% Estimate spatial normalisation parameters
	% by matching first GM image with GM template
	matname = fullfile(pth1,[nam1 '_sn.mat']);
	linfun(['Normalize: ' nam1]);
	prm = spm_normalise(VG1,VT(1),matname,'','',dnrm.estimate);

	% Create a spatially normalised version of the original images
	for j = 1:length(V)
		[pth,nam,ext] = fileparts(V(j).fname);

		VN = spm_write_sn(V(j),prm,dnrm.write);
		% Segment the spatially normalised image
		if (FWHM_diff > 0)
		  if (j == 1)
            img0 = double(VN.dat);
		  else
    		linfun(['Bias correct difference: ' nam]);
		    % calculate masked difference to first image
		    mdiff = (double(VN.dat)-img0);
		    mdiff = mdiff - mean(mdiff(:));
		    % smooth masked difference image
		    smdiff = zeros(size(mdiff));
		    spm_smooth(mdiff,smdiff,FWHM_diff);
		    % and subtract estimated smoothed difference
		    VN.dat = single(round(double(VN.dat) - smdiff));
		  end
		end

		linfun(['Segment: ' nam]);
		VT = spm_segment(VN,eye(4),dseg);

		VN.fname = fullfile(pth,['w' nam ext]);
		spm_write_vol(VN,VN.dat);
		% make average of first warped images
		if j==1, n = n + double(VN.dat); end

    	% filenames for segmented images
    	if dseg.write.mrf
	       nam2 = [nam '_HMRF'];
        else
	       nam2 = nam;
        end
        VT(1).fname = fullfile(pth,['Gw' nam2 ext]);
        VT(2).fname = fullfile(pth,['Ww' nam2 ext]);
        VT(3).fname = fullfile(pth,['Cw' nam2 ext]);
		
		% CSF is better estimated using gray and white matter values
		t1 = ((double(VT(1).dat) + double(VT(2).dat) + double(VT(3).dat)) > 0.5*255);
		t2 = (1 - (double(VT(1).dat) + double(VT(2).dat))/255);
		VT(3).dat = uint8(round((t1.*t2)*255));
		clear t1 t2
	    
		% and finally save segmented images
		for k = 1:max_seg
			spm_write_vol(VT(k),double(VT(k).dat)/255);
		end
	end
	clear VT
	
end;

% close log file
fclose(fid);

% save average image
n = n/length(V);
VN.fname   = 'avgT1.img';
VN.descrip = 'Average of warped T1 images';
spm_write_vol(VN,n);
return

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

function V = smooth_vol(V,Y,s)
% smooth volume with gaussian kernel
% based on spm_smooth

% compute parameters for spm_conv_vol
vx = sqrt(sum(V.mat(1:3,1:3).^2));
s  = s./vx;					% voxel anisotropy
s  = max(s,ones(size(s)));			% lower bound on FWHM
s  = s/sqrt(8*log(2));				% FWHM -> Gaussian parameter

x  = round(6*s(1)); x = [-x:x];
y  = round(6*s(2)); y = [-y:y];
z  = round(6*s(3)); z = [-z:z];
x  = exp(-(x).^2/(2*(s(1)).^2)); 
y  = exp(-(y).^2/(2*(s(2)).^2)); 
z  = exp(-(z).^2/(2*(s(3)).^2));
x  = x/sum(x);
y  = y/sum(y);
z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;

spm_conv_vol(Y,Y,x,y,z,-[i,j,k]);

return
%_______________________________________________________________________

