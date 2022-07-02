function cg_create_template
% Creation of customized template and prior images for vbm analysis
%
% Run this script by selecting all raw (not normalized) images of your sample
% to be analyzed. This function creates a new template and prior images
% which are smoothed with fwhm of 8 mm.
%
% You can run the whole procedure in one step or only calculate the new
% template and the priors from already segmented data. The second option
% is useful if you calculated segmented images using the first option and
% you don't want to include all subjects to caluclat the template and priors.
%
% output images:
%---------------
%	T1.img		smoothed template
%	gray.img:	smoothed gray matter prior image
%	white.img:	smoothed white matter prior image
%	csf.img:	smoothed CSF prior image
%
% Additionally all temporary images (warped and segmented images) are saved
% and can be removed after calculation of template and prior images:
%
%	w*.img		warped raw image
%	w*_seg1.img	gray matter image
%	w*_seg2.img	white matter image
%	w*_seg3.img	CSF image
%_______________________________________________________________________
% @(#)cg_create_template.m	1.22 Christian Gaser 2006/07/06
% based on code snippets by John Ashburner

global defaults

dseg				= defaults.segment;
dseg.write.wrt_cor		= 0;		% don't write bias corrected image
dseg.estimate.affreg.weight	= fullfile(spm('Dir'),'apriori','brainmask.mnc');	% always use brainmask
dseg.write.cleanup   		= 1;		% clean up partitions

dnrm				= defaults.normalise;
dnrm.write.vox			= [1 1 1];
dnrm.write.bb(1,3)		= -70;		% enlarge z-axis to include more cerebellum
dnrm.estimate.graphics		= 0;		% Crashes otherwise

estimate = spm_input('What to do?',1,'m',...
	'segment images and create template and prior images|use previously segmented images to create template and priors',[1 0], 1);

VG0  = spm_vol(fullfile(spm('Dir'),'templates','T1.mnc')); 

% create graphics window
fg=spm_figure('FindWin','Graphics');
if isempty(fg), fg=spm_figure('Create','Graphics'); end

% estimate segmentation
if estimate
	% change default parameters for normalization
	dnrm.estimate.cutoff	= get_cutoff(dnrm.estimate.cutoff);
	if dnrm.estimate.cutoff ~= Inf
		dnrm.estimate.reg    = get_reg(dnrm.estimate.reg);
		dnrm.estimate.nits   = get_nits(dnrm.estimate.nits);
	end

    dseg.write.mrf = spm_input(['HMRF weighting?'],'+1','m',...
        ['large HMRF weighting (0.45)|medium HMRF weighting (0.3)|small HMRF weighting (0.15)|no HMRF (0)'],...
        [0.45 0.3 0.15 0], 2);

	% get raw images
	V   = spm_vol(spm_get(Inf,'*.IMAGE','Select images'));
	VG1 = spm_vol(deblank(dseg.estimate.priors(1,:)));
	
else	% average images to new template and prior images

	% first select normalized anatomical images
	V   = spm_vol(spm_get(Inf,'w*.IMAGE','Select normalized anatomical images'));
	n = length(V);

	% ... and select same # of segmented images
	V1   = spm_vol(spm_get(n,'w*seg1.IMAGE','Select gray matter images'));
	V2   = spm_vol(spm_get(n,'w*seg2.IMAGE','Select white matter images'));
	V3   = spm_vol(spm_get(n,'w*seg3.IMAGE','Select CSF images'));
end

% Initialise the totals to zero
d = [prod(VG0.dim(1:2)) VG0.dim(3)];
n = zeros(d);
g = zeros(d);
w = zeros(d);
c = zeros(d);
s = zeros(d)+eps;

spm_progress_bar('Init',length(V),'Writing templates','subjects completed');
for i=1:length(V)
  if estimate
	[pth,nam,ext] = fileparts(V(i).fname);

	% Estimate spatial normalisation parameters
	% by matching GM with a GM template
	dseg0 = dseg;
	dseg0.write.mrf = 0;   % don't use MRF in raw space because of huge memory demands
	VT  = spm_segment(V(i),VG0,dseg0);
	VT  = VT(1);
	prm = spm_normalise(VG1,VT,'','','',dnrm.estimate);
	clear VT

	% Create a spatially normalised version of the original image
	VN  = spm_write_sn(V(i),prm,dnrm.write);
	VN.fname = fullfile(pth,['w' nam ext]);
	spm_write_vol(VN,VN.dat);

	% Segment the spatially normalised image
	VT  = spm_segment(VN,eye(4),dseg);
	
	% filenames for segmented images
	if dseg.write.mrf
	   nam2 = [nam '_HMRF'];
	else
	   nam2 = nam;
	end

	VT(1).fname = fullfile(pth,['w' nam2 '_seg1' ext]);
	VT(2).fname = fullfile(pth,['w' nam2 '_seg2' ext]);
	VT(3).fname = fullfile(pth,['w' nam2 '_seg3' ext]);

	% CSF is better estimated using gray and white matter values
	t1 = ((double(VT(1).dat) + double(VT(2).dat) + double(VT(3).dat)) > 0.5*255);
	t2 = (1 - (double(VT(1).dat) + double(VT(2).dat))/255);
	VT(3).dat = uint8(round((t1.*t2)*255));

	% and finally save segmented images
	for j=1:3, spm_write_vol(VT(j),double(VT(j).dat)/255); end

  else
  	VN = V(i);
  	VT = [V1(i) V2(i) V3(i)];
  end
  
  % Subsample to a lower resolution, and add the
  % current images to the totals.
  for p=1:VG0.dim(3),
    M   = VN.mat\VG0.mat*spm_matrix([0 0 p]);
    tn  = spm_slice_vol(VN   ,M,VG0.dim(1:2),1);
    tg  = spm_slice_vol(VT(1),M,VG0.dim(1:2),1);
    tw  = spm_slice_vol(VT(2),M,VG0.dim(1:2),1);
    tc  = spm_slice_vol(VT(3),M,VG0.dim(1:2),1);
    
    msk = find(finite(tn));
    tmp = n(:,p); tmp(msk) = tmp(msk) + tn(msk); n(:,p) = tmp;
    tmp = g(:,p); tmp(msk) = tmp(msk) + tg(msk); g(:,p) = tmp;
    tmp = w(:,p); tmp(msk) = tmp(msk) + tw(msk); w(:,p) = tmp;
    tmp = c(:,p); tmp(msk) = tmp(msk) + tc(msk); c(:,p) = tmp;
    tmp = s(:,p); tmp(msk) = tmp(msk) + 1;       s(:,p) = tmp;
  end;
  
  spm_progress_bar('Set',i);  
end;
spm_progress_bar('Clear');

% Write out the averages (dividing the totals by the number
% of finite observations).
VO = struct('fname','',...
            'dim',[VG0.dim(1:3) spm_type('uint8')],...
            'mat',VG0.mat,...
            'descrip','');

n=reshape(n./s,VG0.dim(1:3));
g=reshape(g./s,VG0.dim(1:3));
w=reshape(w./s,VG0.dim(1:3));
c=reshape(c./s,VG0.dim(1:3));

VO.fname =    'T1.img' ; VO.descrip = ['Customised Template conv 8mm n=' num2str(length(V))]; VO = smooth_vol(VO,n,8); spm_write_vol(VO,VO.dat);
VO.fname =  'gray.img' ; VO.descrip =   ['Grey matter prior conv 8mm n=' num2str(length(V))]; VO = smooth_vol(VO,g,8); spm_write_vol(VO,VO.dat);
VO.fname = 'white.img' ; VO.descrip =  ['White matter prior conv 8mm n=' num2str(length(V))]; VO = smooth_vol(VO,w,8); spm_write_vol(VO,VO.dat);
VO.fname =   'csf.img' ; VO.descrip =           ['CSF prior conv 8mm n=' num2str(length(V))]; VO = smooth_vol(VO,c,8); spm_write_vol(VO,VO.dat);

if spm_input('Delete temporary files?',1,'yes|no',[1,0],1)
    if estimate
        for i=1:length(V)
        	[pth,nam,ext] = fileparts(V(i).fname);
        	fname = fullfile(pth,['w' nam ext]);
            cleanup(fname);
        	fname = fullfile(pth,['w' nam '_seg1' ext]);
            cleanup(fname);
        	fname = fullfile(pth,['w' nam '_seg2' ext]);
            cleanup(fname);
        	fname = fullfile(pth,['w' nam '_seg3' ext]);
            cleanup(fname);
        end
    else
        cleanup(str2mat(V.fname));
        cleanup(str2mat(V1.fname));
        cleanup(str2mat(V2.fname));
        cleanup(str2mat(V3.fname));
    end
end

return

%-----------------------------------------------------------------------
function cleanup(P)
% Delete temporary files.
%
for i=1:size(P,1)
    [fpath,fname,ext] = fileparts(deblank(P(i,:)));
    spm_unlink(fullfile(fpath,[fname '.hdr']));
    spm_unlink(fullfile(fpath,[fname '.img']));
end

return;

%-----------------------------------------------------------------------

function V = smooth_vol(V,Y,s)
% smooth volume with gaussian kernel
% based on spm_smooth

% compute parameters for spm_conv_vol
%-----------------------------------------------------------------------
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
V.dat = Y;
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
% FORMAT cutoff = get_cutoof(cutoff)

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

