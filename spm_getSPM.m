function [SPM,xSPM] = spm_getSPM
% computes a specified and thresholded SPM/PPM following parameter estimation
% FORMAT [SPM,xSPM] = spm_getSPM;
%
% xSPM      - structure containing SPM, distribution & filtering details
% .swd      - SPM working directory - directory containing current SPM.mat
% .title    - title for comparison (string)
% .Z        - minimum of statistic images {filtered on u and k}
% .n        - conjunction number <= number of contrasts        
% .STAT     - distribution {Z, T, X, F or P}     
% .df       - degrees of freedom [df{interest}, df{residual}]
% .STATstr  - description string     
% .Ic       - indices of contrasts (in SPM.xCon)
% .Im       - indices of masking contrasts (in xCon)
% .pm       - p-value for masking (uncorrected)
% .Ex       - flag for exclusive or inclusive masking
% .u        - height threshold
% .k        - extent threshold {voxels}
% .XYZ      - location of voxels {voxel coords}
% .XYZmm    - location of voxels {mm}
% .S        - search Volume {voxels}
% .R        - search Volume {resels}
% .FWHM     - smoothness {voxels}     
% .M        - voxels -> mm matrix
% .iM       - mm -> voxels matrix
% .VOX      - voxel dimensions {mm} - column vector
% .DIM      - image dimensions {voxels} - column vector
% .Vspm     - Mapped statistic image(s)
% .Ps       - list of P values for voxels at SPM.xVol.XYZ (used by FDR)
%
% Required fields of SPM
%
% xVol   - structure containing details of volume analysed
%
% xX     - Design Matrix structure
%        - (see spm_spm.m for structure)
%
% xCon   - Contrast definitions structure array
%        - (see also spm_FcUtil.m for structure, rules & handling)
% .name  - Contrast name
% .STAT  - Statistic indicator character ('T', 'F' or 'P')
% .c     - Contrast weights (column vector contrasts)
% .X0    - Reduced design matrix data (spans design space under Ho)
%          Stored as coordinates in the orthogonal basis of xX.X from spm_sp
%          (Matrix in SPM99b)  Extract using X0 = spm_FcUtil('X0',...
% .iX0   - Indicates how contrast was specified:
%          If by columns for reduced design matrix then iX0 contains the
%          column indices. Otherwise, it's a string containing the
%          spm_FcUtil 'Set' action: Usually one of {'c','c+','X0'}
% .X1o   - Remaining design space data (X1o is orthogonal to X0)
%          Stored as coordinates in the orthogonal basis of xX.X from spm_sp
%          (Matrix in SPM99b)  Extract using X1o = spm_FcUtil('X1o',...
% .eidf  - Effective interest degrees of freedom (numerator df)
%        - Or effect-size threshold for Posterior probability
% .Vcon  - Name of contrast (for 'T's) or ESS (for 'F's) image
% .Vspm  - Name of SPM image
%
% In addition, the xCon.mat file is updated. For newly evaluated
% contrasts, SPM images (spmT_????.{img,hdr}) are written, along with
% contrast (con_????.{img,hdr}) images for SPM{T}'s, or Extra
% Sum-of-Squares images (ess_????.{img,hdr}) for SPM{F}'s.
% 
% The contrast images are the weighted sum of the parameter images,
% where the weights are the contrast weights, and are uniquely
% estimable since contrasts are checked for estimability by the
% contrast manager. These contrast images (for appropriate contrasts)
% are suitable summary images of an effect at this level, and can be
% used as input at a higher level when effecting a random effects
% analysis. (Note that the ess_????.{img,hdr} and
% SPM{T,F}_????.{img,hdr} images are not suitable input for a higher
% level analysis.) See spm_RandFX.man for further details.
%
%_______________________________________________________________________
%
% spm_getSPM prompts for an SPM and applies thresholds {u & k}
% to a point list of voxel values (specified with their locations {XYZ})
% This allows the SPM be displayed and characterized in terms of regionally 
% significant effects by subsequent routines.
% 
% For general linear model Y = XB + E with data Y, design matrix X,
% parameter vector B, and (independent) errors E, a contrast c'B of the
% parameters (with contrast weights c) is estimated by c'b, where b are
% the parameter estimates given by b=pinv(X)*Y.
% 
% Either single contrasts can be examined or conjunctions of different
% contrasts. Contrasts are estimable linear combinations of the
% parameters, and are specified using the SPM contrast manager
% interface [spm_conman.m]. SPMs are generated for the null hypotheses
% that the contrast is zero (or zero vector in the case of
% F-contrasts). See the help for the contrast manager [spm_conman.m]
% for a further details on contrasts and contrast specification.
% 
% A conjunction assesses the conjoint expression of multiple effects. The
% conjunction SPM is the minimum of the component SPMs defined by the
% multiple contrasts.  Inference on the minimum statistics can be
% performed in different ways.  Inference on the Conjunction Null (one or
% more of the effects null) is accomplished by assessing the minimum as
% if it were a single statistic; one rejects the conjunction null in
% favor of the alternative that k=nc, that the number of active effects k
% is equal to the number of contrasts nc.  No assumptions are needed on
% the dependence between the tests.
%
% Another approach is to make inference on the Global Null (all effects
% null).  Rejecting the Global Null of no (u=0) effects real implies an
% alternative that k>0, that one or more effects are real.   A third
% Intermediate approach, is to use a null hypothesis of no more than u
% effects are real.  Rejecting the intermediate null that k<=u implies an
% alternative that k>u, that more than u of the effects are real.  
% 
% The Global and Intermediate nulls use results for minimum fields which
% require the SPMs to be identically distributed and independent. Thus,
% all component SPMs must be either SPM{t}'s, or SPM{F}'s with the same
% degrees of freedom. Independence is roughly guaranteed for large
% degrees of freedom (and independent data) by ensuring that the
% contrasts are "orthogonal". Note that it is *not* the contrast weight
% vectors per se that are required to be orthogonal, but the subspaces of
% the data space implied by the null hypotheses defined by the contrasts
% (c'pinv(X)). Furthermore, this assumes that the errors are
% i.i.d. (i.e. the estimates are maximum likelihood or Gauss-Markov. This
% is the default in spm_spm).  
%
% To ensure approximate independence of the component SPMs in the case of
% the global or intermediate null, non-orthogonal contrasts are serially
% orthogonalised in the order specified, possibly generating new
% contrasts, such that the second is orthogonal to the first, the third
% to the first two, and so on.  Note that significant inference on the
% global null only allows one to conclude that one or more of the effects
% are real.  Significant inference on the conjunction null allows one to
% conclude that all of the effects are real.
%
% Masking simply eliminates voxels from the current contrast if they
% do not survive an uncorrected p value (based on height) in one or
% more further contrasts.  No account is taken of this masking in the
% statistical inference pertaining to the masked contrast.
% 
% The SPM is subject to thresholding on the basis of height (u) and the
% number of voxels comprising its clusters {k}. The height threshold is
% specified as above in terms of an [un]corrected p value or
% statistic.  Clusters can also be thresholded on the basis of their
% spatial extent. If you want to see all voxels simply enter 0.  In this
% instance the 'set-level' inference can be considered an 'omnibus test'
% based on the number of clusters that obtain.
%
% BAYESIAN INFERENCE AND PPMS - POSTERIOR PROBABILITY MAPS
%
% If conditional estimates are available (and your contrast is a T
% contrast) then you are asked whether the inference should be 'Bayesian'
% or 'classical' (using GRF).  If you choose Bayesian the contrasts are of
% conditional (i.e. MAP) estimators and the inference image is a
% posterior probability map (PPM).  PPMs encode the probability that the
% contrast exceeds a specified threshold.  This threshold is stored in
% the xCon.eidf.  Subsequent plotting and tables will use the conditional
% estimates and associated posterior or conditional probabilities.
% 
% see spm_results_ui.m for further details of the SPM results section.
% see also spm_contrasts.m
%_______________________________________________________________________
% @(#)spm_getSPM.m	2.51 Andrew Holmes, Karl Friston & Jean-Baptiste Poline 03/05/22
% $Id: spm_getSPM.m,v 1.15 2005/01/21 14:33:20 nichols Exp $ UM Biostat


SCCSid = '2.51';

%-GUI setup
%-----------------------------------------------------------------------
SPMid  = spm('SFnBanner',mfilename,SCCSid);
spm_help('!ContextHelp',mfilename)

%-Select SPM.mat & note SPM results directory
%-----------------------------------------------------------------------
swd    = spm_str_manip(spm_get(1,'SPM.mat','Select SPM.mat'),'H');

%-Preliminaries...
%=======================================================================

%-Load SPM.mat
%-----------------------------------------------------------------------
load(fullfile(swd,'SPM.mat'));
SPM.swd = swd;

%-Get volumetric data from SPM.mat
%-----------------------------------------------------------------------
try
	xX   = SPM.xX;				%-Design definition structure
	XYZ  = SPM.xVol.XYZ;			%-XYZ coordinates
	S    = SPM.xVol.S;			%-search Volume {voxels}
	R    = SPM.xVol.R;			%-search Volume {resels}
	M    = SPM.xVol.M(1:3,1:3);		%-voxels to mm matrix
	VOX  = sqrt(diag(M'*M))';		%-voxel dimensions
catch
    
	% check the model has been estimated
	%---------------------------------------------------------------
	str = {	'This model has not been estimated.';...
		'Would you like to estimate it now?'};
	if spm_input(str,1,'bd','yes|no',[1,0],1)
		[SPM] = spm_spm(SPM);
	else
		return
	end
end


%-Contrast definitions
%=======================================================================

%-Load contrast definitions (if available)
%-----------------------------------------------------------------------
try
	xCon = SPM.xCon;
catch
	xCon = {};
end


%=======================================================================
% - C O N T R A S T S ,  S P M   C O M P U T A T I O N ,   M A S K I N G
%=======================================================================

%-Get contrasts
%-----------------------------------------------------------------------
[Ic,xCon] = spm_conman(xX,xCon,'T&F',Inf,...
		'	Select contrasts...',' for conjunction',1);


nc       = length(Ic);  % Number of contrasts

%-Allow user to extend the null hypothesis for conjunctions
%
% n: conjunction number
% u: Null hyp is k<=u effects real; Alt hyp is k>u effects real 
%    (NB Here u is from Friston et al 2004 paper, not statistic thresh).
%                  u         n      
% Conjunction Null nc-1      1     |    u = nc-n
% Intermediate     1..nc-2   nc-u  |    #effects under null <= u
% Global Null      0         nc    |    #effects under alt  > u,  >= u+1
%----------------------------------+-------------------------------------
if (nc > 1)
  if nc==2
    But='Conjunction|Global';      Val=[1 nc]; 
  else
    But='Conj''n|Intermed|Global'; Val=[1 NaN nc]; 
  end
  n = spm_input('Null hyp. to assess?','+1','b',But,Val,1);
  if isnan(n)
    if nc==3,
      n = nc-1;
    else
      n = nc-spm_input('Effects under null ','0','n1','1',nc-1);
    end
  end
else
  n = 1;
end


%-Enforce orthogonality of multiple contrasts for conjunction
% (Orthogonality within subspace spanned by contrasts)
%-----------------------------------------------------------------------
if nc>1 & n>1 & ~spm_FcUtil('|_?',xCon(Ic), xX.xKXs)

    OrthWarn = 0;

    %-Successively orthogonalise
    %-NB: This loop is peculiarly controlled to account for the
    %     possibility that Ic may shrink if some contrasts diasppear
    %     on orthogonalisation (i.e. if there are colinearities)
    %-------------------------------------------------------------------
    i = 1; while(i < nc), i = i + 1;
    
	%-Orthogonalise (subspace spanned by) contrast i wirit previous
	%---------------------------------------------------------------
	oxCon = spm_FcUtil('|_',xCon(Ic(i)), xX.xKXs, xCon(Ic(1:i-1)));
    
	%-See if this orthogonalised contrast has already been entered
	% or is colinear with a previous one. Define a new contrast if
	% neither is the case.
	%---------------------------------------------------------------
	d     = spm_FcUtil('In',oxCon,xX.xKXs,xCon);

	if spm_FcUtil('0|[]',oxCon,xX.xKXs)

	    %-Contrast was colinear with a previous one - drop it
	    %-----------------------------------------------------------
	    Ic(i) = [];
	    i     = i - 1;

	elseif any(d)

	    %-Contrast unchanged or already defined - note index
	    %-----------------------------------------------------------
	    Ic(i) = min(d);

	else

	    OrthWarn = OrthWarn + 1;

	    %-Define orthogonalised contrast as new contrast
	    %-----------------------------------------------------------
	    oxCon.name = [xCon(Ic(i)).name,' (orth. w.r.t {',...
    	    	sprintf('%d,',Ic(1:i-2)), sprintf('%d})',Ic(i-1))];
    	    xCon  = [xCon, oxCon];
	    Ic(i) = length(xCon); 
	end

    end % while...
    
    if OrthWarn
      warning(sprintf('Contrasts changed!  %d contrasts orthogonalized ',...
		      'to allow conjunction inf.', OrthWarn))
    end
end % if (orthogonalize)


%-Get contrasts for masking
%-----------------------------------------------------------------------
if spm_input('mask with other contrast(s)','+1','y/n',[1,0],2)
	[Im,xCon] = spm_conman(xX,xCon,'T&F',-Inf,...
		'Select contrasts for masking...',' for masking',1);

	%-Threshold for mask (uncorrected p-value)
	%---------------------------------------------------------------
	pm = spm_input('uncorrected mask p-value','+1','r',0.05,1,[0,1]);

	%-Inclusive or exclusive masking
	%---------------------------------------------------------------
	Ex = spm_input('nature of mask','+1','b','inclusive|exclusive',[0,1]);
else
	Im = [];
	pm = [];
	Ex = [];
end


%-Create/Get title string for comparison
%-----------------------------------------------------------------------
if nc == 1
	str  = xCon(Ic).name;
else
	str  = [sprintf('contrasts {%d',Ic(1)),sprintf(',%d',Ic(2:end)),'}'];
	if n==nc
	    str = [str ' (global null)'];
	elseif n==1
	    str = [str ' (conj. null)'];
	else
	    str = [str sprintf(' (Ha: k>=%d)',(nc-n)+1)];
	end
end
if Ex
	mstr = 'masked [excl.] by';
else
	mstr = 'masked [incl.] by';
end
if length(Im) == 1
	str  = sprintf('%s (%s %s at p=%g)',str,mstr,xCon(Im).name,pm);

elseif ~isempty(Im)
	str  = [sprintf('%s (%s {%d',str,mstr,Im(1)),...
		sprintf(',%d',Im(2:end)),...
		sprintf('} at p=%g)',pm)];
end
titlestr     = spm_input('title for comparison','+1','s',str);



%-Bayesian or classical Inference?
%-----------------------------------------------------------------------
if isfield(SPM,'PPM') & xCon(Ic(1)).STAT == 'T'

    if nc == 1 & isempty(xCon(Ic).Vcon)

	if spm_input('Inference',1,'b',{'Bayesian','classical'},[1 0]);

	% set STAT to 'P'
        %---------------------------------------------------------------
	xCon(Ic).STAT = 'P';

	%-Get Bayesian threshold (Gamma) stored in xCon(Ic).eidf
	% The default is one conditional s.d. of the contrast
        %---------------------------------------------------------------
	str           = 'threshold {default: prior s.d.}';
	Gamma         = sqrt(xCon(Ic).c'*SPM.PPM.Cb*xCon(Ic).c);
	xCon(Ic).eidf = spm_input(str,'+1','e',sprintf('%0.2f',Gamma));

	end
    end
end


%-Compute & store contrast parameters, contrast/ESS images, & SPM images
%=======================================================================
SPM.xCon = xCon;
SPM      = spm_contrasts(SPM,unique([Ic,Im]));
xCon     = SPM.xCon;
VspmSv   = cat(1,xCon(Ic).Vspm);
STAT     = xCon(Ic(1)).STAT;

%-Check conjunctions - Must be same STAT w/ same df
%-----------------------------------------------------------------------
if (nc > 1) & (any(diff(double(cat(1,xCon(Ic).STAT)))) | ...
	      any(abs(diff(cat(1,xCon(Ic).eidf))) > 1))
	error('illegal conjunction: can only conjoin SPMs of same STAT & df')
end


%-Degrees of Freedom and STAT string describing marginal distribution
%-----------------------------------------------------------------------
df          = [xCon(Ic(1)).eidf xX.erdf];
if nc>1
	if n>1
	    str = sprintf('^{%d \\{Ha:k\\geq%d\\}}',nc,(nc-n)+1);
	else
	    str = sprintf('^{%d \\{Ha:k=%d\\}}',nc,(nc-n)+1);
	end
else
	str = '';
end

switch STAT
case 'T'
	STATstr = sprintf('%c%s_{%.0f}','T',str,df(2));
case 'F'
	STATstr = sprintf('%c%s_{%.0f,%.0f}','F',str,df(1),df(2));
case 'P'
	STATstr = sprintf('%s^{%0.2f}','PPM',df(1));
end


%-Compute (unfiltered) SPM pointlist for masked conjunction requested
%=======================================================================
fprintf('\t%-32s: %30s\n','SPM computation','...initialising')         %-#


%-Compute conjunction as minimum of SPMs
%-----------------------------------------------------------------------
Z         = Inf;
for i     = Ic
	Z = min(Z,spm_get_data(xCon(i).Vspm,XYZ));
end

% P values for False Discovery FDR rate computation (all search voxels)
%=======================================================================
switch STAT
case 'T'
	Ps = (1 - spm_Tcdf(Z,df(2))).^n;
case 'P'
	Ps = (1 - Z).^n;
case 'F'
	Ps = (1 - spm_Fcdf(Z,df)).^n;
end


%-Compute mask and eliminate masked voxels
%-----------------------------------------------------------------------
for i = Im
	fprintf('%s%30s',sprintf('\b')*ones(1,30),'...masking')

	Mask = spm_get_data(xCon(i).Vspm,XYZ);
	um   = spm_u(pm,[xCon(i).eidf,xX.erdf],xCon(i).STAT);
	if Ex
		Q = Mask <= um;
	else
		Q = Mask >  um;
	end
	XYZ       = XYZ(:,Q);
	Z         = Z(Q);
	if isempty(Q)
		fprintf('\n')                                           %-#
    		warning(sprintf('No voxels survive masking at p=%4.2f',pm))
		break
	end
end

%-clean up interface
%-----------------------------------------------------------------------
fprintf('\t%-32s: %30s\n','SPM computation','...done')         %-#
spm('Pointer','Arrow')



%=======================================================================
% - H E I G H T   &   E X T E N T   T H R E S H O L D S
%=======================================================================

%-Height threshold - classical inference
%-----------------------------------------------------------------------
u      = -Inf;
k      = 0;
if STAT ~= 'P'


    %-Get height threshold
    %-------------------------------------------------------------------
    str = 'FWE|FDR|none';
    % str = 'FWE|none';      % Use this line to disable FDR threshold
    switch spm_input('p value adjustment to control','+1','b',str,[],1)


	case 'FWE' % family-wise false positive rate
        %---------------------------------------------------------------
	u  = spm_input('p value (family-wise error)','+0','r',0.05,1,[0,1]);
	u  = spm_uc(u,df,STAT,R,n,S);

	case 'FDR' % False discovery rate
	%---------------------------------------------------------------	
	u  = spm_input('p value (false discovery rate)','+0','r',0.05,1,[0,1]);
	u  = spm_uc_FDR(u,df,STAT,n,VspmSv,0);

	otherwise  %-NB: no adjustment
	% p for conjunctions is p of the conjunction SPM
   	%---------------------------------------------------------------
	u  = spm_input(['threshold {',STAT,' or p value}'],'+0','r',0.001,1);
	if u <= 1; u = spm_u(u^(1/n),df,STAT); end

    end

%-Height threshold - Bayesian inference
%-----------------------------------------------------------------------
elseif STAT == 'P'

	u  = spm_input(['p value threshold for PPM'],'+0','r',.95,1);

end % (if STAT)

%-Calculate height threshold filtering
%-------------------------------------------------------------------
Q      = find(Z > u);

%-Apply height threshold
%-------------------------------------------------------------------
Z      = Z(:,Q);
XYZ    = XYZ(:,Q);
if isempty(Q)
	warning(sprintf('No voxels survive height threshold u=%0.2g',u))
end


%-Extent threshold (disallowed for conjunctions)
%-----------------------------------------------------------------------
noniso = 0;
if ~isempty(XYZ) & nc == 1

	pk     = spm_input('extent threshold {k or p-value}','+1','r',0,1,[0,Inf]);
	if (pk < 1) & (pk > 0)
    	extent_FWE = spm_input('p value (extent)','+1','b','uncorrected|FWE corrected',[0 1],1);
	end

    if (STAT ~= 'P') & (STAT == 'F')
        tmpdf  = df;                    %-df (F-stat)
        tmpdfw = [tmpdf(2) tmpdf(2)];   %-df for RPV field
    elseif (STAT ~= 'P') & (STAT == 'T')
        tmpdf  = df(2); %-df (t-stat)
        tmpdfw = tmpdf; %-df for RPV field
    end

	if isfield(SPM.xVol,'VRpv')
		if exist(SPM.xVol.VRpv.fname)
	    	noniso = spm_input('Correct for non-stationarity smoothness?','+1','b','no|yes',[0 1],2);
	    end
    end

    %-if stationary, then RPV image has infinite df
    if ~noniso
        tmpdfw  = Inf*tmpdfw;
    end

	FWHM = SPM.xVol.FWHM;
	v2r   = 1/prod(FWHM(~isinf(FWHM)));
	if (pk < 1) & (pk > 0)
		if extent_FWE
			Pk = 1;
			k = 0;
			while (Pk >= pk & k<S)
				k = k + 1;
				[tmp1 Pk tmp2 Pn] = stat_thresh(R,S,1,[tmpdf;tmpdfw],0.05,u,k*v2r);
			end
		else
			Pn = 1;
			k = 0;
			while (Pn >= pk & k<S)
				k = k + 1;
				[tmp1 Pk tmp2 Pn] = stat_thresh(R,S,1,[tmpdf;tmpdfw],0.05,u,k*v2r);
			end
		end
	else
		k = pk;
	end

    %-Calculate extent threshold filtering
    %-------------------------------------------------------------------
    if noniso
		fprintf('Use local RPV values to correct for non-stationarity of smoothness.\n');
	    A     = spm_clusters(XYZ);
		[N Z2 XYZ2 A2 V2R] = spm_max_nS(Z,XYZ,SPM.xVol.VRpv);
		N = V2R/v2r;
    	Q     = [];
    	for i = 1:max(A)
        	j = find(A == i);
        	if N(min(find(A2 == i))) >= k; Q = [Q j]; end
    	end
    else
	    A     = spm_clusters(XYZ);
    	Q     = [];
    	for i = 1:max(A)
        	j = find(A == i);
        	if length(j) >= k; Q = [Q j]; end
    	end
	end

    % ...eliminate voxels
    %-------------------------------------------------------------------
    Z     = Z(:,Q);
    XYZ   = XYZ(:,Q);
    if isempty(Q)
	warning(sprintf('No voxels survive extent threshold k=%0.2g',k))
    end

else

    k = 0;

end % (if ~isempty(XYZ))


%=======================================================================
% - E N D
%=======================================================================
fprintf('\t%-32s: %30s\n','SPM computation','...done')         %-#

%-Assemble output structures of unfiltered data
%=======================================================================
xSPM   = struct('swd',		swd,...
		'title',	titlestr,...
		'Z',		Z,...
		'n',		n,...
		'STAT',		STAT,...
		'df',		df,...
		'STATstr',	STATstr,...
		'Ic',		Ic,...
		'Im',		Im,...
		'pm',		pm,...
		'Ex',		Ex,...
		'u',		u,...
		'k',		k,...
		'XYZ',		XYZ,...
		'XYZmm',	SPM.xVol.M(1:3,:)*[XYZ; ones(1,size(XYZ,2))],...
		'S',		SPM.xVol.S,...
		'R',		SPM.xVol.R,...
		'FWHM',		SPM.xVol.FWHM,...
		'M',		SPM.xVol.M,...
		'iM',		SPM.xVol.iM,...
		'DIM',		SPM.xVol.DIM,...
		'VOX',		VOX,...
		'Vspm',		VspmSv,...
		'Ps',		Ps);

% RESELS per voxel (density) if it exists
%-----------------------------------------------------------------------
if noniso
	xSPM.VRpv = SPM.xVol.VRpv;
else
	try, xSPM.VRpv  = SPM.VRpv;       end
end
