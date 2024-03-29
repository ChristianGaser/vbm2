function varargout = cg_spm_ui(varargin)
% Setting up the general linear model for independent data
% FORMATs (given in Programmers Help)
%_______________________________________________________________________
%
% cg_spm_ui.m configures the design matrix (describing the general
% linear model), data specification, and other parameters necessary for
% the statistical analysis. These parameters are saved in a
% configuration file (SPMcfg.mat) in the current directory, and are
% passed on to spm_spm.m which estimates the design. Inference on these
% estimated parameters is then handled by the SPM results section.
%
% A separate program (spm_spm_fmri_ui.m) handles design configuration
% for fMRI time series, though this program can be used for fMRI data
% when observations can be regarded as independent.
%
% ----------------------------------------------------------------------
%
% Various data and parameters need to be supplied to specify the design:
%       * the image files
%       * indicators of the corresponding condition/subject/group
%       * any covariates, nuisance variables, or design matrix partitions
%       * the type of global normalisation (if any)
%       * grand mean scaling options
%       * thresholds and masks defining the image volume to analyse
%
% The interface supports a comprehensive range of options for all these
% parameters, which are described below in the order in which the
% information is requested. Rather than ask for all these parameters,
% cg_spm_ui.m uses a "Design Definition", a structure describing the
% options and defaults appropriate for a particular analysis. Thus,
% once the user has chosen a design, a subset of the following prompts
% will be presented.
%
% If the pre-specified design definitions don't quite
% have the combination of options you want, you can pass a custom design
% structure D to be used as parameter: cg_spm_ui('cfg',D). The format of
% the design structure and option definitions are given in the programmers
% help, at the top of the main body of the code.
%
%                           ----------------
%
% Design class & Design type
% ==========================
%
% Unless a design definition is passed to cg_spm_ui.m as a parameter,
% the user is prompted first to select a design class, and then to
% select a design type from that class.
%
% The designs are split into three classes:
%   i) Basic stats: basic models for simple statistics
%      These specify designs suitable for simple voxel-by-voxel analyses.
%       - one-sample t-test
%       - two-sample t-test
%       - paired t-test
%       - one way Anova
%       - one way Anova (with constant)
%       - one way Anova (within subject)
%       - simple regression (equivalent to correlation)
%       - multiple regression
%       - multiple regression  (with constant)
%       - basic AnCova (ANalysis of COVAriance)
%         (essentially a two-sample t-test with a nuisance covariate)
%
%  ii) PET models: models suitable for analysis of PET/SPECT experiments
%       - Single-subject: conditions & covariates
%       - Single-subject: covariates only
%
%       - Multi-subj: conditions & covariates
%       - Multi-subj: cond x subj  interaction & covariates
%       - Multi-subj: covariates only
%       - Multi-group: conditions & covariates
%       - Multi-group: covariates only
%
%       - Population main effect: 2 cond's, 1 scan/cond (paired t-test)
%       - Dodgy population main effect: >2 cond's, 1 scan/cond
%       - Compare-populations: 1 scan/subject (two sample t-test)
%       - Compare-populations: 1 scan/subject (AnCova)
%
%       - The Full Monty... (asks you everything!)
%
% iii) SPM96 PET models: models used in SPM96 for PET/SPECT
%      These models are provided for backward compatibility, but as they
%      don't include some of the advanced modelling features, we recommend
%      you switch to the new (SPM99) models at the earliest opportunity.
%       - SPM96:Single-subject: replicated conditions
%       - SPM96:Single-subject: replicated conditions & covariates
%       - SPM96:Single-subject: covariates only
%       - SPM96:Multi-subject: different conditions
%       - SPM96:Multi-subject: replicated conditions
%       - SPM96:Multi-subject: different conditions & covariates
%       - SPM96:Multi-subject: replicated conditions & covariates
%       - SPM96:Multi-subject: covariates only
%       - SPM96:Multi-group: different conditions
%       - SPM96:Multi-group: replicated conditions
%       - SPM96:Multi-group: different conditions & covariates
%       - SPM96:Multi-group: replicated conditions & covariates
%       - SPM96:Multi-group: covariates only
%       - SPM96:Compare-groups: 1 scan per subject
%
%
% Random effects, generalisability, population inference...
% =========================================================
%
% Note that SPM only considers a single component of variance, the
% residual error variance. When there are repeated measures, all
% analyses with SPM are fixed effects analyses, and inference only
% extends to the particular subjects under consideration (at the times
% they were imaged).
%
% In particular, the multi-subject and multi-group designs ignore the
% variability in response from subject to subject. Since the
% scan-to-scan (within-condition, within-subject variability is much
% smaller than the between subject variance which is ignored), this can
% lead to detection of group effects that are not representative of the
% population(s) from which the subjects are drawn. This is particularly
% serious for multi-group designs comparing two groups. If inference
% regarding the population is required, a random effects analysis is
% required.
%
% However, random effects analyses can be effected by appropriately
% summarising the data, thereby collapsing the model such that the
% residual variance for the new model contains precisely the variance
% components needed for a random effects analysis. In many cases, the
% fixed effects models here can be used as the first stage in such a
% two-stage procedure to produce appropriate summary data, which can
% then be used as raw data for a second-level analysis. For instance,
% the "Multi-subj: cond x subj  interaction & covariates" design can be
% used to write out an image of the activation for each subject. A
% simple t-test on these activation images then turns out to be
% equivalent to a mixed-effects analysis with random subject and
% subject by condition interaction effects, inferring for the
% population based on this sample of subjects (strictly speaking the
% design would have to be balanced, with equal numbers of scans per
% condition per subject, and also only two conditions per subject). For
% additional details, see spm_RandFX.man.
%
%                           ----------------
%
% Selecting image files & indicating conditions
% =============================================
%
% You may now be prompted to specify how many studies, subjects and
% conditions you have, and then will be promted to select the scans.
%
% The data should all have the same orientation and image and voxel size.
%
% File selection is handled by spm_get.m - the help for which describes
% efficient use of the interface.
%
% You may be asked to indicate the conditions for a set of scans, with
% a prompt like "[12] Enter conditions? (2)". For this particular
% example you need to indicate for 12 scans the corresponding
% condition, in this case from 2 conditions. Enter a vector of
% indicators, like '0 1 0 1...', or a string of indicators, like
% '010101010101' or '121212121212', or 'rararararara'. (This
% "conditions" input is handled by spm_input.m, where comprehensive
% help can be found.)
%
%                           ----------------
%
% Covariate & nuisance variable entry
% ===================================
%
% * If applicable, you'll be asked to specify covariates and nuisance
% variables. Unlike SPM94/5/6, where the design was partitioned into
% effects of interest and nuisance effects for the computation of
% adjusted data and the F-statistic (which was used to thresh out
% voxels where there appeared to be no effects of interest), SPM99 does
% not partition the design in this way. The only remaining distinction
% between effects of interest (including covariates) and nuisance
% effects is their location in the design matrix, which we have
% retained for continuity.  Pre-specified design matrix partitions can
% be entered. (The number of covariates / nuisance variables specified,
% is actually the number of times you are prompted for entry, not the
% number of resulting design matrix columns.) You will be given the
% opportunity to name the covariate.
% 
% * Factor by covariate interactions: For covariate vectors, you may be
% offered a choice of interaction options. (This was called "covariate
% specific fits" in SPM95/6.) The full list of possible options is:
%       - <none>
%       - with replication
%       - with condition (across group)
%       - with subject (across group)
%       - with group
%       - with condition (within group)
%       - with subject (within group)
%
% * Covariate centering: At this stage may also be offered "covariate
% centering" options. The default is usually that appropriate for the
% interaction chosen, and ensures that main effects of the interacting
% factor aren't affected by the covariate. You are advised to choose
% the default, unless you have other modelling considerations. The full
% list of possible options is:
%       - around overall mean
%       - around replication means
%       - around condition means (across group)
%       - around subject means (across group)
%       - around group means
%       - around condition means (within group)
%       - around subject means (within group)
%       - <no centering>
%
%                           ----------------
%
% Global options
% ==============
%
% Depending on the design configuration, you may be offered a selection
% of global normalisation and scaling options:
%
% * Method of global flow calculation
%       - SPM96:Compare-groups: 1 scan per subject
%       - None (assumming no other options requiring the global value chosen)
%       - User defined (enter your own vector of global values)
%       - SPM standard: mean voxel value (within per image fullmean/8 mask)
%
% * Grand mean scaling : Scaling of the overall grand mean simply
% scales all the data by a common factor such that the mean of all the
% global values is the value specified. For qualitative data, this puts
% the data into an intuitively accessible scale without altering the
% statistics. When proportional scaling global normalisation is used
% (see below), each image is seperately scaled such that it's global
% value is that specified (in which case the grand mean is also
% implicitly scaled to that value). When using AnCova or no global
% normalisation, with data from different subjects or sessions, an
% intermediate situation may be appropriate, and you may be given the
% option to scale group, session or subject grand means seperately. The
% full list of possible options is:
%       - scaling of overall grand mean
%       - caling of replication grand means
%       - caling of condition grand means (across group)
%       - caling of subject grand means (across group)
%       - caling of group grand means
%       - caling of condition (within group) grand means
%       - caling of subject (within group) grand means
%       - implicit in PropSca global normalisation)
%       - no grand Mean scaling>'
%
% * Global normalisation option : Global nuisance effects are usually
% accounted for either by scaling the images so that they all have the
% same global value (proportional scaling), or by including the global
% covariate as a nuisance effect in the general linear model (AnCova).
% Much has been written on which to use, and when. Basically, since
% proportional scaling also scales the variance term, it is appropriate
% for situations where the global measurement predominantly reflects
% gain or sensitivity. Where variance is constant across the range of
% global values, linear modelling in an AnCova approach has more
% flexibility, since the model is not restricted to a simple
% proportional regression.
%
% Considering AnCova global normalisation, since subjects are unlikely
% to have the same relationship between global and local measurements,
% a subject-specific AnCova ("AnCova by subject"), fitting a different
% slope and intercept for each subject, would be preferred to the
% single common slope of a straight AnCova. (Assumming there's enough
% scans per subject to estimate such an effect.) This is basically an
% interaction of the global covariate with the subject factor. You may
% be offered various AnCova options, corresponding to interactions with
% various factors according to the design definition: The full list of
% possible options is:
%       - AnCova
%       - AnCova by replication
%       - AnCova by condition (across group)
%       - AnCova by subject (across group)
%       - AnCova by group
%       - AnCova by condition (within group)
%       - AnCova by subject (within group)
%       - Proportional scaling
%       - <no global normalisation>
%
% Since differences between subjects may be due to gain and sensitivity
% effects, AnCova by subject could be combined with "grand mean scaling
% by subject" to obtain a combination of between subject proportional
% scaling and within subject AnCova.
%
% * Global centering: Lastly, for some designs using AnCova, you will
% be offered a choice of centering options for the global covariate. As
% with covariate centering, this is only relevant if you have a
% particular interest in the parameter estimates. Usually, the default
% of a centering corresponding to the AnCova used is chosen. The full
% list of possible options is:
%       - around overall mean
%       - around replication means
%       - around condition means (across group)
%       - around subject means (across group)
%       - around group means
%       - around condition means (within group)
%       - around subject means (within group)
%       - <no centering>
%       - around user specified value
%       - (as implied by AnCova)
%       - GM (The grand mean scaled value)
%       - (redundant: not doing AnCova)
%
%
%
% Note that this is a logical ordering for the global options, which is
% not the order used by the interface due to algorithm constraints. The
% interface asks for the options in this order:
%       - Global normalisation
%       - Grand mean scaling options
%         (if not using proportional scaling global normalisation)
%       - Value for grand mean scaling  proportional scaling GloNorm
%         (if appropriate)
%       - Global centering options
%       - Value for global centering (if "user-defined" chosen)
%       - Method of calculation
%
%                           ----------------
%
% Masking options
% ===============
%
% The mask specifies the voxels within the image volume which are to be
% assessed. SPM supports three methods of masking. The volume analysed
% is the intersection of all masks:
%
%   i) Threshold masking : "Analysis threshold"
%       - images are thresholded at a given value and only voxels at
%         which all images exceed the threshold are included in the
%         analysis.
%       - The threshold can be absolute, or a proportion of the global
%         value (after scaling), or "-Inf" for no threshold masking.
%       - (This was called "Grey matter threshold" in SPM94/5/6)
%
%  ii) Implicit masking
%       - An "implicit mask" is a mask implied by a particular voxel
%         value. Voxels with this mask value are excluded from the
%         analysis.
%       - For image data-types with a representation of NaN
%         (see spm_type.m), NaN's is the implicit mask value, (and
%         NaN's are always masked out).
%       - For image data-types without a representation of NaN, zero is
%         the mask value, and the user can choose whether zero voxels
%          should be masked out or not.
%
% iii) Explicit masking
%       - Explicit masks are other images containing (implicit) masks
%         that are to be applied to the current analysis.
%       - All voxels with value NaN (for image data-types with a
%         representation of NaN), or zero (for other data types) are
%         excluded from the analysis.
%       - Explicit mask images can have any orientation and voxel/image
%         size. Nearest neighbour interpolation of a mask image is used if
%         the voxel centers of the input images do not coincide with that
%         of the mask image.
%
%
%                           ----------------
%
% Non-sphericity correction
% =========================
%
% In some instances the i.i.d. assumptions about the errors do not hold.
% For example, in 2nd-level analyses serial correlations can be expressed
% as correlations among contrasts.  If you request a non-sphericity
% correction you will be asked to specify whether the repeated measures
% are correlated over replications.  It is important to ensure that the
% number of repeated measures (e.g. conditions) is not large in relation
% to the number of replications (e.g. subjects).  If you say 'no' to correlated
% repeated measures, then SPM will assume unequal variances for each level
% of each factor.
% The non-sphericity option is offered even if there is only one factor.
% This allows you to model unequal variances among its levels.
% The ensuing covariance components will be estimated using ReML in spm_spm
% (assuming the same for all voxels) and used to adjust the statistics and
% degrees of freedom during inference.  By default spm_spm will use wieghted
% least squares to produce Gauss-Markov or Maximum likelihood estimators
% using the non-sphericity strcuture specified at this stage. The components
% will be found in xX.xVi and enter the estimation procedure exactly as the 
% serial correlations in fMRI models.
% 
% 
%                           ----------------
%
% Multivariate analyses
% =====================
%
% Mulitvariate analyses with n-variate response variables are supported
% and automatically invoke a ManCova and CVA in spm_spm.  Multivariate
% designs are, at the moment limited to Basic and PET designs.
%
% ----------------------------------------------------------------------
%
% Variables saved in the SPM stucture
%
% xY.VY         - nScan x 1 struct array of memory mapped images
%                 (see spm_vol for definition of the map structure)
% xX            - structure describing design matrix
% xX.D          - design definition structure
%                 (See definition in main body of cg_spm_ui.m)
% xX.I          - nScan x 4 matrix of factor level indicators
%                 I(n,i) is the level of factor i corresponding to image n
% xX.sF         - 1x4 cellstr containing the names of the four factors
%                 xX.sF{i} is the name of factor i
% xX.X          - design matrix
% xX.xVi        - correlation constraints for non-spericity correction
% xX.iH         - vector of H partition (condition effects) indices,
%                 identifying columns of X correspoding to H
% xX.iC         - vector of C partition (covariates of interest) indices
% xX.iB         - vector of B partition (block effects) indices
% xX.iG         - vector of G partition (nuisance variables) indices
% xX.name     - p x 1 cellstr of effect names corresponding to columns
%                 of the design matrix
% 
% xC            - structure array of covariate details
% xC(i).rc      - raw (as entered) i-th covariate
% xC(i).rcname  - name of this covariate (string)
% xC(i).c       - covariate as appears in design matrix (after any scaling,
%                 centering of interactions)
% xC(i).cname   - cellstr containing names for effects corresponding to
%                 columns of xC(i).c
% xC(i).iCC     - covariate centering option
% xC(i).iCFI    - covariate by factor interaction option
% xC(i).type    - covariate type: 1=interest, 2=nuisance, 3=global
% xC(i).cols    - columns of design matrix corresponding to xC(i).c
% xC(i).descrip - cellstr containing a description of the covariate
% 
% xGX           - structure describing global options and values
% xGX.iGXcalc   - global calculation option used
% xGX.sGXcalc   - string describing global calculation used
% xGX.rg        - raw globals (before scaling and such like)
% xGX.iGMsca    - grand mean scaling option
% xGX.sGMsca    - string describing grand mean scaling
% xGX.GM        - value for grand mean (/proportional) scaling
% xGX.gSF       - global scaling factor (applied to xGX.rg)
% xGX.iGC       - global covariate centering option
% xGX.sGC       - string describing global covariate centering option
% xGX.gc        - center for global covariate
% xGX.iGloNorm  - Global normalisation option
% xGX.sGloNorm  - string describing global normalisation option
% 
% xM            - structure describing masking options
% xM.T          - Threshold masking value (-Inf=>None,
%                 real=>absolute, complex=>proportional (i.e. times global) )
% xM.TH         - nScan x 1 vector of analysis thresholds, one per image
% xM.I          - Implicit masking (0=>none, 1=>implicit zero/NaN mask)
% xM.VM         - struct array of explicit mask images
%                 (empty if no explicit masks)
% xM.xs         - structure describing masking options
%                 (format is same as for xsDes described below)
% 
% xsDes         - structure of strings describing the design:
%                 Fieldnames are essentially topic strings (use "_"'s for
%                 spaces), and the field values should be strings or cellstr's
%                 of information regarding that topic. spm_DesRep.m
%                 uses this structure to produce a printed description
%                 of the design, displaying the fieldnames (with "_"'s 
%                 converted to spaces) in bold as topics, with
%                 the corresponding text to the right
% 
% SPMid         - String identifying SPM and program versions
%
%                           ----------------
%
% NB: The SPM.mat file is not very portable: It contains
% memory-mapped handles for the images, which hardcodes the full file
% pathname and datatype. Therefore, subsequent to creating the
% SPM.mat, you cannot move the image files, and cannot move the
% entire analysis to a system with a different byte-order (even if the
% full file pathnames are retained. Further, the image scalefactors
% will have been pre-scaled to effect any grand mean or global
% scaling.
%_______________________________________________________________________
% modified version of spm_spm_ui 2.49 Andrew Holmes 03/03/20
% @(#)cg_spm_ui.m	1.02 Christiaqn Gaser 2006/07/31
SCCSid  = '1.02';

%-Condition arguments
%-----------------------------------------------------------------------
if (nargin==0), Action = 'CFG'; else, Action = varargin{1}; end


switch lower(Action), case 'cfg'
%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
% cg_spm_ui('CFG',D)
if nargin<2, D = []; else, D = varargin{2}; end

%-GUI setup
%-----------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: Setup analysis',0);
spm_help('!ContextHelp',mfilename)


%-Ask about overwriting files from previous analyses...
%-----------------------------------------------------------------------
if exist(fullfile('.','SPM.mat'))
	str = {	'Current directory contains existing SPM file:',...
		'Continuing will overwrite existing file!'};
	if spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename);
		fprintf('%-40s: %30s\n\n',...
			'Abort...   (existing SPM file)',spm('time'))
		spm_clf(Finter)
		return
	end
end


%-Option definitions
%-----------------------------------------------------------------------
%-Generic factor names
sF = {'sF1','sF2','sF3','sF4'};

%-Covariate by factor interaction options
sCFI = {'<none>';...							%-1
	'with sF1';'with sF2';'with sF3';'with sF4';...			%-2:5
	'with sF2 (within sF4)';'with sF3 (within sF4)'};		%-6,7

%-DesMtx argument components for covariate by factor interaction options
% (Used for CFI's Covariate Centering (CC), GMscale & Global normalisation)
CFIforms = {	'[]',		'C',	'{}';...			%-1
		'I(:,1)',	'FxC',	'{D.sF{1}}';...			%-2
		'I(:,2)',	'FxC',	'{D.sF{2}}';...			%-3
		'I(:,3)',	'FxC',	'{D.sF{3}}';...			%-4
		'I(:,4)',	'FxC',	'{D.sF{4}}';...			%-5
		'I(:,[4,2])',	'FxC',	'{D.sF{4},D.sF{2}}';...		%-6
		'I(:,[4,3])',	'FxC',	'{D.sF{4},D.sF{3}}'	};	%-7

%-Centre (mean correction) options for covariates & globals            (CC)
% (options 9-12 are for centering of global when using AnCova GloNorm) (GC)
sCC = {		'around overall mean';...				%-1
		'around sF1 means';...					%-2
		'around sF2 means';...					%-3
		'around sF3 means';...					%-4
		'around sF4 means';...					%-5
		'around sF2 (within sF4) means';...			%-6
		'around sF3 (within sF4) means';...			%-7
		'<no centering>';...					%-8
		'around user specified value';...			%-9
		'(as implied by AnCova)';...				%-10
		'GM';...						%-11
		'(redundant: not doing AnCova)'}';			%-12
%-DesMtx I forms for covariate centering options
CCforms = {'ones(nScan,1)',CFIforms{2:end,1},''}';


%-Global normalization options (options 1-7 match CFIforms)       (GloNorm)
sGloNorm = {	'AnCova';...						%-1
		'AnCova by sF1';...					%-2
		'AnCova by sF2';...					%-3
		'AnCova by sF3';...					%-4
		'AnCova by sF4';...					%-5
		'AnCova by sF2 (within sF4)';...			%-6
		'AnCova by sF3 (within sF4)';...			%-7
		'proportional scaling';...				%-8
		'<no global normalisation>'};				%-9

%-Grand mean scaling options                                        (GMsca)
sGMsca = {	'scaling of overall grand mean';...			%-1
		'scaling of sF1 grand means';...			%-2
		'scaling of sF2 grand means';...			%-3
		'scaling of sF3 grand means';...			%-4
		'scaling of sF4 grand means';...			%-5
		'scaling of sF2 (within sF4) grand means';...		%-6
		'scaling of sF3 (within sF4) grand means';...		%-7
		'(implicit in PropSca global normalisation)';...	%-8
		'<no grand Mean scaling>'	};			%-9
%-NB: Grand mean scaling by subject is redundent for proportional scaling


%-Global calculation options                                       (GXcalc)
sGXcalc  = {	'omit';...						%-1
		'user specified';...					%-2
		'mean voxel value (within per image fullmean/8 mask)'};	%-3



%=======================================================================
%-D E S I G N   P A R A M E T E R S
%=======================================================================
%-Get design type
%-----------------------------------------------------------------------
if isempty(D)
	D = cg_spm_ui( ...
		char(spm_input('Select design class...','+1','m',...
		{'Cross-sectional','Cross-sectional different style','Cross sectional asymmetry','Longitudinal'},...
		{'cohort_stats','cohort_stats2','asymmetry_stats','longitudinal_stats'},1)));
end

D = D(spm_input('Select design type...','+1','m',{D.DesName}'));


%-Set factor names for this design
%-----------------------------------------------------------------------
sCC      = sf_estrrep(sCC,[sF',D.sF']);
sCFI     = sf_estrrep(sCFI,[sF',D.sF']);
sGloNorm = sf_estrrep(sGloNorm,[sF',D.sF']);
sGMsca   = sf_estrrep(sGMsca,[sF',D.sF']);

%-Get filenames & factor indicies
%-----------------------------------------------------------------------
[P,I]    = cg_spm_ui('Files&Indices',D.sF,D.n,D.b.aTime);
nScan    = size(I,1);						%-#obs

%-Additional design parameters
%-----------------------------------------------------------------------
bL       = any(diff(I,1),1); 	%-Multiple factor levels?
	    % NB: bL(2) might be thrown by user specified f1 levels
	    %     (D.b.aTime & D.n(2)>1) - assumme user is consistent?
bFI      = [bL(1),bL(2:3)&~bL(4),bL(4),bL([2,3])&bL(4)];
	    %-Allowable interactions for covariates
	    %-Only offer interactions with multi-level factors, and
	    % don't offer by F2|F3 if bL(4)!

%-Build Condition (H) and Block (B) partitions
%=======================================================================
eval(['[H,Hnames] = spm_DesMtx(',D.Hform,');'])
if rank(H)==nScan, error('unestimable condition effects'), end
eval(['[B,Bnames] = spm_DesMtx(',D.Bform,');'])
if rank(B)==nScan, error('unestimable block effects'), end

%-Drop a constant H partition if B partition can model constant
if size(H,2)>0 & all(H(:)==1) & (rank([H B])==rank(B))
	H = []; Hnames = {};
	warning('Dropping redundant constant H partition')
end


%-Covariate partition(s): interest (C) & nuisance (G) excluding global
%=======================================================================
nC = D.nC;			%-Default #covariates
C  = {[],[]}; Cnames = {{},{}};	%-Covariate DesMtx partitions & names
xC = [];			%-Struct array to hold raw covariates


dcname = {'CovInt','NusCov'};	%-Default root names for covariates
dstr   = {'covariate','nuisance variable'};

GUIpos = spm_input('!NextPos');
nc     = [0,0];
for i  = 1:2			% 1:covariates of interest, 2:nuisance variables

    if isinf(nC(i)), nC(i)=spm_input(['# ',dstr{i},'s'],GUIpos,'w1'); end

    while nc(i) < nC(i)

	%-Create prompt, get covariate, get covariate name
        %---------------------------------------------------------------
	if nC(i)==1, str=dstr{i}; else, str=sprintf('%s %d',dstr{i},nc(i)+1); end
        c = spm_input(str,GUIpos,'r',[],[nScan,Inf]);
        if any(isnan(c(:))), break, end		%-NaN is dummy value to exit
	nc(i)  = nc(i)+1;			%-#Covariates (so far)
	if nC(i)>1,	tstr = sprintf('%s^{%d}',dcname{i},nc(i));
	else,		tstr = dcname{i}; end
       	cname  = spm_input([str,' name?'],'+1','s',tstr);
       	rc     = c;				%-Save covariate value
	rcname = cname;				%-Save covariate name

        %-Interaction option? (if single covariate vector entered)?
        %---------------------------------------------------------------
        if size(c,2) == 1
       	    if length(D.iCFI{i})>1
       		%-User choice of interaction options, default is negative
       		%-Only offer interactions for appropriate factor combinations
		iCFI = intersect(abs(D.iCFI{i}),find([1,bFI]));
		dCFI = max([1,intersect(iCFI,-D.iCFI{i}(D.iCFI{i}<0))]);
        	iCFI = spm_input([str,': interaction?'],'+1','m',...
			sCFI(iCFI),iCFI,find(iCFI==dCFI));
	    else
		iCFI = abs(D.iCFI{i});		%-AutoSelect default option
	    end
	else
	    iCFI = 1;
	end

        %-Centre covariate(s)? (Default centring to correspond to CFI)
        % Always offer "no centering" as default for design matrix blocks
        %---------------------------------------------------------------
	DiCC = D.iCC{i};
	if size(c,2)>1, DiCC = union(DiCC,-8); end
        if length(DiCC)>1
        	%-User has a choice of centering options
		%-Only offer factor specific for appropriate factor combinations
		iCC = intersect(abs(DiCC),find([1,bFI,1]) );
        	%-Default is max -ve option in D, overridden by iCFI if CFI
		if iCFI == 1, dCC = -DiCC(DiCC<0); else, dCC = iCFI; end
		dCC = max([1,intersect(iCC,dCC)]);
		iCC = spm_input([str,': centre?'],'+1','m',...
			sCC(iCC),iCC,find(iCC==dCC));
        else
        	iCC = abs(DiCC);	%-AutoSelect default option
        end
	%-Centre within factor levels as appropriate
        if any(iCC == [1:7]), c = c - spm_meanby(c,eval(CCforms{iCC})); end

        %-Do any interaction (only for single covariate vectors)
        %---------------------------------------------------------------
        if iCFI > 1				%-(NB:iCFI=1 if size(c,2)>1)
       		tI        = [eval(CFIforms{iCFI,1}),c];
		tConst    = CFIforms{iCFI,2};
		tFnames   = [eval(CFIforms{iCFI,3}),{cname}];
		[c,cname] = spm_DesMtx(tI,tConst,tFnames);
	elseif size(c,2)>1			%-Design matrix block
		[null,cname] = spm_DesMtx(c,'X',cname);
	else
		cname = {cname};
	end

	%-Store raw covariate details in xC struct for reference
	%-Pack c into appropriate DesMtx partition
        %---------------------------------------------------------------
	%-Construct description string for covariate
	str = {sprintf('%s: %s',str,rcname)};
	if size(rc,2)>1, str = {sprintf('%s (block of %d covariates)',...
		str{:},size(rc,2))}; end
	if iCC < 8, str=[str;{['used centered ',sCC{iCC}]}]; end
	if iCFI> 1, str=[str;{['fitted as interaction ',sCFI{iCFI}]}]; end

	tmp       = struct(	'rc',rc,	'rcname',rcname,...
				'c',c,		'cname',{cname},...
				'iCC',iCC,	'iCFI',iCFI,...
				'type',i,...
				'cols',[1:size(c,2)] + ...
						size([H,C{1}],2) +  ...
						size([B,C{2}],2)*(i-1),...
				'descrip',{str}				);
	if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end
	C{i}      = [C{i},c];
	Cnames{i} = [Cnames{i}; cname];

    end	% (while)

end % (for)
clear c tI tConst tFnames
spm_input('!SetNextPos',GUIpos);

%-Unpack into C & G design matrix sub-partitions
G = C{2}; Gnames = Cnames{2};
C = C{1}; Cnames = Cnames{1};


%-Options...
%=======================================================================
%-Global normalization options                                 (GloNorm)
%-----------------------------------------------------------------------
if length(D.iGloNorm)>1
	%-User choice of global normalisation options, default is negative
	%-Only offer factor specific for appropriate factor combinations
	iGloNorm = intersect(abs(D.iGloNorm),find([1,bFI,1,1]));
	dGloNorm = max([0,intersect(iGloNorm,-D.iGloNorm(D.iGloNorm<0))]);
	iGloNorm = spm_input('GloNorm: Select global normalisation','+1','m',...
	    	sGloNorm(iGloNorm),iGloNorm,find(iGloNorm==dGloNorm));
else
	iGloNorm = abs(D.iGloNorm);
end


%-Grand mean scaling options                                     (GMsca)
%-----------------------------------------------------------------------
if iGloNorm==8
	iGMsca=8;	%-grand mean scaling implicit in PropSca GloNorm
elseif length(D.iGMsca)==1
	iGMsca = abs(D.iGMsca);
else
	%-User choice of grand mean scaling options
	%-Only offer factor specific for appropriate factor combinations
	iGMsca = intersect(abs(D.iGMsca),find([1,bFI,0,1]));
        %-Default is max -ve option in D, overridden by iGloNorm if AnCova
        if iGloNorm==9, dGMsca=-D.iGMsca(D.iGMsca<0); else, dGMsca=iGloNorm; end
	dGMsca = max([0,intersect(iGMsca,dGMsca)]);
	iGMsca = spm_input('GMsca: grand mean scaling','+1','m',...
	    	sGMsca(iGMsca),iGMsca,find(iGMsca==dGMsca));
end


%-Value for PropSca / GMsca                                         (GM)
%-----------------------------------------------------------------------
if iGMsca == 9                          %-Not scaling (GMsca or PropSca)
	GM = 0;                         %-Set GM to zero when not scaling
else                                    %-Ask user value of GM
	if iGloNorm==8
		str = 'PropSca global mean to';
	else
		str = [strrep(sGMsca{iGMsca},'scaling of','scale'),' to'];
	end
	GM = spm_input(str,'+1','r',D.GM,1);
	%-If GM is zero then don't GMsca! or PropSca GloNorm
	if GM==0, iGMsca=9; if iGloNorm==8, iGloNorm=9; end, end
end

%-Sort out description strings for GloNorm and GMsca
%-----------------------------------------------------------------------
sGloNorm = sGloNorm{iGloNorm};
sGMsca   = sGMsca{iGMsca};
if iGloNorm==8
	sGloNorm = sprintf('%s to %-4g',sGloNorm,GM);
elseif iGMsca<8
	sGMsca   = sprintf('%s to %-4g',sGMsca,GM);
end


%-Global centering (for AnCova GloNorm)                             (GC)
%-----------------------------------------------------------------------
%-Specify the centering option for the global covariate for AnCova
%-Basically, if 'GMsca'ling then should centre to GM (iGC=11). Otherwise,
% should centre in similar fashion to AnCova (i.e. by the same factor(s)),
% such that models are seperable (iGC=10). This is particularly important
% for subject specific condition effects if then passed on to a second-level
% model. (See also spm_adjmean_ui.m) SPM96 (& earlier) used to just centre
% GX around its (overall) mean (iGC=1).

%-This code allows more general options to be specified (but is a bit complex)
%-Setting D.iGC=[-10,-11] gives the standard choices above

%-If not doing AnCova then GC is irrelevant
if ~any(iGloNorm == [1:7])
	iGC = 12;
	gc  = [];
else
	%-Annotate options 10 & 11 with specific details
	%---------------------------------------------------------------
	%-Tag '(as implied by AnCova)' with actual AnCova situation
	sCC{10} = [sCC{iGloNorm},' (<= ',sGloNorm,')'];
	%-Tag 'GM' case with actual GM & GMsca case
	sCC{11} = sprintf('around GM=%g (i.e. %s after grand mean scaling)',...
		GM,strrep(sCC{iGMsca},'around ',''));

	%-Constuct vector of allowable iGC
	%---------------------------------------------------------------
	%-Weed out redundent factor combinations from pre-set allowable options
	iGC = intersect(abs(D.iGC),find([1,bFI,1,1,1,1]));
	%-Omit 'GM' option if didn't GMsca (iGMsca~=8 'cos doing AnCova)
	if any(iGMsca==[8,9]), iGC = setdiff(iGC,11); end
	%-Omit 'GM' option if same as '(as implied by AnCova)'
	if iGloNorm==iGMsca, iGC = setdiff(iGC,11); end

	%-If there's a choice, set defaults (if any), & get answer
	%---------------------------------------------------------------
	if length(iGC)>1
		dGC = max([0,intersect(iGC,-D.iGC(D.iGC<0))]);
		str = 'Centre global covariate';
		if iGMsca<8, str = [str,' (after grand mean scaling)']; end
		iGC = spm_input(str,'+1','m',sCC(iGC),iGC,find(iGC==dGC));
	elseif isempty(iGC)
		error('Configuration error: empty iGC')
	end

	%-If 'user specified' then get value
	%---------------------------------------------------------------
	if iGC==9
		gc     = spm_input('Centre globals around','+0','r',D.GM,1);
		sCC{9} = sprintf('%s of %g',sCC{iGC},gc);
	else
		gc  = 0;
	end
end


%-Thresholds & masks defining voxels to analyse                   (MASK)
%=======================================================================
GUIpos = spm_input('!NextPos');

%-Analysis threshold mask
%-----------------------------------------------------------------------
%-Work out available options:
% -Inf=>None, real=>absolute, complex=>proportional, (i.e. times global)
M_T = D.M_.T; if isempty(M_T), M_T = [-Inf, 100, 0.8*sqrt(-1)]; end
M_T = {	'none',		M_T(min(find(isinf(M_T))));...
	'absolute',	M_T(min(find(isfinite(M_T)&(M_T==real(M_T)))));...
	'relative',	M_T(min(find(isfinite(M_T)&(M_T~=real(M_T)))))	};

%-Work out available options
%-If there's a choice between proportional and absolute then ask
%-----------------------------------------------------------------------
q = ~[isempty(M_T{1,2}), isempty(M_T{2,2}), isempty(M_T{3,2})];

if all(q(2:3))
	tmp = spm_input('Threshold masking',GUIpos,'b',M_T(q,1),find(q));
	q(setdiff([1:3],tmp))=0;
end

%-Get mask value - note that at most one of q(2:3) is true
%-----------------------------------------------------------------------
if ~any(q)				%-Oops - nothing specified!
	M_T = -Inf;
elseif all(q==[1,0,0])			%-no threshold masking
	M_T = -Inf;
else					%-get mask value
	if q(1),	args = {'br1','None',-Inf,abs(M_T{1+find(q(2:3)),2})};
	else,		args = {'r',abs(M_T{1+find(q(2:3)),2})}; end
	if q(2)
		M_T = spm_input('threshold',GUIpos,args{:});
	elseif q(3)
		M_T = spm_input('threshold (relative to global)',GUIpos,...
								args{:});
		if isfinite(M_T) & isreal(M_T), M_T=M_T*sqrt(-1); end
	else
		error('Shouldn''t get here!')
	end
end

%-Make a description string
%-----------------------------------------------------------------------
if isinf(M_T)
	xsM.Analysis_threshold = 'None (-Inf)';
elseif isreal(M_T)
	xsM.Analysis_threshold = sprintf('images thresholded at %6g',M_T);
else
	xsM.Analysis_threshold = sprintf(['images thresholded at %6g ',...
		'times global'],imag(M_T));
end


%-Implicit masking: Ignore zero voxels in low data-types?
%-----------------------------------------------------------------------
% (Implicit mask is NaN in higher data-types.)
type = getfield(spm_vol(P{1,1}),'dim')*[0,0,0,1]';
if ~spm_type(type,'nanrep')
	switch D.M_.I
	case Inf,    M_I = spm_input('Implicit mask (ignore zero''s)?',...
			'+1','y/n',[1,0],1);		%-Ask
	case {0,1}, M_I = D.M_.I;			%-Pre-specified
	otherwise,  error('unrecognised D.M_.I type')
	end

	if M_I, xsM.Implicit_masking = 'Yes: zero''s treated as missing';
	else,   xsm.Implicit_masking = 'No'; end
else
	M_I = 1;
	xsM.Implicit_masking = 'Yes: NaN''s treated as missing';
end


%-Explicit mask images (map them later...)
%-----------------------------------------------------------------------
switch(D.M_.X)
case Inf,   M_X = spm_input('explicitly mask images?','+1','y/n',[1,0],2);
case {0,1}, M_X = D.M_.X;
otherwise,  error('unrecognised D.M_.X type')
end
if M_X, M_P = spm_get(Inf,'*.img',{'select mask images'}); else, M_P = {}; end


%-Global calculation                                            (GXcalc)
%=======================================================================
iGXcalc = abs(D.iGXcalc);
%-Only offer "omit" option if not doing any GloNorm, GMsca or PropTHRESH
if ~(iGloNorm==9 & iGMsca==9 & (isinf(M_T)|isreal(M_T)))
	iGXcalc = intersect(iGXcalc,[2:size(sGXcalc,1)]);
end
if isempty(iGXcalc)
	error('no GXcalc options')
elseif length(iGXcalc)>1
	%-User choice of global calculation options, default is negative
	dGXcalc = max([1,intersect(iGXcalc,-D.iGXcalc(D.iGXcalc<0))]);
	iGXcalc = spm_input('Global calculation','+1','m',...
	    	sGXcalc(iGXcalc),iGXcalc,find(iGXcalc==dGXcalc));
else
	iGXcalc = abs(D.iGXcalc);
end

if iGXcalc==2				%-Get user specified globals
	g = spm_input('globals','+0','r',[],[nScan,1]);
end
sGXcalc = sGXcalc{iGXcalc};


% Non-sphericity correction
%=======================================================================

% if there are multilevel factors, ask for correction
%-----------------------------------------------------------------------
if length(find(max(I) > 1)) > 1
	xVi.iid  = spm_input('non-sphericity correction?','+1','y/n',[0,1],1);
else
	xVi.iid  = 1;
end


if xVi.iid

	% i.i.d. assumptions where xVi.V = 1
	%---------------------------------------------------------------
	xVi.V    = speye(nScan);

else
	% otherwise, we have repeated measures design 
	%===============================================================
	nL      = max(I);		% number of levels
	mL      = find(nL > 1);		% multilevel factors
	xVi.I   = I;
	xVi.sF  = D.sF;
	xVi.var = sparse(1,4);
	xVi.dep = sparse(1,4);


	% eliminate replication factor from mL
	%---------------------------------------------------------------
	for i = 1:4
		mstr{i} = sprintf('%s (%i)',D.sF{i},nL(i));
	end
	str   = 'replications are over?';
	rep   = spm_input(str,'+1','m',mstr(mL),1:length(mL));

	% and ask whether repeated measures are independent
	%---------------------------------------------------------------
	str   = 'correlated repeated measures';
	dep   = spm_input(str,'+1','b',{'yes','no'},[1 0],2);


	%-Place covariance components Q{:} in xVi.Vi
	%---------------------------------------------------------------
	mL(rep)     = [];
	xVi.var(mL) = 1;
	xVi.dep(mL) = dep;
	xVi         = spm_non_sphericity(xVi);

end


%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
spm('FigName','Stats: configuring',Finter,CmdLine);
spm('Pointer','Watch');


%-Images & image info: Map Y image files and check consistency of
% dimensions and orientation / voxel size
%=======================================================================
fprintf('%-40s: ','Mapping files')                                   %-#
VY    = spm_vol(char(P));


%-Check compatability of images (Bombs for single image)
%-----------------------------------------------------------------------
if any(any(diff(cat(1,VY(:).dim),1,1),1) & [1,1,1,0]) 
	error('images do not all have the same dimensions')
end
if any(any(any(diff(cat(3,VY(:).mat),1,3),3)))
	error('images do not all have same orientation & voxel size')
end

fprintf('%30s\n','...done')                                          %-#


%-Global values, scaling and global normalisation
%=======================================================================
%-Compute global values
%-----------------------------------------------------------------------
switch iGXcalc, case 1
	%-Don't compute => no GMsca (iGMsca==9) or GloNorm (iGloNorm==9)
	g = [];
case 2
	%-User specified globals
case 3
	%-Compute as mean voxel value (within per image fullmean/8 mask)
	g     = zeros(nScan,1 );
	fprintf('%-40s: %30s','Calculating globals',' ')             %-#
	for i = 1:nScan
		str = sprintf('%3d/%-3d',i,nScan);
		fprintf('%s%30s',sprintf('\b')*ones(1,30),str)%-#
		g(i) = spm_global(VY(i));
	end
	fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')       %-#
otherwise
	error('illegal iGXcalc')
end
rg = g;


fprintf('%-40s: ','Design configuration')                            %-#


%-Scaling: compute global scaling factors gSF required to implement proportional
% scaling global normalisation (PropSca) or grand mean scaling (GMsca),
% as specified by iGMsca (& iGloNorm)
%-----------------------------------------------------------------------
switch iGMsca, case 8
	%-Proportional scaling global normalisation
	if iGloNorm~=8, error('iGloNorm-iGMsca(8) mismatch for PropSca'), end
	gSF    = GM./g;
	g      = GM*ones(nScan,1);
case {1,2,3,4,5,6,7}
	%-Grand mean scaling according to iGMsca
	gSF    = GM./spm_meanby(g,eval(CCforms{iGMsca}));
	g      = g.*gSF;
case 9
	%-No grand mean scaling
	gSF    = ones(nScan,1);
otherwise
	error('illegal iGMsca')
end


%-Apply gSF to memory-mapped scalefactors to implement scaling
%-----------------------------------------------------------------------
for i = 1:nScan
	VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*gSF(i);
end


%-AnCova: Construct global nuisance covariates partition (if AnCova)
%-----------------------------------------------------------------------
if any(iGloNorm == [1:7])

	%-Centre global covariate as requested
	%---------------------------------------------------------------
	switch iGC, case {1,2,3,4,5,6,7}	%-Standard sCC options
		gc = spm_meanby(g,eval(CCforms{iGC}));
	case 8					%-No centering
		gc = 0;
	case 9					%-User specified centre
		%-gc set above
	case 10					%-As implied by AnCova option
		gc = spm_meanby(g,eval(CCforms{iGloNorm}));
	case 11					%-Around GM
		gc = GM;
	otherwise				%-unknown iGC
		error('unexpected iGC value')
	end
	
	
	%-AnCova - add scaled centred global to DesMtx `G' partition
	%---------------------------------------------------------------
	rcname     = 'global'; 
	tI         = [eval(CFIforms{iGloNorm,1}),g - gc];
	tConst     = CFIforms{iGloNorm,2};
	tFnames    = [eval(CFIforms{iGloNorm,3}),{rcname}];
	[f,gnames]  = spm_DesMtx(tI,tConst,tFnames);
	clear tI tConst tFnames
	
	%-Save GX info in xC struct for reference
	%---------------------------------------------------------------
	str     = {sprintf('%s: %s',dstr{2},rcname)};
	if any(iGMsca==[1:7]), str=[str;{['(after ',sGMsca,')']}]; end
	if iGC ~= 8, str=[str;{['used centered ',sCC{iGC}]}]; end
	if iGloNorm > 1
		str=[str;{['fitted as interaction ',sCFI{iGloNorm}]}]; 
	end
	tmp  = struct(	'rc',rg.*gSF,		'rcname',rcname,...
			'c',f,			'cname'	,{gnames},...
			'iCC',iGC,		'iCFI'	,iGloNorm,...
			'type',			3,...
			'cols',[1:size(f,2)] + size([H C B G],2),...
				'descrip',		{str}		);

	G = [G,f]; Gnames = [Gnames; gnames];
	if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end


elseif iGloNorm==8 | iGXcalc>1

	%-Globals calculated, but not AnCova: Make a note of globals
	%---------------------------------------------------------------
	if iGloNorm==8
		str = { 'global values: (used for proportional scaling)';...
			'("raw" unscaled globals shown)'};
	elseif isfinite(M_T) & ~isreal(M_T)
		str = { 'global values: (used to compute analysis threshold)'};
	else
		str = { 'global values: (computed but not used)'};
	end

	rcname ='global';
	tmp     = struct(	'rc',rg,	'rcname',rcname,...
				'c',{[]},	'cname'	,{{}},...
				'iCC',0,	'iCFI'	,0,...
				'type',		3,...
				'cols',		{[]},...
				'descrip',	{str}			);

	if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end
end


%-Save info on global calculation in xGX structure
%-----------------------------------------------------------------------
xGX = struct(...
	'iGXcalc',iGXcalc,	'sGXcalc',sGXcalc,	'rg',rg,...
	'iGMsca',iGMsca,	'sGMsca',sGMsca,	'GM',GM,'gSF',gSF,...
	'iGC',	iGC,		'sGC',	sCC{iGC},	'gc',	gc,...
	'iGloNorm',iGloNorm,	'sGloNorm',sGloNorm);



%-Construct masking information structure and compute actual analysis
% threshold using scaled globals (rg.*gSF)
%-----------------------------------------------------------------------
if isreal(M_T),	M_TH =      M_T  * ones(nScan,1);	%-NB: -Inf is real
else,		M_TH = imag(M_T) * (rg.*gSF); end

if ~isempty(M_P)
	VM  = spm_vol(char(M_P));
	xsM.Explicit_masking = [{'Yes: mask images :'};{VM.fname}'];
else
	VM  = [];
	xsM.Explicit_masking = 'No';
end
xM     = struct('T',M_T, 'TH',M_TH, 'I',M_I, 'VM',{VM}, 'xs',xsM);


%-Construct full design matrix (X), parameter names and structure (xX)
%=======================================================================
X      = [H C B G];
tmp    = cumsum([size(H,2), size(C,2), size(B,2), size(G,2)]);
xX     = struct(	'X',		X,...
			'iH',		[1:size(H,2)],...
			'iC',		[1:size(C,2)] + tmp(1),...
			'iB',		[1:size(B,2)] + tmp(2),...
			'iG',		[1:size(G,2)] + tmp(3),...
			'name',		{[Hnames; Cnames; Bnames; Gnames]},...
			'I',		I,...
			'sF',		{D.sF});


%-Design description (an nx2 cellstr) - for saving and display
%=======================================================================
tmp = {	sprintf('%d condition, +%d covariate, +%d block, +%d nuisance',...
		size(H,2),size(C,2),size(B,2),size(G,2));...
	sprintf('%d total, having %d degrees of freedom',...
		size(X,2),rank(X));...
	sprintf('leaving %d degrees of freedom from %d images',...
		size(X,1)-rank(X),size(X,1))				};
xsDes = struct(	'Design',			{D.DesName},...
		'Global_calculation',		{sGXcalc},...
		'Grand_mean_scaling',		{sGMsca},...
		'Global_normalisation',		{sGloNorm},...
		'Parameters',			{tmp}			);


fprintf('%30s\n','...done')                                          %-#



%-Assemble SPM structure
%=======================================================================
SPM.xY.P	= P;			% filenames
SPM.xY.VY	= VY;			% mapped data
SPM.nscan	= size(xX.X,1);		% scan number
SPM.xX		= xX;			% design structure
SPM.xC		= xC;			% covariate structure
SPM.xGX		= xGX;			% global structure
SPM.xVi		= xVi;			% non-sphericity structure
SPM.xM		= xM;			% mask structure
SPM.xsDes	= xsDes;		% description
SPM.SPMid	= SPMid;		% version

%-Save SPM.mat and set output argument
%-----------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')                        %-#
save SPM SPM;
fprintf('%30s\n','...SPM.mat saved')                                 %-#
varargout = {SPM};

%-Display Design report
%=======================================================================
fprintf('%-40s: ','Design reporting')                                %-#
fname     = cat(1,{SPM.xY.VY.fname}');
spm_DesRep('DesMtx',SPM.xX,fname,SPM.xsDes)
fprintf('%30s\n','...done')     


%-End: Cleanup GUI
%=======================================================================
spm_clf(Finter)
spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                     %-#
spm('FigName','Stats: configured',Finter,CmdLine);
spm('Pointer','Arrow')
fprintf('\n\n')



case 'files&indices'
%=======================================================================
% - Get files and factor indices
%=======================================================================
% [P,I] = cg_spm_ui('Files&Indices',DsF,Dn,DbaTime,nV)
% DbaTime=D.b.aTime; Dn=D.n; DsF=D.sF;
if nargin<5, nV = 1; else, nV = varargin{5}; end
if nargin<4, DbaTime = 1; else, DbaTime = varargin{4}; end
if nargin<3, Dn  = [Inf,Inf,Inf,Inf]; else, Dn=varargin{3}; end
if nargin<2, DsF = {'Fac1','Fac2','Fac3','Fac4'}; else, DsF=varargin{2}; end

%-Initialise variables
%-----------------------------------------------------------------------
i4 = [];		% factor 4 index (usually group)
i3 = [];		% factor 3 index (usually subject), per f4
i2 = [];		% factor 2 index (usually condition), per f3/f4
i1 = [];		% factor 1 index (usually replication), per f2/f3/f4
P  = {};		% cell array of string filenames

%-Accrue filenames and factor level indicator vectors
%-----------------------------------------------------------------------
bMV = nV>1;
if isinf(Dn(4)), n4 = spm_input(['#',DsF{4},'''s'],'+1','n1');
	else, n4 = Dn(4); end
bL4 = n4>1;

ti2 = '';
GUIpos = spm_input('!NextPos');
for j4  = 1:n4
    spm_input('!SetNextPos',GUIpos);
    sF4P=''; if bL4, sF4P=[DsF{4},' ',int2str(j4),': ']; end
    if isinf(Dn(3)), n3=spm_input([sF4P,'#',DsF{3},'''s'],'+1','n1');
	    else, n3 = Dn(3); end
    bL3 = n3>1;
    
    if DbaTime & Dn(2)>1
	%disp('NB:selecting in time order - manually specify conditions')
	%-NB: This means f2 levels might not be 1:n2
	GUIpos2 = spm_input('!NextPos');
	for j3 = 1:n3
	    sF3P=''; if bL3, sF3P=[DsF{3},' ',int2str(j3),': ']; end
	    str = [sF4P,sF3P];
	    tP  = {};
	    n21 = Dn(2)*Dn(1);
	    for v=1:nV
	    	vstr=''; if bMV, vstr=sprintf(' (var-%d)',v); end
	    	ttP = spm_get(n21,'.img',{[str,'select images',vstr]});
	    	n21 = length(ttP);
	    	tP  = [tP,ttP];
	    end
	    ti2 = spm_input([str,' ',DsF{2},'?'],GUIpos2,'c',ti2',n21,Dn(2));
	    %-Work out i1 & check
	    [tl2,null,j] = unique(ti2);
	    tn1 = zeros(size(tl2)); ti1 = zeros(size(ti2));
	    for i=1:length(tl2)
		    tn1(i)=sum(j==i); ti1(ti2==tl2(i))=1:tn1(i); end
	    if isfinite(Dn(1)) & any(tn1~=Dn(1))
		%-#i1 levels mismatches specification in Dn(1)
		error(sprintf('#%s not %d as pre-specified',DsF{1},Dn(1)))
	    end
	    P  = [P;tP];
	    i4 = [i4; j4*ones(n21,1)];
	    i3 = [i3; j3*ones(n21,1)];
	    i2 = [i2; ti2];
	    i1 = [i1; ti1];
	end

    else

	if isinf(Dn(2))
	    n2 = spm_input([sF4P,'#',DsF{2},'''s'],'+1','n1');
	else
	    n2 = Dn(2);
	end
	bL2 = n2>1;

	if n2==1 & Dn(1)==1 %-single scan per f3 (subj)
	    %disp('NB:single scan per f3')
	    str = [sF4P,'select images, ',DsF{3},' 1-',int2str(n3)];
	    tP = {};
	    for v=1:nV
	    	vstr=''; if bMV, vstr=sprintf(' (var-%d)',v); end
	    	ttP = spm_get(n3,'.img',{[str,vstr]});
	    	tP = [tP,ttP];
	    end
	    P   = [P;tP];
	    i4  = [i4; j4*ones(n3,1)];
	    i3  = [i3; [1:n3]'];
	    i2  = [i2; ones(n3,1)];
	    i1  = [i1; ones(n3,1)];
	else
	    %-multi scan per f3 (subj) case
	    %disp('NB:multi scan per f3')
	    for j3 = 1:n3
		sF3P=''; if bL3, sF3P=[DsF{3},' ',int2str(j3),': ']; end
		if Dn(1)==1
			%-No f1 (repl) within f2 (cond)
			%disp('NB:no f1 within f2')
			str = [sF4P,sF3P,'select images: ',DsF{2},...
				 ' 1-',int2str(n2)];
			tP = {};
			for v=1:nV
				vstr=''; if bMV, vstr=sprintf(' (var-%d)',v); end
				ttP = spm_get(n2,'.img',{[str,vstr]});
				tP = [tP,ttP];
			end
			P   = [P;tP];
			i4  = [i4; j4*ones(n2,1)];
			i3  = [i3; j3*ones(n2,1)];
			i2  = [i2; [1:n2]'];
			i1  = [i1; ones(n2,1)];
		else
		    %-multi f1 (repl) within f2 (cond)
		    %disp('NB:f1 within f2')
		    for j2  = 1:n2
			sF2P='';
			if bL2, sF2P=[DsF{2},' ',int2str(j2),': ']; end
			str = [sF4P,sF3P,sF2P,' select images...'];
			tP  = {};
			n1  = Dn(1);
			for v=1:nV
				vstr=''; if bMV, vstr=sprintf(' (var-%d)',v); end
				ttP = spm_get(n1,'.img',{[str,vstr]});
				n1  = length(ttP);
				tP  = [tP,ttP];
			end
			P   = [P;tP];
			i4  = [i4; j4*ones(n1,1)];
			i3  = [i3; j3*ones(n1,1)];
			i2  = [i2; j2*ones(n1,1)];
			i1  = [i1; [1:n1]'];
		    end                         % (for j2)
		end                             % (if Dn(1)==1)
	    end                                 % (for j3)
	end                                     % (if  n2==1 &...)
    end                                         % (if DbaTime & Dn(2)>1)
end                                             % (for j4)
varargout = {P,[i1,i2,i3,i4]};


case 'cohort_stats'

    D = struct(...
            'DesName','Cross-sectional design (1 scan per subject)',...
            'n',	[Inf Inf 1 1],	'sF',{{'repl','group','',''}},...
            'Hform',		'I(:,2),''-'',''group''',...
            'Bform',		'I(:,3),''-'',''\mu''',...
            'nC',[Inf,Inf],'iCC',{{[1,3,4,8],[1,4,8]}},'iCFI',{{[1,3,-4],[1,-4]}},...
            'iGXcalc',1,'iGMsca',9,'GM',[],...
            'iGloNorm',9,'iGC',12,...
            'M_',struct('T',0.1,'I',0,'X',Inf),...
            'b',struct('aTime',0));
            
varargout = {D};

case 'cohort_stats2'

    D = struct(...
            'DesName','Cross-sectional design (1 scan per subject)',...
            'n',	[Inf Inf 1 1],	'sF',{{'repl','group','',''}},...
            'Hform',		'I(:,2),''-'',''group''',...
            'Bform',		'I(:,3),''-'',''\mu''',...
            'nC',[Inf,Inf],'iCC',{{[1,3,4,8],[1,4,8]}},'iCFI',{{[1,3,-4],[1,-4]}},...
            'iGXcalc',1,'iGMsca',9,'GM',[],...
            'iGloNorm',9,'iGC',12,...
            'M_',struct('T',0.1,'I',0,'X',Inf),...
            'b',struct('aTime',1));

varargout = {D};

case 'asymmetry_stats'

    D = struct(...
            'DesName','Cross-sectional design asymmetry',...
            'n',	[Inf Inf 1 1],	'sF',{{'repl','group','',''}},...
            'Hform',		'I(:,2),''-'',''group''',...
            'Bform',		'I(:,3),''-'',''\mu''',...
            'nC',[Inf,Inf],'iCC',{{[1,3,4,8],[1,4,8]}},'iCFI',{{[1,3,-4],[1,-4]}},...
            'iGXcalc',1,'iGMsca',9,'GM',[],...
            'iGloNorm',9,'iGC',12,...
            'M_',struct('T',-Inf,'I',Inf,'X',Inf),...
            'b',struct('aTime',0));
            
varargout = {D};

case 'longitudinal_stats'
            
    D = struct(...
            'DesName','Longitudinal design: (more than 1 time point per subject)',...
            'n',[Inf Inf Inf Inf],	'sF',{{'repl','cond','subject','group'}},...
            'Hform',		'I(:,[4,2]),''-'',{''group'',''repl''}',...
            'Bform',		'I(:,[4,3]),''-'',{''group'',''subj''}',...
            'nC',[Inf,Inf],'iCC',{{[5:8],[5,7,8]}},'iCFI',{{[1,5,6,-7],[1,5,-7]}},...
            'iGXcalc',1,'iGMsca',9,'GM',[],...
            'iGloNorm',9,'iGC',12,...
            'M_',struct('T',0.1,'I',0,'X',Inf),...
            'b',struct('aTime',1));

varargout = {D};

otherwise
%=======================================================================
% - U N K N O W N   A C T I O N
%=======================================================================
warning(['Illegal Action string: ',Action])

%=======================================================================
% - E N D
%=======================================================================
end

fprintf('Use the ''Estimate'' button to estimate parameters\n');


%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function str = sf_estrrep(str,srstr)
%=======================================================================
for i = 1:size(srstr,1)
	str = strrep(str,srstr{i,1},srstr{i,2});
end
