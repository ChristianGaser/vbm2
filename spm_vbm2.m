function spm_vbm2
% VBM2 Toolbox wrapper to call vbm functions
%_______________________________________________________________________
% @(#)spm_vbm2.m	v1.08 Christian Gaser 2007/06/25

tbs = spm('tbs');
for i=1:length(tbs)
	if strcmp(lower(tbs(i).name),'vbm2')
		addpath(tbs(i).dir,'-begin');
	end
end

SPMid = spm('FnBanner',mfilename,'v1.08');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','VBM2');
spm_help('!ContextHelp',mfilename);
spm_help('!Disp','vbm2.man','',Fgraph,'     Voxel-based morphometry toolbox for SPM2');

fig = spm_figure('GetWin','Interactive');
h0  = uimenu(fig,...
	'Label',	'VBM',...
	'Separator',	'on',...
	'Tag',		'VBM',...
	'HandleVisibility','on');
h1  = uimenu(h0,...
	'Label',	'Create customized template',...
	'Separator',	'off',...
	'Tag',		'Create template',...
	'CallBack',	'cg_create_template;',...
	'HandleVisibility','on');
h2  = uimenu(h0,...
	'Label',	'Segment',...
	'Separator',	'on',...
	'Tag',		'Segmentation',...
	'HandleVisibility','on');
h21  = uimenu(h2,...
	'Label',	'Cross-sectional data (1 time point per subject)',...
	'Separator',	'off',...
	'Tag',		'Segment cross-sectional data',...
	'CallBack',	'cg_vbm_optimized;',...
	'HandleVisibility','on');
h22  = uimenu(h2,...
	'Label',	'Cross-sectional data asymmetry: calculate lateralization index',...
	'Separator',	'off',...
	'Tag',		'Segment lateralization index',...
	'CallBack',	'cg_vbm_lat_index;',...
	'HandleVisibility','on');
h23  = uimenu(h2,...
	'Label',	'Longitudinal data (more than 1 time point per subject)',...
	'Separator',	'off',...
	'Tag',		'Segment longitudinal data',...
	'CallBack',	'cg_vbm_longitudinal_bias;',...
	'HandleVisibility','on');
h3  = uimenu(h0,...
	'Label',	'Statistical analysis',...
	'Separator',	'on',...
	'Tag',		'Statistical analysis',...
	'HandleVisibility','on');
h31  = uimenu(h3,...
	'Label',	'Cross-sectional data (1 time point per subject)',...
	'Separator',	'off',...
	'Tag',		'VBM analyze cross-sectional data',...
	'CallBack','SPM = cg_spm_ui(''cfg'',cg_spm_ui(''cohort_stats''));',...
	'HandleVisibility','on');
h32  = uimenu(h3,...
	'Label',	'Cross-sectional data asymmetry',...
	'Separator',	'off',...
	'Tag',		'VBM analyze cross-sectional data asymmetry',...
	'CallBack','SPM = cg_spm_ui(''cfg'',cg_spm_ui(''asymmetry_stats''));',...
	'HandleVisibility','on');
h33  = uimenu(h3,...
	'Label',	'Longitudinal data (more than 1 time point per subject)',...
	'Separator',	'off',...
	'Tag',		'VBM analyze longitudinal data',...
	'CallBack','SPM = cg_spm_ui(''cfg'',cg_spm_ui(''longitudinal_stats''));',...
	'HandleVisibility','on');
h4  = uimenu(h0,...
	'Label',	'Results of statistical analysis',...
	'Separator',	'off',...
	'Tag',		'Show results of statistical analysis',...
	'CallBack','[hReg,xSPM,SPM] = spm_results_ui;',...
	'HandleVisibility','on');
h5  = uimenu(h0,...
	'Label',	'Tools',...
	'Separator',	'on',...
	'Tag',		'Tools',...
	'HandleVisibility','on');
h51  = uimenu(h5,...
	'Label',	'Check sample homogeneity using standard deviation across sample',...
	'Separator',	'off',...
	'Tag',		'Check sample homogeneity using standard deviation across sample',...
	'CallBack',	'cg_check_sample_sd;',...
	'HandleVisibility','on');
h52  = uimenu(h5,...
	'Label',	'Display one slice for all images',...
	'Separator',	'off',...
	'Tag',		'Display one slice for all images',...
	'CallBack',	'cg_showslice_all;',...
	'HandleVisibility','on');
h53  = uimenu(h5,...
	'Label',	'Calculate raw volumes for GM/WM/CSF',...
	'Separator',	'off',...
	'Tag',		'VBM Calculate raw volumes for GM/WM/CSF',...
	'CallBack',	'[gm, wm, csf, names] = cg_read_vbm_volumes;fprintf(''The variables gm, wm, csf, and names are now saved in the workspace.\n'');',...
	'HandleVisibility','on');
h54  = uimenu(h5,...
	'Label',	'Threshold and transform spmT-maps',...
	'Separator',	'off',...
	'Tag',		'Threshold and transform spmT-maps',...
	'CallBack',	'cg_spmT2x;',...
	'HandleVisibility','on');