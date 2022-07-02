function [gm, wm, csf, vbm_names] = cg_read_vbm_volumes
%CG_READ_VBM_VOLUMES
% read raw volumes from vbm_volumes.txt
%
% If you have applied cg_vbm_optimized a file vbm_volumes.txt is saved
% containing the raw volumes for gray/white matter and CSF. You can read
% raw volumes with
% [gm wm csf] = cg_read_vbm_volumes
%
% To get sure that the order in vbm_volumes.txt is according to the files that will be used in the 
% design matrix you can order the volume values by selecting the respective files. This is helpful if 
% the order in vbm_volumes.txt is different or the number of files differs.
%
% You can use these variables either as nuisance in an AnCova model or as user-specified globals with 
% the "global calculation" option. Depending on your hypothesis and/or your data you can just use gray 
% matter ("gm") or calculate the sum of gray/white matter with "gm+wm". The use of raw volumes as 
% nuisance or globals is only recommended for modulated data. These data are corrected for size changes 
% due to spatial  normalization and are thought to be in raw (un-normalized) space. In contrast, un-
% modulated data are yet corrected for differences in size due to spatial normalization to a 
% reference brain and there is no need to correct for these differences again.
%_______________________________________________________________________
% @(#)cg_read_vbm_volumes.m	1.05 Christian Gaser 2006/06/19

Pvbm = spm_get(1,'vbm*txt','Select file');

% get volumes
out = dlmread(Pvbm,'\t',1,1);
vbm_names = [];

% get filenames
P = spm_get(Inf,'IMAGE','Select segmented files which should be included in statistics');
if ~isempty(P)
    % get filenames and check number of files
    n_files = size(P,1);
    if n_files > size(out,1)
        error('Insufficient number of files in vbm_volumes.txt');
    end
    % get filenames in vbm_volumes.txt
    name = textread(Pvbm,'%s');
    % skip first line with 5 entries and collect every 4th element
    for i=6:4:length(name)
        vbm_names = strvcat(vbm_names,name{i});
    end
    order = ones(n_files,1);
    % read last files first
    for i=n_files:-1:1
        for j=size(vbm_names,1):-1:1
            if findstr(deblank(P(i,:)),vbm_names(j,:))
                order(i) = j;
                break;
            end
        end
    end
    out = out(order,:);
    vbm_names = vbm_names(order,:);
end

gm  = out(:,1);
wm  = out(:,2);
csf = out(:,3);

