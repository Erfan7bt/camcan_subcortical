cortex=load('/data/erfan/brainstorm_db/addMix_5000/anat/020_CC120065/tess_cortex_concat_02_fix.mat');
headmodel=load('/data/erfan/brainstorm_db/addMix_5000/data/020_CC120065/@default_study/headmodel_mix_openmeeg.mat');

leftcortex=headmodel.GridAtlas.Scouts(10);
rightcortex=headmodel.GridAtlas.Scouts(11);
nROI_cortex=68;
nROI_subcortex=17;
dk=cortex.Atlas(4).Scouts;

ind_roi={};
ind_cortex=[];

for iroi=1:2:68
ind_roi{iroi}=dk(iroi).Vertices;
[~,ind_roi{iroi}]=ismember(ind_roi{iroi},leftcortex.Vertices);
ind_roi{iroi}=leftcortex.GridRows(ind_roi{iroi});
ind_cortex=cat(2,ind_cortex,ind_roi{iroi});
end


for iroi=2:2:68
ind_roi{iroi}=dk(iroi).Vertices;
[~,ind_roi{iroi}]=ismember(ind_roi{iroi},rightcortex.Vertices);
ind_roi{iroi}=rightcortex.GridRows(ind_roi{iroi});
ind_cortex=cat(2,ind_cortex,ind_roi{iroi});
end

ind_roi_sub={};
ind_subcortex=[];

%remove left and right cortex from headmodel.GridAtlas.Scouts
headmodel.GridAtlas.Scouts(10:11)=[];

for iroi=1:length(headmodel.GridAtlas.Scouts)
ind_roi_sub{iroi}=headmodel.GridAtlas.Scouts(iroi).GridRows;
ind_subcortex=cat(2,ind_subcortex,ind_roi_sub{iroi});
end

ind=cat(2,ind_cortex,ind_subcortex);

for iROI = 1:nROI_subcortex
    [~, ind_roi_subcortex{iROI}, ~] = intersect(ind, ind_roi_sub{iROI});
end

for iROI = 1:nROI_cortex

[~, ind_roi_cortex{iROI}, ~] = intersect(ind, ind_roi{iROI});
end


ind_roi_all=cat(2,ind_roi_subcortex,ind_roi_cortex);
%% add the subcortical 