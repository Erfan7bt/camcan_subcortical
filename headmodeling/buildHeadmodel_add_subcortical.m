[ret, name] = system('hostname');

if startsWith(name,'ra')

    addpath '/data/erfan/brainstorm3'

    BrainstormDbDir = '/data/erfan/brainstorm_db/';
    
    camcan_dir = '/home/erfan/camcan/';
else 
    
    BrainstormDbDir = '/home/erfan/Thesis/brainstorm_db/';
    
    camcan_dir = '/home/erfan/Thesis/CamCan_erfan/';

    addpath /home/erfan/Thesis/brainstorm3/
end 

fig_folder = [camcan_dir '/Results/Figures/headmodeling/'];

results_folder = [camcan_dir 'Results/headmodeling/']; mkdir(results_folder)
thick=4;

asegVertices_1=500;
asegVertices_2=1250;

protocol_name='addMix_5000';  

brainstorm  

bst_set('BrainstormDbDir', BrainstormDbDir)
    
db_reload_database()

iProtocol= bst_get('Protocol',protocol_name);

bst_set('iProtocol', iProtocol)

proto_db = bst_get('ProtocolSubjects');

for isbj=length(proto_db.Subject)
 tic
    sbj= proto_db.Subject(isbj).Name;
    mkdir([fig_folder sbj])
%%
    for isubj = 1:length(proto_db.Subject)
    if isequal(proto_db.Subject(isubj).Name, sbj)
      isbj_db = isubj;  % should be identical to isbj
    end
    end
    
    sbj_db = bst_get('Subject', isbj_db);
    
    MriFile = [sbj_db.Anatomy(sbj_db.iAnatomy).FileName];
    subCortexKeeps_All   = {'Accumbens L','Accumbens R','Amygdala L', 'Amygdala R', ...
        'Brainstem','Caudate L','Caudate R','Cerebellum L','Cerebellum R', 'Pallidum L',...
        'Pallidum R','Hippocampus L', 'Hippocampus R','Putamen L','Putamen R', 'Thalamus L', 'Thalamus R'};
%% Remove Existing BEM shells 
    
    for isurf = 1:length(sbj_db.Surface)
    if startsWith(sbj_db.Surface(isurf).Comment, 'bem') 
      isDeleted = file_delete([BrainstormDbDir '/addMix/anat/' sbj_db.Surface(isurf).FileName], 1,-1 );
    end
    end
    
    db_reload_database()
    
    sbj_db = bst_get('Subject', isbj_db);
    %% Less Vertices for the subcortical aseg Atlas
    newAsegFile_high=[sbj  '/tess_aseg.mat'];


for isurf = 1:length(sbj_db.Surface)
if strfind(sbj_db.Surface(isurf).FileName, 'tess_aseg.mat')     
  [newAsegFile2, iSurface2, I, J] = tess_downsize(file_fullpath(sbj_db.Surface(isurf).FileName), asegVertices_2, 'reducepatch');
end
end

 for isurf = 1:length(sbj_db.Surface)
if strfind(sbj_db.Surface(isurf).FileName, 'tess_aseg.mat')     
  [newAsegFile, iSurface, I, J] = tess_downsize(file_fullpath(sbj_db.Surface(isurf).FileName), asegVertices_1, 'reducepatch');
end
end
    %%
    panel_scout('SetCurrentSurface', newAsegFile);
    sScouts = panel_scout('GetScouts');
    [~, iScouts] = ismember(subCortexKeeps_All, {sScouts.Label});
    panel_scout('SetSelectedScouts', iScouts);
    newAsegFile = panel_scout('NewSurface', 1);
    %% Merge the low aseg and low cortex 
    for isurf = 1:length(sbj_db.Surface)
    if strfind(sbj_db.Surface(isurf).FileName, 'tess_cortex_mid_low_2000V.mat')
       defaultCortex=sbj_db.Surface(isurf).FileName;    
    end
 end

for isurf = 1:length(sbj_db.Surface)
    if strfind(sbj_db.Surface(isurf).FileName, 'tess_cortex_mid_high.mat')
       defaultCortex_high=sbj_db.Surface(isurf).FileName;    
    end
end

for isurf = 1:length(sbj_db.Surface)
if strfind(sbj_db.Surface(isurf).FileName, 'tess_cortex_mid_low.mat')
   defaultCortex2=sbj_db.Surface(isurf).FileName;    
end
end

TessFiles={newAsegFile, defaultCortex};
TessFiles2={newAsegFile2, defaultCortex2};
TessFiles_high={newAsegFile_high, defaultCortex_high};

[mixedFile2,iSurface]=tess_concatenate(TessFiles,'Cortex_mix_2000_500','Cortex');
[mixedFile2,iSurface]=tess_concatenate(TessFiles2,'Cortex_mix_5000_1250','Cortex');
[mixedFile_high,iSurface]=tess_concatenate(TessFiles_high,'Cortex_mix_high','Cortex');
    %%
    % Atlas with structures
    atlasName = 'Structures';
    % Display mixed cortex
    hFigMix = view_surface(mixedFile2);
    [~, sSurf] = panel_scout('GetScouts');
    iAtlas = find(strcmpi(atlasName, {sSurf.Atlas.Name}));
    panel_scout('SetCurrentAtlas', iAtlas, 1);
    panel_surface('SelectHemispheres', 'struct');
    out_figure_image(hFigMix, [fig_folder sbj '/Mixed_cortex.tif'])
    pause(1);
    % Unload everything
    bst_memory('UnloadAll', 'Forced');
    %% ===== LOCATIONS AND ORIENTATIONS CONSTRAINTS =====
    % Select atlas with structures
    panel_scout('SetCurrentSurface', mixedFile2);
    [~, sSurf] = panel_scout('GetScouts');
    iAtlas = find(strcmpi(atlasName, {sSurf.Atlas.Name}));
    panel_scout('SetCurrentAtlas', iAtlas, 1);
    % Create source model atlas
    panel_scout('CreateAtlasInverse');
    % Set modeling options
    sScouts = panel_scout('GetScouts');
    % Set location and orientation constraints
    for iScout = 1 : length(sScouts)
        panel_scout('SetLocationConstraint', iScout, 'Deep brain');
        % Set orientation constraint
        panel_scout('SetOrientationConstraint', iScout, '');
    end
    bst_memory('UnloadAll', 'Forced');
    %%
    db_reload_database()
    sbj_db = bst_get('Subject', isbj_db);
 for isurf = 1:length(sbj_db.Surface)
    if strfind(sbj_db.Surface(isurf).FileName, 'tess_concat.mat')
      db_surface_default(isbj_db, 'Cortex', isurf);      
    end
  end
    %% create new BEM shells
    sFiles = [];
    
    % Process: Generate BEM surfaces
    sFiles = bst_process('CallProcess', 'process_generate_bem', sFiles, [], ...
        'subjectname', sbj, ...
        'nscalp',      1922, ...
        'nouter',      1922, ...
        'ninner',      1922, ...
        'thickness',   thick, ...
        'method',      'brainstorm');  % Brainstorm
     %% PLOT BEM
     SurfaceFile = [BrainstormDbDir protocol_name '/anat/' sbj '/tess_cortex_mid_low_2000V.mat'];
     SurfaceFile2 = [BrainstormDbDir protocol_name '/anat/' sbj '/tess_innerskull_bem_1922V.mat'];
     SurfaceFile3 = [BrainstormDbDir protocol_name '/anat/' sbj '/tess_outerskull_bem_1922V.mat'];
     SurfaceFile4 = [BrainstormDbDir protocol_name '/anat/' sbj '/tess_head_bem_1922V.mat'];
     [hFig, iDS, iFig] = view_surface(SurfaceFile,0,[1 .6 0],[],0);

     [hFig, iDS, iFig] = view_surface(SurfaceFile2,0.1);
     [hFig, iDS, iFig] = view_surface(SurfaceFile3);
     [hFig, iDS, iFig] = view_surface(SurfaceFile4);
     out_figure_image(hFig, [fig_folder sbj '/BEM_surfaces_1.tif']);
    
     % View from different angles:
     az = 0; el = 0; view(az, el)
     out_figure_image(hFig, [fig_folder sbj '/BEM_surfaces_2.tif']);
    
     az = -90; el = 0; view(az, el)
     out_figure_image(hFig, [fig_folder sbj '/BEM_surfaces_3.tif']);
    
     az = -180; el = 0; view(az, el)
     out_figure_image(hFig, [fig_folder sbj '/BEM_surfaces_4.tif']);
    
     az = -270; el = 0; view(az, el)
     out_figure_image(hFig, [fig_folder sbj '/BEM_surfaces_5.tif']);
    
     bst_memory('UnloadAll', 'Forced')% Close all the existing figures.
    
     %% PLOT BEM IN MRI:
    
     [hFig, iDS, iFig] = view_mri(MriFile, SurfaceFile)
     out_figure_image(hFig, [fig_folder sbj '/mri_with_cortex.tif']);
     bst_memory('UnloadAll', 'Forced')
    
    
     [hFig, iDS, iFig] = view_mri(MriFile, SurfaceFile2)
     out_figure_image(hFig, [fig_folder sbj '/mri_with_innerskull.tif']);
     bst_memory('UnloadAll', 'Forced')
    
     [hFig, iDS, iFig] = view_mri(MriFile, SurfaceFile3)
     out_figure_image(hFig, [fig_folder sbj '/mri_with_outerskull.tif']);
     bst_memory('UnloadAll', 'Forced')
    
     [hFig, iDS, iFig] = view_mri(MriFile, SurfaceFile4)
     out_figure_image(hFig, [fig_folder sbj '/mri_with_head.tif']);
     bst_memory('UnloadAll', 'Forced')
     
    %% Force in the vertices out of inner skull if needed
    db_reload_database()
    sbj_db = bst_get('Subject', isbj_db);
    bst_progress('start', 'Fix cortex surface', 'Loading surfaces...');
    % Load surface file
    TessMat = in_tess_bst(mixedFile2);
    TessMat.Faces    = double(TessMat.Faces);
    TessMat.Vertices = double(TessMat.Vertices);
    % Load envelope file
    EnvMat = in_tess_bst([sbj '/tess_innerskull_bem_1922V.mat']);
    EnvMat.Faces    = double(EnvMat.Faces);
    EnvMat.Vertices = double(EnvMat.Vertices);
    % Compute best fitting sphere from envelope
    bfs_center = bst_bfs(EnvMat.Vertices);
    
    % Center the two surfaces on the center of the sphere
    vCortex = bst_bsxfun(@minus, TessMat.Vertices, bfs_center(:)');
    vInner = bst_bsxfun(@minus, EnvMat.Vertices, bfs_center(:)');
    % Convert to spherical coordinates
    [thCortex, phiCortex, rCortex] = cart2sph(vCortex(:,1), vCortex(:,2), vCortex(:,3));
    % Look for points of the cortex inside the innerskull
    iVertOut = find(~inpolyhd(vCortex, vInner, EnvMat.Faces));
    % If no points outside, nothing to do
    if ~isempty(iVertOut)  
        disp(['forcing the cortex in the skull for the subject' sbj])
         for isurf = 1:length(sbj_db.Surface)
            if strfind(sbj_db.Surface(isurf).FileName, 'tess_innerskull_bem_1922V.mat')
               [mixedFile_fixed, iSurface] = tess_force_envelope(mixedFile2, sbj_db.Surface(isurf).FileName);
        
            end
         end

           bst_memory('UnloadAll', 'Forced')
   
        db_reload_database()
        sbj_db = bst_get('Subject', isbj_db);
        db_surface_type(mixedFile_fixed, 'Cortex');
        db_reload_database()
    end
    %%
    ChannelFile = [BrainstormDbDir protocol_name '/data/' sbj '/@default_study/channel_vectorview306_acc1.mat'];
    HeadFile = [BrainstormDbDir protocol_name '/anat/' sbj '/tess_head_mask.mat'];

hFig = view_surface(HeadFile);
[hFig, iDS, iFig] = view_helmet(ChannelFile, hFig)
hFig = view_channels(ChannelFile, 'MEG', 1, 0, hFig);
out_figure_image(hFig, [fig_folder sbj '/sensors_on_orig_head_1.tif']);
 
% View from different angles:
az = 0; el = 0; view(az, el)
out_figure_image(hFig, [fig_folder sbj '/sensors_on_orig_head_2.tif']);

az = -90; el = 0; view(az, el)
out_figure_image(hFig, [fig_folder sbj '/sensors_on_orig_head_3.tif']);

az = -180; el = 0; view(az, el)
out_figure_image(hFig, [fig_folder sbj '/sensors_on_orig_head_4.tif']);

az = -270; el = 0; view(az, el)
out_figure_image(hFig, [fig_folder sbj '/sensors_on_orig_head_5.tif']);

bst_memory('UnloadAll', 'Forced')% Close all the existing figures.
    %%
    sRecFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname',   sbj, ...
        'condition',     '', ...
        'tag',           '', ...
        'includebad',    0, ...
        'includeintra',  0, ...
        'includecommon', 0);
    
    sFiles = bst_process('CallProcess', 'process_headmodel', sRecFiles, [], ...
        'Comment',     '', ...
        'sourcespace', 3, ...  % Custom source model
        'meg',         4, ...  % OpenMEEG BEM
        'eeg',         3, ...  % OpenMEEG BEM
        'ecog',        2, ...  % OpenMEEG BEM
        'seeg',        2, ...  % OpenMEEG BEM
        'openmeeg',    struct(...
             'BemSelect',    [1, 1, 1], ...
             'BemCond',      [1, 0.0125, 1], ...
             'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
             'BemFiles',     {{}}, ...
             'isAdjoint',    0, ...
             'isAdaptative', 1, ...
             'isSplit',      0, ...
             'SplitLength',  []), ...
        'channelfile', '');
    toc
end