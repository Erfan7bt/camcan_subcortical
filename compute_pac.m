sourceFolder_me = '/home/erfan/camcan_erfan/Results/source/';
sourceFolder_ver=[sourceFolder_me 'veronika/'];
sourceFolders={sourceFolder_ver, sourceFolder_me};
allSourceData =[];
addpath ../PAC/
addpath ../PAC/external/
addpath ../haufe/

% Create a new figure to hold all the subplots
figure('Position', [100, 100, 1200, 1200]);

for j=1:2
    sourceFolder=sourceFolders{j};
    % Get a list of all directories in the source folder
    dirList = dir(sourceFolder);
    dirList = dirList([dirList.isdir]);
    dirList = dirList(~ismember({dirList.name}, {'.', '..','veronika'}));

    % Preallocate the matrix to store all the source ROI data
    filt.low=[9 11];
    filt.high=[44 50];
    fs=100;

    % Loop through each directory
    for i = 1:numel(dirList)
        % Get the current directory path
        currentDir = fullfile(sourceFolder, dirList(i).name);

        % Load the source_roi_data field
        load(fullfile(currentDir, 'source_rec_results.mat'), 'source_roi_data');

        if  size(source_roi_data,1)<86
            source_roi_data = source_roi_data(18:end,:,1:20);
            [b_orig, b_anti, b_orig_norm,b_anti_norm] = fp_pac_bispec(source_roi_data,fs,filt);
            % Save the plot for b_anti_norm 
            subplot(2, numel(dirList), i + (j-1)*numel(dirList));
            imagesc(b_anti);
            colormap hot
%             colorbar;
            name=dirList(i).name;
            title(name);
        else 
            source_roi_data = source_roi_data(1:3:end,:,1:20);
            [b_orig, b_anti, b_orig_norm,b_anti_norm] = fp_pac_bispec(source_roi_data,fs,filt);
            subplot(2, numel(dirList), i + (j-1)*numel(dirList));
            imagesc(b_anti);
            colormap hot
%             colorbar;
            name=dirList(i).name;
            title(name);
        end

        % Save the figure
        if i == numel(dirList) && j == 2
            % Save the final plot with all figures in one file
            saveas(gcf, '/home/erfan/camcan_erfan/Results/all_figures.png');
        end
    end
end
