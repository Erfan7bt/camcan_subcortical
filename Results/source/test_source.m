% Define the source folder path
sourceFolder_me = '/home/erfan/camcan_erfan/Results/source/';
sourceFolder_ver=[sourceFolder_me 'veronika/'];
sourceFolders={sourceFolder_ver, sourceFolder_me};
allSourceData =[];

for j=1:2
    sourceFolder=sourceFolders{j};
% Get a list of all directories in the source folder
dirList = dir(sourceFolder);
dirList = dirList([dirList.isdir]);
dirList = dirList(~ismember({dirList.name}, {'.', '..','veronika'}));


% Preallocate the matrix to store all the source ROI data


% Loop through each directory
for i = 1:numel(dirList)
    % Get the current directory path
    currentDir = fullfile(sourceFolder, dirList(i).name);
    
    % Load the source_roi_data field
    load(fullfile(currentDir, 'source_rec_results.mat'), 'source_roi_data');

  
    if  size(source_roi_data,1)<86
        source_roi_data = source_roi_data(18:end,:,1:20);
        %make the data 2d by reshaping
        source_roi_data = reshape(source_roi_data, 1,  size(source_roi_data,1) *size(source_roi_data,2)*size(source_roi_data,3));
        %print the source_roi_data size
        source_roi_data=abs(source_roi_data);
        disp(size(source_roi_data));
    
        % Append the current source ROI data to the matrix
        allSourceData(:, i) = source_roi_data;
    else 
        source_roi_data = source_roi_data(1:3:end,:,1:20);
        %make the data 2d by reshaping
        source_roi_data = reshape(source_roi_data, 1,  size(source_roi_data,1) *size(source_roi_data,2)*size(source_roi_data,3));
        %print the source_roi_data size
          source_roi_data=abs(source_roi_data);
        disp(size(source_roi_data));
    
        % Append the current source ROI data to the matrix
        allSourceData(:, i+9) = source_roi_data;
    end
    
end
end
% Calculate the covariance matrix
covarianceMatrix = corrcoef(allSourceData)
covarianceMatrix (10:18,1:9)
% Plot the covariance matrix
figure;
imagesc(covarianceMatrix);
colormap hot
colorbar;
title('Covariance Matrix - All Source ROI Data');

% Save the plot
saveas(gcf, fullfile(sourceFolder, 'covariance_plot.png'));
