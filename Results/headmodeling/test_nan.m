% Define the root directory
rootDir = '/home/erfan/camcan_erfan/Results/headmodeling';

% Get a list of all subdirectories
subDirs = dir(rootDir);
subDirs = subDirs([subDirs.isdir]);
subDirs = subDirs(~ismember({subDirs.name}, {'.', '..'}));
L={};
% Loop through each subdirectory
for i = 1:numel(subDirs)
    subDirPath = fullfile(rootDir, subDirs(i).name);
    
    % Load the bs_results.mat file
    bsResultsFile = fullfile(subDirPath, 'bs_results.mat');
    load(bsResultsFile, 'leadfield');
    disp(size(leadfield));
    L{i}=leadfield;
    % Check the number of NaN values in the leadfield
    numNaN = sum(isnan(leadfield(:)));
    
    % Display the result
    fprintf('Subdirectory: %s\n', subDirPath);
    fprintf('Percentage of NaN values: %.2f%%\n', numNaN / numel(leadfield) * 100);
    
end
% Calculate the correlation matrix
disp(size(L));