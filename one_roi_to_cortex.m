% Load data
load(['./Results/source/sub-CC120065/source_rec_results.mat'],'source_roi_data','labels','regions_cortex');

froi_id = 11; % Hippocampus R
cortex_id = 18:85;

fs = 100;          % Sampling frequency
freq = 1:50-1;     % Frequency range
numFreqs = length(freq);

big = zeros(length(cortex_id), numFreqs, numFreqs);

% Parallel for loop for cortex_id
for ids_idx = 1:length(cortex_id)
    ids = cortex_id(ids_idx);
    disp(ids);
    tic;
    
    % Load the data for froi_id and the current cortex id
    data = source_roi_data([froi_id, ids], :, 1:150);

    % Parallelize the frequency computations
    for low_freqs = freq
        for high_freqs = low_freqs:numFreqs
            filt.low = low_freqs;
            filt.high = high_freqs;
            
            % Call the fp_pac_bispec function
            [b_orig, b_anti, b_orig_norm, b_anti_norm] = fp_pac_bispec(data, fs, filt);
            
            % Update the 'big' matrix with the computed values
            big(ids_idx, low_freqs, high_freqs) = b_anti(1, 2);
            big(ids_idx, high_freqs, low_freqs) = b_anti(2, 1);
        end
    end
    
    toc;
end

for ids_idx = 1:length(cortex_id)
    ids = cortex_id(ids_idx);
    fig=figure('Visible','off');

    imagesc(squeeze(big(ids_idx,:,:)));
    colormap hot
    colorbar;
    name=labels{ids};
    title(name);
%     clim([0 1])
    saveas(fig,['./hipp_cort/' name '.png'])
    close all
end

% group the ids of regions_cortex that are in the same region
% and compute the mean of the bicoherence values
% for each pair of regions
regions = unique(regions_cortex);

for region_idx = 1:length(regions)
    region = regions{region_idx};
    ids = find(strcmp(regions_cortex, region));
    
    % Compute the mean of the bicoherence values
    mean_bicoherence = mean(big(ids, :, :), 1);
    
    fig=figure('Visible','off');
    imagesc(squeeze(mean_bicoherence));
    colormap hot
    colorbar;
    title(region);
%     clim([0 1])
    saveas(fig,['./hipp_cort/' region '.png'])
    close all
end

