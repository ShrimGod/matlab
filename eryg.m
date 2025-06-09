
%% Advanced Standing Wave Removal with Temporal Coherence Analysis
% This improved method identifies and removes standing waves by analyzing
% their temporal stability - the key characteristic that distinguishes them
% from propagating ultrasonic signals.
%
% Key improvements:
% 1. Temporal coherence analysis across multiple frames
% 2. Proper spatial windowing for standing wave patterns
% 3. Animation capabilities for verification
% 4. Multiple metrics to validate removal effectiveness
%
% Author: Advanced NDT Analysis Framework
% Date: 2024

clear; clc; close all;

%% ============================================================================
%% LOAD DATA AND PARAMETERS
%% ============================================================================

fprintf('=== ADVANCED STANDING WAVE REMOVAL ===\n');
fprintf('Loading ultrasonic data...\n');

% Load data
[fname, location] = uigetfile({'*.mat'; '*.*'}, 'Select Ultrasonic Data File');
if fname == 0
    error('No file selected. Exiting.');
end
load(fullfile(location, fname));
clear location fname;

% Verify data dimensions and create coordinate vectors if needed
[nx, ny, nt] = size(Data);
fprintf('Data dimensions: %d x %d x %d (x, y, time)\n', nx, ny, nt);

if ~exist('fs', 'var')
    fs = 50e6; % Default 50 MHz sampling rate
    fprintf('Using default sampling frequency: %.1f MHz\n', fs/1e6);
end

if ~exist('x', 'var') || ~exist('y', 'var') || ~exist('t', 'var')
    x = 1:nx;
    y = 1:ny;
    t = (0:nt-1) / fs * 1e6; % Time in microseconds
    fprintf('Generated default coordinate vectors\n');
end

%% ============================================================================
%% ALGORITHM PARAMETERS
%% ============================================================================

% Temporal analysis parameters
temporal_window = 20;        % Number of frames to analyze for temporal stability
overlap_ratio = 0.5;        % Overlap between temporal windows
stability_threshold = 0.85;  % Correlation threshold for standing wave detection

% Spatial analysis parameters
spatial_window = 7;         % Larger window to capture standing wave patterns
freq_range = [2e6, 8e6];   % Extended frequency range for analysis

% Filtering parameters
suppression_factor = 0.05;  % Reduce standing waves to 5% of original amplitude
transition_width = 10;      % Smooth transition to avoid artifacts

% Visualization parameters
animation_fps = 10;         % Frames per second for animations
time_range = [100, 140];    % Time range (μs) for detailed analysis

%% ============================================================================
%% TEMPORAL COHERENCE ANALYSIS
%% ============================================================================

fprintf('\n=== PERFORMING TEMPORAL COHERENCE ANALYSIS ===\n');
fprintf('This identifies regions that remain stationary over time.\n');

% Calculate the number of temporal analysis windows
n_windows = floor((nt - temporal_window) / (temporal_window * (1 - overlap_ratio))) + 1;
temporal_coherence_map = zeros(nx, ny, n_windows);

% Initialize progress tracking
tic;
total_operations = n_windows * nx * ny;
operations_done = 0;

fprintf('Analyzing %d temporal windows...\n', n_windows);

for w = 1:n_windows
    % Define the time window indices
    start_idx = round(1 + (w-1) * temporal_window * (1 - overlap_ratio));
    end_idx = min(start_idx + temporal_window - 1, nt);
    
    if end_idx > nt
        break;
    end
    
    fprintf('Window %d/%d: frames %d-%d\n', w, n_windows, start_idx, end_idx);
    
    % Analyze each spatial position
    for i = 1:nx
        for j = 1:ny
            % Extract the temporal signal at this position
            temporal_signal = squeeze(Data(i, j, start_idx:end_idx));
            
            % Calculate temporal stability using autocorrelation
            % Standing waves will show high correlation between frames
            correlations = zeros(length(temporal_signal)-1, 1);
            
            for k = 1:length(temporal_signal)-1
                % Compare consecutive frames
                r = corrcoef(temporal_signal(k), temporal_signal(k+1));
                correlations(k) = r(1,2);
            end
            
            % High mean correlation indicates temporal stability (standing wave)
            temporal_coherence_map(i, j, w) = mean(correlations);
            
            operations_done = operations_done + 1;
            if mod(operations_done, 1000) == 0
                elapsed = toc;
                progress = operations_done / total_operations * 100;
                fprintf('  Progress: %.1f%% (%.1f seconds elapsed)\n', progress, elapsed);
            end
        end
    end
end

% Combine temporal windows to create overall standing wave probability map
standing_wave_probability = mean(temporal_coherence_map, 3);
fprintf('Temporal analysis completed in %.1f seconds\n', toc);

%% ============================================================================
%% FREQUENCY-DOMAIN TEMPORAL ANALYSIS
%% ============================================================================

fprintf('\n=== FREQUENCY-DOMAIN TEMPORAL TRACKING ===\n');
fprintf('Analyzing phase evolution over time at key frequencies.\n');

% Setup frequency analysis
nfft = 2^nextpow2(nt);
freq_axis = (0:nfft-1) * fs / nfft;
freq_indices = find(freq_axis >= freq_range(1) & freq_axis <= freq_range(2));

% Analyze temporal phase evolution
phase_stability_map = zeros(nx, ny);

tic;
for i = 1:nx
    if mod(i, 10) == 0
        fprintf('Processing row %d/%d\n', i, nx);
    end
    
    for j = 1:ny
        % Extract full temporal signal
        signal = squeeze(Data(i, j, :));
        
        % Divide into segments for phase tracking
        segment_length = 50;  % Frames per segment
        n_segments = floor(nt / segment_length);
        
        phase_variations = [];
        
        for seg = 1:n_segments-1
            % Extract consecutive segments
            seg1_start = (seg-1) * segment_length + 1;
            seg1_end = seg * segment_length;
            seg2_start = seg * segment_length + 1;
            seg2_end = min((seg+1) * segment_length, nt);
            
            segment1 = signal(seg1_start:seg1_end);
            segment2 = signal(seg2_start:seg2_end);
            
            % FFT of each segment
            fft1 = fft(segment1, nfft);
            fft2 = fft(segment2, nfft);
            
            % Extract phases in frequency range
            phases1 = angle(fft1(freq_indices));
            phases2 = angle(fft2(freq_indices));
            
            % Calculate phase difference (accounting for wrapping)
            phase_diff = angle(exp(1i * (phases2 - phases1)));
            
            % Standing waves should have minimal phase change
            phase_variations = [phase_variations; std(phase_diff)];
        end
        
        % Low phase variation indicates standing wave
        phase_stability_map(i, j) = 1 / (mean(phase_variations) + 0.01);
    end
end
fprintf('Frequency-domain analysis completed in %.1f seconds\n', toc);

%% ============================================================================
%% COMBINED STANDING WAVE DETECTION
%% ============================================================================

fprintf('\n=== COMBINING DETECTION METHODS ===\n');

% Normalize both maps to [0, 1]
standing_wave_probability_norm = (standing_wave_probability - min(standing_wave_probability(:))) / ...
    (max(standing_wave_probability(:)) - min(standing_wave_probability(:)));

phase_stability_norm = (phase_stability_map - min(phase_stability_map(:))) / ...
    (max(phase_stability_map(:)) - min(phase_stability_map(:)));

% Combine both metrics with weighting
% Temporal coherence is more reliable, so we weight it higher
combined_standing_wave_map = 0.7 * standing_wave_probability_norm + 0.3 * phase_stability_norm;

% Create binary mask using adaptive thresholding
threshold = stability_threshold * max(combined_standing_wave_map(:));
standing_wave_mask = combined_standing_wave_map > threshold;

% Apply morphological operations to clean up the mask
standing_wave_mask = imclose(standing_wave_mask, ones(3, 3));
standing_wave_mask = imopen(standing_wave_mask, ones(3, 3));

fprintf('Detected standing waves in %.1f%% of spatial positions\n', ...
    100 * sum(standing_wave_mask(:)) / numel(standing_wave_mask));

%% ============================================================================
%% ADVANCED FILTERING WITH SMOOTH TRANSITIONS
%% ============================================================================

fprintf('\n=== APPLYING ADVANCED FILTERING ===\n');
fprintf('Using frequency-selective suppression with smooth transitions.\n');

Data_filtered = zeros(size(Data));

tic;
for i = 1:nx
    if mod(i, 10) == 0
        elapsed = toc;
        fprintf('Filtering progress: %d/%d rows (%.1f seconds)\n', i, nx, elapsed);
    end
    
    for j = 1:ny
        signal = squeeze(Data(i, j, :));
        
        if standing_wave_mask(i, j)
            % This position contains standing waves
            
            % Perform spectral subtraction
            fft_signal = fft(signal, nfft);
            
            % Estimate standing wave spectrum from stable periods
            % Use the most stable temporal window
            [~, most_stable_window] = max(temporal_coherence_map(i, j, :));
            stable_start = round(1 + (most_stable_window-1) * temporal_window * (1 - overlap_ratio));
            stable_end = min(stable_start + temporal_window - 1, nt);
            
            stable_segment = signal(stable_start:stable_end);
            standing_wave_spectrum = fft(stable_segment, nfft);
            
            % Subtract standing wave component with frequency-dependent weighting
            for k = 1:length(freq_indices)
                freq_idx = freq_indices(k);
                
                % Calculate suppression weight based on phase stability
                weight = 1 - suppression_factor;  % How much to keep
                
                % Apply suppression
                fft_signal(freq_idx) = fft_signal(freq_idx) - ...
                    (1 - weight) * standing_wave_spectrum(freq_idx);
                
                % Maintain conjugate symmetry
                if freq_idx > 1 && freq_idx <= nfft/2
                    fft_signal(nfft - freq_idx + 2) = conj(fft_signal(freq_idx));
                end
            end
            
            % Reconstruct signal
            filtered_signal = real(ifft(fft_signal));
            Data_filtered(i, j, :) = filtered_signal(1:nt);
        else
            % No standing waves detected - preserve original signal
            Data_filtered(i, j, :) = signal;
        end
    end
end

% Apply spatial smoothing at boundaries to prevent artifacts
fprintf('Applying boundary smoothing...\n');
for k = 1:nt
    slice = Data_filtered(:, :, k);
    mask_dilated = imdilate(standing_wave_mask, ones(transition_width));
    mask_eroded = imerode(standing_wave_mask, ones(transition_width/2));
    transition_region = mask_dilated & ~mask_eroded;
    
    if any(transition_region(:))
        % Smooth transition region
        slice_smooth = imgaussfilt(slice, 2);
        slice(transition_region) = slice_smooth(transition_region);
        Data_filtered(:, :, k) = slice;
    end
end

fprintf('Filtering completed in %.1f seconds\n', toc);

%% ============================================================================
%% VERIFICATION METRICS
%% ============================================================================

fprintf('\n=== CALCULATING VERIFICATION METRICS ===\n');

% Metric 1: Temporal Variance Reduction
temporal_variance_original = zeros(nx, ny);
temporal_variance_filtered = zeros(nx, ny);

for i = 1:nx
    for j = 1:ny
        % Calculate variance over time
        temporal_variance_original(i, j) = var(squeeze(Data(i, j, :)));
        temporal_variance_filtered(i, j) = var(squeeze(Data_filtered(i, j, :)));
    end
end

% Standing wave regions should show significant variance reduction
variance_reduction = 1 - temporal_variance_filtered ./ (temporal_variance_original + eps);
variance_reduction_sw = variance_reduction(standing_wave_mask);
variance_reduction_other = variance_reduction(~standing_wave_mask);

fprintf('\nTemporal Variance Reduction:\n');
fprintf('- In standing wave regions: %.1f%% ± %.1f%%\n', ...
    100*mean(variance_reduction_sw), 100*std(variance_reduction_sw));
fprintf('- In other regions: %.1f%% ± %.1f%%\n', ...
    100*mean(variance_reduction_other), 100*std(variance_reduction_other));

% Metric 2: Frequency Content Analysis
fprintf('\nAnalyzing frequency content changes...\n');

% Select representative positions
sw_positions = find(standing_wave_mask(:));
if length(sw_positions) > 10
    sw_positions = sw_positions(randperm(length(sw_positions), 10));
end

avg_spectrum_original = zeros(nfft/2, 1);
avg_spectrum_filtered = zeros(nfft/2, 1);

for idx = sw_positions'
    [i, j] = ind2sub([nx, ny], idx);
    
    original_spectrum = abs(fft(squeeze(Data(i, j, :)), nfft));
    filtered_spectrum = abs(fft(squeeze(Data_filtered(i, j, :)), nfft));
    
    avg_spectrum_original = avg_spectrum_original + original_spectrum(1:nfft/2);
    avg_spectrum_filtered = avg_spectrum_filtered + filtered_spectrum(1:nfft/2);
end

avg_spectrum_original = avg_spectrum_original / length(sw_positions);
avg_spectrum_filtered = avg_spectrum_filtered / length(sw_positions);

% Metric 3: Signal-to-Noise Ratio Improvement
% Define signal and noise regions based on time
[~, signal_start] = min(abs(t - time_range(1)));
[~, signal_end] = min(abs(t - time_range(2)));
noise_region = signal_end+50:min(signal_end+200, nt);

snr_improvement = zeros(nx, ny);
for i = 1:nx
    for j = 1:ny
        % Original SNR
        signal_power_orig = mean(Data(i, j, signal_start:signal_end).^2);
        noise_power_orig = mean(Data(i, j, noise_region).^2);
        snr_orig = 10*log10(signal_power_orig / (noise_power_orig + eps));
        
        % Filtered SNR
        signal_power_filt = mean(Data_filtered(i, j, signal_start:signal_end).^2);
        noise_power_filt = mean(Data_filtered(i, j, noise_region).^2);
        snr_filt = 10*log10(signal_power_filt / (noise_power_filt + eps));
        
        snr_improvement(i, j) = snr_filt - snr_orig;
    end
end

fprintf('\nSignal-to-Noise Ratio Improvement:\n');
fprintf('- Average: %.1f dB\n', mean(snr_improvement(:)));
fprintf('- In standing wave regions: %.1f dB\n', mean(snr_improvement(standing_wave_mask)));

%% ============================================================================
%% VISUALIZATION - STATIC PLOTS
%% ============================================================================

fprintf('\n=== GENERATING ANALYSIS PLOTS ===\n');

% Figure 1: Detection Results
figure('Name', 'Standing Wave Detection Results', 'Position', [50, 50, 1400, 900]);

subplot(2, 3, 1);
imagesc(x, y, standing_wave_probability.');
title('Temporal Coherence Map');
xlabel('X Position'); ylabel('Y Position');
colorbar; colormap(gca, 'hot'); axis equal tight;
cb = colorbar; ylabel(cb, 'Temporal Stability');

subplot(2, 3, 2);
imagesc(x, y, phase_stability_map.');
title('Phase Stability Map');
xlabel('X Position'); ylabel('Y Position');
colorbar; colormap(gca, 'hot'); axis equal tight;
cb = colorbar; ylabel(cb, 'Phase Stability');

subplot(2, 3, 3);
imagesc(x, y, combined_standing_wave_map.');
title('Combined Detection Map');
xlabel('X Position'); ylabel('Y Position');
colorbar; colormap(gca, 'parula'); axis equal tight;
cb = colorbar; ylabel(cb, 'Standing Wave Probability');

subplot(2, 3, 4);
imagesc(x, y, standing_wave_mask.');
title('Standing Wave Mask (Binary)');
xlabel('X Position'); ylabel('Y Position');
colormap(gca, 'gray'); axis equal tight;

subplot(2, 3, 5);
imagesc(x, y, variance_reduction.' * 100);
title('Temporal Variance Reduction (%)');
xlabel('X Position'); ylabel('Y Position');
colorbar; colormap(gca, 'jet'); axis equal tight;
cb = colorbar; ylabel(cb, 'Reduction (%)');
caxis([0, 100]);

subplot(2, 3, 6);
freq_plot = freq_axis(1:nfft/2) / 1e6;
plot(freq_plot, avg_spectrum_original, 'b-', 'LineWidth', 1.5); hold on;
plot(freq_plot, avg_spectrum_filtered, 'r-', 'LineWidth', 1.5);
xlabel('Frequency (MHz)'); ylabel('Average Amplitude');
title('Frequency Content in Standing Wave Regions');
legend('Original', 'Filtered', 'Location', 'best');
xlim([0, fs/2e6]);
grid on;

% Figure 2: Time Domain Comparison
figure('Name', 'Time Domain Analysis', 'Position', [100, 100, 1400, 600]);

% Find a representative standing wave position
[max_sw_prob, max_idx] = max(combined_standing_wave_map(:));
[sw_x, sw_y] = ind2sub([nx, ny], max_idx);

subplot(1, 2, 1);
plot(t, squeeze(Data(sw_x, sw_y, :)), 'b-', 'LineWidth', 1.5); hold on;
plot(t, squeeze(Data_filtered(sw_x, sw_y, :)), 'r-', 'LineWidth', 1.5);
xlabel('Time (μs)'); ylabel('Amplitude');
title(sprintf('A-scan at Standing Wave Position (%d, %d)', sw_x, sw_y));
legend('Original', 'Filtered', 'Location', 'best');
grid on;
xlim([t(1), t(end)]);

% Find a position without standing waves for comparison
non_sw_mask = ~standing_wave_mask;
[i_non_sw, j_non_sw] = find(non_sw_mask, 1);

subplot(1, 2, 2);
plot(t, squeeze(Data(i_non_sw, j_non_sw, :)), 'b-', 'LineWidth', 1.5); hold on;
plot(t, squeeze(Data_filtered(i_non_sw, j_non_sw, :)), 'r-', 'LineWidth', 1.5);
xlabel('Time (μs)'); ylabel('Amplitude');
title(sprintf('A-scan at Non-Standing Wave Position (%d, %d)', i_non_sw, j_non_sw));
legend('Original', 'Filtered', 'Location', 'best');
grid on;
xlim([t(1), t(end)]);

%% ============================================================================
%% ANIMATION GENERATION
%% ============================================================================

fprintf('\n=== GENERATING ANIMATIONS ===\n');
fprintf('Creating side-by-side comparison animations...\n');

% Select time range for animation
[~, anim_start] = min(abs(t - time_range(1)));
[~, anim_end] = min(abs(t - time_range(2)));
anim_frames = anim_start:2:anim_end;  % Use every other frame for speed

% Determine common color scale
all_data = [Data(:); Data_filtered(:)];
cmin = prctile(all_data, 1);
cmax = prctile(all_data, 99);

% Create animation figure
fig_anim = figure('Name', 'Standing Wave Removal Animation', ...
    'Position', [150, 150, 1200, 500]);

% Preallocate movie frames
n_anim_frames = length(anim_frames);
mov = struct('cdata', [], 'colormap', []);

fprintf('Recording %d frames...\n', n_anim_frames);
for k = 1:n_anim_frames
    frame_idx = anim_frames(k);
    
    % Original data
    subplot(1, 3, 1);
    imagesc(x, y, Data(:, :, frame_idx).');
    title(sprintf('Original Data - t = %.2f μs', t(frame_idx)));
    xlabel('X Position'); ylabel('Y Position');
    colorbar; caxis([cmin, cmax]);
    axis equal tight;
    
    % Filtered data
    subplot(1, 3, 2);
    imagesc(x, y, Data_filtered(:, :, frame_idx).');
    title(sprintf('Filtered Data - t = %.2f μs', t(frame_idx)));
    xlabel('X Position'); ylabel('Y Position');
    colorbar; caxis([cmin, cmax]);
    axis equal tight;
    
    % Difference
    subplot(1, 3, 3);
    diff_data = Data(:, :, frame_idx) - Data_filtered(:, :, frame_idx);
    imagesc(x, y, diff_data.');
    title('Removed Components (Standing Waves)');
    xlabel('X Position'); ylabel('Y Position');
    colorbar; colormap(gca, 'seismic');
    caxis([-max(abs(diff_data(:))), max(abs(diff_data(:)))]);
    axis equal tight;
    
    drawnow;
    mov(k) = getframe(fig_anim);
    
    if mod(k, 10) == 0
        fprintf('  Frame %d/%d\n', k, n_anim_frames);
    end
end

% Save animation as video file
fprintf('Saving animation as AVI file...\n');
video_writer = VideoWriter('standing_wave_removal_animation.avi');
video_writer.FrameRate = animation_fps;
open(video_writer);
writeVideo(video_writer, mov);
close(video_writer);
fprintf('Animation saved as: standing_wave_removal_animation.avi\n');

%% ============================================================================
%% WATERFALL PLOT FOR TEMPORAL EVOLUTION
%% ============================================================================

figure('Name', 'Temporal Evolution Analysis', 'Position', [200, 200, 1000, 800]);

% Select a line through a standing wave region
line_y = sw_y;  % Use the y-position of maximum standing wave

% Create waterfall plot
subplot(2, 1, 1);
waterfall_frames = 1:10:nt;  % Every 10th frame
waterfall_data_orig = zeros(length(waterfall_frames), nx);
waterfall_data_filt = zeros(length(waterfall_frames), nx);

for idx = 1:length(waterfall_frames)
    frame = waterfall_frames(idx);
    waterfall_data_orig(idx, :) = Data(:, line_y, frame);
    waterfall_data_filt(idx, :) = Data_filtered(:, line_y, frame);
end

% Original data waterfall
surf(x, t(waterfall_frames), waterfall_data_orig);
shading interp; colormap(jet);
xlabel('X Position'); ylabel('Time (μs)'); zlabel('Amplitude');
title(sprintf('Original Data Evolution at Y=%d', line_y));
view(-45, 30);

subplot(2, 1, 2);
surf(x, t(waterfall_frames), waterfall_data_filt);
shading interp; colormap(jet);
xlabel('X Position'); ylabel('Time (μs)'); zlabel('Amplitude');
title(sprintf('Filtered Data Evolution at Y=%d', line_y));
view(-45, 30);

%% ============================================================================
%% SUMMARY REPORT
%% ============================================================================

fprintf('\n=== STANDING WAVE REMOVAL SUMMARY ===\n');
fprintf('=====================================\n\n');

fprintf('Detection Results:\n');
fprintf('- Total spatial positions analyzed: %d\n', nx * ny);
fprintf('- Standing wave positions detected: %d (%.1f%%)\n', ...
    sum(standing_wave_mask(:)), 100 * sum(standing_wave_mask(:)) / numel(standing_wave_mask));
fprintf('- Maximum temporal coherence: %.3f\n', max(standing_wave_probability(:)));
fprintf('- Maximum phase stability: %.3f\n', max(phase_stability_map(:)));

fprintf('\nFiltering Performance:\n');
fprintf('- Average variance reduction in SW regions: %.1f%%\n', ...
    100 * mean(variance_reduction_sw));
fprintf('- Average variance reduction in other regions: %.1f%%\n', ...
    100 * mean(variance_reduction_other));
fprintf('- Average SNR improvement: %.1f dB\n', mean(snr_improvement(:)));
fprintf('- Maximum SNR improvement: %.1f dB\n', max(snr_improvement(:)));

fprintf('\nFrequency Analysis:\n');
[~, peak_freq_idx] = max(avg_spectrum_original);
peak_freq = freq_axis(peak_freq_idx) / 1e6;
fprintf('- Dominant frequency in SW regions: %.2f MHz\n', peak_freq);
reduction_at_peak = 100 * (1 - avg_spectrum_filtered(peak_freq_idx) / avg_spectrum_original(peak_freq_idx));
fprintf('- Amplitude reduction at dominant frequency: %.1f%%\n', reduction_at_peak);

fprintf('\nOutput Files Generated:\n');
fprintf('- standing_wave_removal_animation.avi\n');
fprintf('- Multiple analysis figures\n');

fprintf('\nRecommendations:\n');
if mean(variance_reduction_sw) > 0.5
    fprintf('✓ Excellent standing wave suppression achieved\n');
else
    fprintf('⚠ Moderate suppression - consider adjusting parameters\n');
end

if mean(variance_reduction_other) < 0.1
    fprintf('✓ Good preservation of non-standing wave signals\n');
else
    fprintf('⚠ Some signal loss in non-SW regions - reduce suppression factor\n');
end

% Save results
fprintf('\nSaving analysis results...\n');
save('standing_wave_removal_results.mat', ...
    'Data_filtered', 'standing_wave_mask', 'combined_standing_wave_map', ...
    'temporal_coherence_map', 'phase_stability_map', 'variance_reduction', ...
    'snr_improvement', '-v7.3');

fprintf('\nAnalysis complete! Results saved to standing_wave_removal_results.mat\n');
fprintf('Animation saved to standing_wave_removal_animation.avi\n');
