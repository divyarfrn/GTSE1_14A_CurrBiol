function compare_cell_lines(types, use_patch, bin_size, astral_or_comet, colors, ylim_counts_max, ...
                legend_index, dirname, num_cells)
%COMPARE_CELL_LINES Loads all cell samples for type1, type2 and compares

% bin_size: one of [0.2, 0.5, 1, 2]
% astral_or_comet: one of 'astral' or 'comet'
%

% Process inputs
% restrict the number of cells to be same across all conditions per experiment
legend_text = {};
legend_text{end+1} = 'U2OS WT - Ctrl RNAi';
legend_text{end+1} = 'U2OS WT - GTSE1 RNAi';
legend_text{end+1} = 'U2OS WT - Kif18B RNAi';
legend_text{end+1} = 'U2OS WT - GTSE1+Kif18B RNAi';
legend_text{end+1} = 'U2OS^{(323 WT.13)} - GTSE1 RNAi';
legend_text{end+1} = 'U2OS^{(323 WT.13)} - GTSE1+Kif18B RNAi';
legend_text{end+1} = 'U2OS^{(323 14A.7)} - GTSE1 RNAi';
legend_text{end+1} = 'U2OS^{(323 14A.7)} - GTSE1+Kif18B RNAi';
legend_text = legend_text(legend_index);

% Setup defaults
if nargin < 2, use_patch = 'percentile'; end
if nargin < 3, bin_size = 0.5; end
if nargin < 4, astral_or_comet = 'astral'; end
if nargin < 5, colors = 'rbgcmyk'; end
if nargin < 6, ylim_counts_max = 4000; end

% Validate
if ~any(bin_size == [0.2, 0.5, 1, 2]), error('Using invalid bin-size'); end
if ~any(strcmp({'astral', 'comet'}, astral_or_comet)), error('Invalid astral_or_comet selection'); end

% Collect astral and comet lengths for all types
fprintf('Now comparing:\n');
cell_lengths = cell(1, length(types));
mat_path = cell(1, length(types));
for k = 1:length(types)
    [cell_lengths{k}, mat_path{k}] = collect_lengths(['../data/' dirname], num_cells, types{k}, astral_or_comet);
end

% Bin it!
MAX_DIST = 20 + bin_size;  % add the last 0.5 to truncate last accumulative bin [16, Inf]
bin_centers = linspace(0 + bin_size/2, MAX_DIST-bin_size/2, MAX_DIST/bin_size);

distrib = cell(1, length(types));
avg_distrib = cell(1, length(types));
for k = 1:length(types)
    num_cells = length(cell_lengths{k});
    % Compute distributions for individual cells
    for c = 1:num_cells
        distrib{k}{c} = hist(cell_lengths{k}{c}, bin_centers);
        % truncate last bin [MAX_DIST, Inf]
        % distrib{k}{c} = distrib{k}{c}(1:end-1);
    end
    % Compute stacked distribution
    all_lengths = cat(1, cell_lengths{k}{:});
    avg_distrib{k} = hist(all_lengths, bin_centers);
    % Average if plotting along with patches/individual lines
    if ~strcmp(use_patch, '')
        avg_distrib{k} = avg_distrib{k}/num_cells;
    end
    % truncate last bin [MAX_DIST, Inf]
    % avg_distrib{k} = avg_distrib{k}(1:end-1);
end
% bin_centers = bin_centers(1:end-1);

% Plot all
figure; hold on;
for k = 1:length(types)
    % plot the main line
    line_style = '-';
    if legend_index(k) > 7, line_style = '--'; end
    avg_h{k} = plot(bin_centers, avg_distrib{k}, 'Color', colors(k), 'LineWidth', 2, 'LineStyle', line_style);

    if any(strcmp({'minmax', 'percentile'}, use_patch))
        % plot patch around it :)
        edge_color = avg_h{k}.Color + (1 - avg_h{k}.Color)*0.8;
        each_line = cat(1, distrib{k}{:});
        patch_x = [bin_centers, fliplr(bin_centers)];
        if strcmp(use_patch, 'minmax')
            patch_y = [min(each_line), fliplr(max(each_line))];  % patch y-axis min-max
        elseif strcmp(use_patch, 'percentile')
            patch_y = [prctile(each_line, 25), fliplr(prctile(each_line, 75))];  % patch y-axis 10-90 percentile
        end
        patch(patch_x, patch_y, 1, ...
              'FaceColor', colors(k), 'FaceAlpha', 0.1, ...
              'EdgeColor', edge_color);
    
    % plot the actual lines
    elseif strcmp(use_patch, 'lines')
        for c = 1:length(distrib{k})
            h = plot(bin_centers, distrib{k}{c}, 'Color', colors(k), 'LineWidth', 0.5);
            h.set('Color', [h.Color, 0.2])
        end
        
    end
end
xlim([0, 20]);
% title(astral_or_comet);
legend([avg_h{:}], legend_text);
xlabel('Lengths (\mum)');
if strcmp(use_patch, ''), ylabel('Counts'); ylim([0, ylim_counts_max]);
else ylabel('Distribution'); ylim([0, 70]);
end
set(gca, 'FontSize', 14);
set(gcf, 'Color', 'w', 'Position', [100, 100, 800, 600]);
grid on;

%% save figure
H = gcf;
saveas(H, fullfile(['..\figures\' dirname '\'], sprintf('%d.fig', H.Number)));
saveas(H, fullfile(['..\figures\' dirname '\'], sprintf('%d.png', H.Number)));
% saveas(H, fullfile(['figures\20191013_Metaphase\'], sprintf('%d.fig', H.Number)));
% saveas(H, fullfile(['figures\20191013_Metaphase\'], sprintf('%d.png', H.Number)));

end

function [cell_lengths, mat_path] = collect_lengths(dirname, num_cells, cell_type, astral_or_comet)

% List all mat files for cells of a specific type
mat_files = dir([fullfile(dirname, cell_type), '.mat']);
fprintf('%s: Found %d cells, ', cell_type, length(mat_files));

% Split type-info to get root directory
type_info = regexp(cell_type, '/', 'split');

% Collect all lengths
cell_lengths = {};
mat_path = {};
for k = 1:length(mat_files)
    % get a fixed number of cells for each condition
    if k > num_cells, continue; end
    % get mat full path
    mat_path{k} = fullfile(dirname, type_info{1}, mat_files(k).name);
    data = load(mat_path{k});
    if strcmp(astral_or_comet, 'astral')
        cell_lengths{k} = data.astral_lengths;
    elseif strcmp(astral_or_comet, 'comet')
        cell_lengths{k} = data.comet_lengths;
    end
end

fprintf('Used: %d cells\n', length(cell_lengths));

end
