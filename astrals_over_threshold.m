function astrals_over_threshold(ROOT_DIR, exp_id, spt)
close all;
clc;

save_dir = 'figures\AstralsOverThreshold\';

cell_lines = dir(ROOT_DIR);
cell_lines(1:2) = [];

thresholds = [7];

% save gtse counts / ctrl counts / experiment name
ctrl_counts = cell(length(cell_lines), length(thresholds));
gtse_counts = cell(length(cell_lines), length(thresholds));
cell_line_name = cell(1, length(cell_lines));

for l = 1:length(cell_lines)
    cell_line_name{l} = cell_lines(l).name;
    fprintf('Processing %s\n', cell_lines(l).name);
    cl_dir = [ROOT_DIR cell_lines(l).name '/']; 

    ctrl_mat_files = dir([cl_dir, 'Ctrl*.mat']);
    ctrl_al = read_astral_lengths(ctrl_mat_files, cl_dir, 'comet');

    gtse_mat_files = dir([cl_dir, 'GTSE*.mat']);
    gtse_al = read_astral_lengths(gtse_mat_files, cl_dir, 'comet');
    % compute counts by thresholding
    for t = 1:length(thresholds)
        ctrl_counts{l, t} = get_astral_counts(ctrl_al, thresholds(t));
        gtse_counts{l, t} = get_astral_counts(gtse_al, thresholds(t));
        % prepare to save
        fig_fname = sprintf('%s_%s_%d', exp_id, cell_lines(l).name, thresholds(t));
        % make plot
        dotplot(ctrl_counts{l, t}, gtse_counts{l, t}, save_dir, fig_fname, thresholds(t));
        % save as csv
        csv_fname = sprintf('%s_%s_%d.csv', exp_id, cell_lines(l).name, thresholds(t));
        counts = nan(max(length(ctrl_counts{l, t}), length(gtse_counts{l, t})), 2);
        counts(1:length(ctrl_counts{l, t}), 1) = ctrl_counts{l, t};
        counts(1:length(gtse_counts{l, t}), 2) = gtse_counts{l, t};
        dlmwrite([save_dir csv_fname], counts);
    end
end

%% plot 4 things together
for t = 1:length(thresholds)
    fig_fname = sprintf('%s_dotplot4_%d', exp_id, thresholds(t));
    % get the 4 lines
    u2os_c = ctrl_counts{strcmp(cell_line_name, 'U2OS'), t};
    u2os_g = gtse_counts{strcmp(cell_line_name, 'U2OS'), t};
    wt13_g = ctrl_counts{strcmp(cell_line_name, 'WT.13'), t};
    a147_g = ctrl_counts{strcmp(cell_line_name, '14A.07'), t};
    % plot and save
    dotplot4things(u2os_c, u2os_g, wt13_g, a147_g, save_dir, fig_fname, thresholds(t))
end

end


function dotplot4things(u2os_c, u2os_g, wt13_g, a147_g, save_dir, fname, t)
% u2os_c = U2OS CTRL
% u2os_g = U2OS GTSE
% wt13_g = WT.13 GTSE
% a417_g = 14A.7 GTSE

figure; hold on;
COUNTS = length(u2os_c);

% draw dots
scatter(0 + ones(COUNTS, 1), u2os_c, 20, 'k', 'filled');
scatter(1 + ones(COUNTS, 1), u2os_g, 20, 'k', 'filled');
scatter(2 + ones(COUNTS, 1), wt13_g, 20, 'k', 'filled');
scatter(3 + ones(COUNTS, 1), a147_g, 20, 'k', 'filled');

% draw average lines
line([0.8, 1.2], [mean(u2os_c), mean(u2os_c)], 'LineWidth', 1, 'Color', 'r');
line([1.8, 2.2], [mean(u2os_g), mean(u2os_g)], 'LineWidth', 1, 'Color', 'r');
line([2.8, 3.2], [mean(wt13_g), mean(wt13_g)], 'LineWidth', 1, 'Color', 'r');
line([3.8, 4.2], [mean(a147_g), mean(a147_g)], 'LineWidth', 1, 'Color', 'r');

% figure properties
set(gca, 'XTick', [1, 2, 3, 4], 'XTickLabel', ...
    {'U2OS Ctrl', 'U2OS GTSE', 'WT.13 GTSE', '14A.7 GTSE'});
grid on;
xlim([0.5, 4.5]);
ylim([0, 2000]);
title(strrep(fname, '_', ' '));

% save figure
H = gcf;
set(H, 'Color', 'w');  % 'Position', [100, 100, 600, 600]);
saveas(H, [save_dir fname '.fig']);
saveas(H, [save_dir fname '.png']);

end


function dotplot(left_counts, right_counts, save_dir, fname, t)
figure; hold on;
COUNTS = length(left_counts);

% draw dots
xsc = ones(COUNTS, 1);
% xsc = create_xsc(nodox_counts);
scatter(xsc, left_counts, 20, 'k', 'filled');
xsc = 2 * ones(COUNTS, 1);
% xsc = 1 + create_xsc(dox_counts);
scatter(xsc, right_counts,   20, 'k', 'filled');
% draw average lines
line([0.8, 1.2], [mean(left_counts), mean(left_counts)], 'LineWidth', 1, 'Color', 'r');
line([1.8, 2.2], [mean(right_counts), mean(right_counts)], 'LineWidth', 1, 'Color', 'r');

% figure properties
% set(gca, 'XTick', [1, 2], 'XTickLabel', {'nodox', 'dox'});
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Ctrl', 'GTSE'});
grid on;
xlim([0.5, 2.5]);
if t < 7, ylim([0, 4000]);
else ylim([0, 2000]);
end
% ylim([0, 100*(12 - t)]);
title(strrep(fname, '_', ' '));

% save figure
H = gcf;
set(H, 'Color', 'w');  % 'Position', [100, 100, 600, 600]);
saveas(H, [save_dir fname '.fig']);
saveas(H, [save_dir fname '.png']);

end

function xsc = create_xsc(counts)
uniq_y = unique(counts);
xsc = ones(length(counts), 1);
if length(uniq_y) == 8 && uniq_y(end) == 14
    uniq_y = [[0, 0]; [1, 2]; [3, 4]; [5, 6]; [14, 14]];
elseif length(uniq_y) == 16 && uniq_y(end) == 42
    uniq_y = [[0, 0]; [1, 2]; [3, 4]; [5, 6]; [8, 9]; [10, 10]; [12, 12];
              [14, 15]; [18, 18]; [29, 29]; [42, 42]];
else
    keyboard;
end

for k = 1:size(uniq_y, 1)
    points_at_y = counts >= uniq_y(k, 1) & counts <= uniq_y(k, 2);
    npoints = sum(points_at_y);
    if npoints == 1
        continue
    elseif npoints <= 3
        xsc(points_at_y) = linspace(0.9, 1.1, npoints);
    else
        xsc(points_at_y) = linspace(0.8, 1.2, npoints);
    end
end
end

function al = read_astral_lengths(fnames, cl_dir, astral_or_comet)

if ~exist('astral_or_comet', 'var')
    astral_or_comet = 'astral';
end

al = cell(1, length(fnames));
for f = 1:length(fnames)
    if strcmp(astral_or_comet, 'astral')
        load([cl_dir fnames(f).name], 'astral_lengths');
        al{f} = astral_lengths;
    elseif strcmp(astral_or_comet, 'comet')
        load([cl_dir fnames(f).name], 'comet_lengths');
        al{f} = comet_lengths;
    else
        error('Invalid argument for astral_or_comet');
    end
end

end

function counts = get_astral_counts(al, t)

counts = zeros(length(al), 1);
for k = 1:length(al)
    counts(k) = sum(al{k} >= t);
end

end
