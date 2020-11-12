function line_scan_cortical(thr, filtmethod)
%LINE_SCAN_CORTICAL
%

% for t = [30, 50, 70]; for m = {'none', 'avg'}; line_scan_cortical(t, m{1}); end; end

% single experiment thresholds
if ~exist('thr', 'var'), thr = 30; end
if ~exist('filtmethod', 'var'), filtmethod = 'none'; end

clc;
close all;

ROOT_DIR = '../data/Linescans/';
FIG_DIR = '../figures/Linescans/';
if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end

DataStruct = struct('exp', '', 'cline', '', 'cnum', 0, 'xend', nan, ...
                    'bg', 0, 'peaks', nan);
exps = dir(ROOT_DIR);
exps(1:2) = [];
exps = only_folders(exps);
e_thr = [30, 30, 30];  %  -- new thresholds for each experiment
filtmethod = 'none';
for e = 1:length(exps)
    thr = e_thr(e);
    exp_dir = [ROOT_DIR exps(e).name '/'];
    %%% background information
    bg_fname = [exp_dir 'background_' exps(e).name '.txt'];
    if ~exist(bg_fname, 'file')
        fprintf('BG file not found: %s', bg_fname);
        keyboard;
    end
    bg = read_background(bg_fname);

    %%% go over cell lines
    clines = dir(exp_dir);
    clines(1:2) = [];
    clines = only_folders(clines);

    for c = 1:length(clines)
        cell_dir = [exp_dir clines(c).name '/'];
        cells = dir([cell_dir, '*.csv']);
        fprintf('%s - %s - %d cells\n', exps(e).name, clines(c).name, length(cells));
        
        %%% process each cell
        for k = 1:length(cells)
%         for k = 3:5
            fprintf('\tProcessing %s\n', cells(k).name);
            cell_fname = [cell_dir cells(k).name];
            [x, y] = read_one_linescan(cell_fname);
            % get cell number
            cell_number = str2double(cells(k).name(end-5:end-4));
            % pick appropriate bg
            row = strcmp(bg(:, 1), clines(c).name) .* ...
                  cell2mat(bg(:, 2)) == cell_number;
            if sum(row) == 0, fprintf(2, '\b, skipped, not found\n');continue; end
            bg_value = bg{row, 3};
            
            % save information to structure
            DataStruct(end+1).exp = exps(e).name;
            DataStruct(end).cline = clines(c).name;
            DataStruct(end).cnum = cell_number;
            DataStruct(end).bg = bg_value;
            DataStruct(end).peaks = nan;
            DataStruct(end).xend = x(end);

            % get peaks
            y = y - bg_value;
            if strcmp(filtmethod, 'none')
                npeak = analyze_linescan(x, y, thr);
                DataStruct(end).peaks = npeak;
            elseif strcmp(filtmethod, 'avg')
                npeak = analyze_linescan(x, filter(ones(1, 7)/7, 1, y), thr);
                DataStruct(end).peaks = npeak;
            else
                error('Invalid filter method');
            end
            
%             ythr = ones(1, length(x)) * thr;
%             subplot(3, 4, 4*(k-3)+c);
%             plot(x, y, x, ythr, 'r');
%             title(sprintf('%s - %d', clines(c).name, cell_number));
%             axis tight;
        end
    end
end

DataStruct = DataStruct(2:end);

% H = gcf;
% suptitle(exps(e).name);
% set(H, 'Position', [0, 0, 1600, 900]);
% fname = sprintf('choosethr_%s', exps(e).name);
% % saveas(H, ['figures\line_scan_cortical\', fname, '.fig']);
% saveas(H, ['figures\line_scan_cortical\', fname, '.png']);
% return;

%% printing
% clc;
% fprintf('Threshold = %d\n', thresholds(t));
% for k = 1:length(DataStruct)
%     fprintf('%15s | %02d | %3d (nof) | %3d (avgf)\n', ...
%         DataStruct(k).cline, DataStruct(k).cnum, DataStruct(k).peaks(t, 1), DataStruct(k).peaks(t, 2));
% end

%% plotting
plot_order = {'U2OS-Ctrl', 'U2OS-GTSE1', 'WT.13-GTSE1', '14A.07-GTSE1'};
MAX = 40;
N = length(plot_order);
peak_dumps = nan(N, MAX);
xend_dumps = nan(N, MAX);
exp_dumps = cell(N, MAX);
last_value = ones(1, N);
for k = 1:length(DataStruct)
    x = find(strcmp(DataStruct(k).cline, plot_order));
    exp_dumps{x, last_value(x)} = DataStruct(k).exp;
    peak_dumps(x, last_value(x)) = DataStruct(k).peaks;
    xend_dumps(x, last_value(x)) = DataStruct(k).xend;
    last_value(x) = last_value(x) + 1;
end
pppl_dumps = peak_dumps ./ xend_dumps;

%% plot separately for each experiment
for e = 1:length(exps)
    which_exp = strcmp(exp_dumps, exps(e).name);
    values = nansum(pppl_dumps .* which_exp, 2) ./ sum(which_exp, 2);
    for k = 1:length(plot_order)
        fprintf('%7s %15s %.4f\n', exps(e).name, plot_order{k}, values(k));
    end
end

%% scatter plot
figure(1); hold on;
for k = 1:length(plot_order)
    scatter(k * ones(1, MAX), peak_dumps(k, :), 20, 'k', 'filled');
    x = [k-0.2, k+0.2];
    y = [nanmean(peak_dumps(k, :)), nanmean(peak_dumps(k, :))];
    plot(x, y, 'r-', 'LineWidth', 2);
end

ylabel('No. of Peaks'); xlabel('Cell lines');
set(gca, 'XTick', 1:length(plot_order), 'XTickLabel', plot_order);
xlim([0, N+1]);

% save figure
H = gcf;
saveas(H, [FIG_DIR, 'peak.fig']);
saveas(H, [FIG_DIR, 'peak.png']);

%% scatter plot for path lengths
figure(2); hold on;
for k = 1:length(plot_order)
    scatter(k * ones(1, MAX), xend_dumps(k, :), 20, 'k', 'filled');
    x = [k-0.2, k+0.2];
    y = [nanmean(xend_dumps(k, :)), nanmean(xend_dumps(k, :))];
    plot(x, y, 'r-', 'LineWidth', 2);
end

ylabel('Path lengths'); xlabel('Cell lines');
set(gca, 'XTick', 1:length(plot_order), 'XTickLabel', plot_order);
xlim([0, N+1]);

% save figure
H = gcf;
saveas(H, [FIG_DIR, 'xend.fig']);
saveas(H, [FIG_DIR, 'xend.png']);

%% scatter plot for path lengths
figure(3); hold on;
for k = 1:length(plot_order)
    scatter(k * ones(1, MAX), pppl_dumps(k, :), 20, 'k', 'filled');
    x = [k-0.2, k+0.2];
    y = [nanmean(pppl_dumps(k, :)), nanmean(pppl_dumps(k, :))];
    plot(x, y, 'r-', 'LineWidth', 2);
end

ylabel('Peaks per path length'); xlabel('Cell lines');
set(gca, 'XTick', 1:length(plot_order), 'XTickLabel', plot_order);
xlim([0, N+1]);

% save figure
H = gcf;
saveas(H, [FIG_DIR, 'pppl.fig']);
saveas(H, [FIG_DIR, 'pppl.png']);

%% save value dumps as excel file
keyboard;


end


function npeak = analyze_linescan(x, y, thresh)

% compute difference between consecutive values
dy = diff(y);
% look for "1, -1" change in sign patterns
dy1 = dy(1:end-1);
dy2 = dy(2:end);
dyc = [sign(dy1); sign(dy2)];
peaks_at = find([0, all(dyc == [1; -1], 1), 0]);
% check if y is actually above threshold at these peak locations
thr_peaks = peaks_at(y(peaks_at) > thresh);
npeak = length(thr_peaks);

end


function plot_peaks(x, y, peaks)
plot(x, y);
keyboard;
end


function [x, y] = read_one_linescan(fname)
x = []; y = [];
fid = fopen(fname, 'r');
tline = fgetl(fid);
% assert(strcmp(tline, 'X,Y'), 'X,Y not init?');
assert(strcmp(tline, 'Distance_(µm),Gray_Value'));
tline = fgetl(fid);
while ischar(tline)
    % X, Y
    info = regexp(strtrim(tline), ',', 'split');
    x(end+1) = str2double(info{1});
    y(end+1) = str2double(info{2});
    tline = fgetl(fid);
end
fclose(fid);
end

function bg = read_background(fname)
bg = {};
fid = fopen(fname, 'r');
tline = fgetl(fid);
while ischar(tline)
    % cell-line-type cell-number b1 b2 b3
    % U2OS-Ctrl 4 124 129 132
    info = regexp(tline, ' ', 'split');
    bg = [bg; {info{1}, str2double(info{2}), mean(cellfun(@(x) str2double(x), info(3:5)))}];
    tline = fgetl(fid);
end
fclose(fid);
end


function keepdirs = only_folders(dirlist)
keepdirs = dirlist(arrayfun(@(x) x.isdir, dirlist));
end

