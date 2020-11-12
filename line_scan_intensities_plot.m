function line_scan_intensities_plot()
clc;
close all;
run_one('G1', false);
run_one('G2', false);
run_one('NEBD', false);
run_one('Prometa', false);
run_one('Cyto', false);
end

function run_one(phase, use_raw_intensities)
data_dir = ['LineScanIntensityProfile/final_data/', phase, '/'];

%% load the data
% read cells
cells = dirnames(data_dir);
% for each cell
all_line_ints = {};
for c = 1:length(cells)
    cell_dir = [data_dir, num2str(cells(c)), '/'];
    sections = dirnames(cell_dir);
    % for each section
    for s = 1:length(sections)
        sec_dir = [cell_dir, num2str(sections(s)), '/'];
        % start loading
        fprintf('Processing phase %s | cell %02d | section %02d\n', phase, cells(c), sections(s));
        uid = sprintf('%s-%02d_%02d', phase, cells(c), sections(s));

        if use_raw_intensities
            % directly use intensities provided by ImageJ
            ints_dir = [sec_dir, 'Raw intensities/'];
            int_csvs = dir([ints_dir, '*.csv']);
            fprintf('Found %d raw intensity csv\n', length(int_csvs));
            this_line_ints = cell(length(int_csvs)/2, 4);
            for k = 1:length(int_csvs)
                line_num = str2double(int_csvs(k).name(1:end-4));
                row = ceil(line_num / 2);
                if mod(line_num, 2) == 0, col = 4;
                else col = 3;
                end
                fname = [ints_dir, int_csvs(k).name];
                data = readtable(fname);
                this_line_ints{row, col} = data.Y;
            end
            all_line_ints = [all_line_ints; this_line_ints];
            
        else
            % load image
            tif_fname = [sec_dir, sprintf('%s.tif', uid)];
            im = load_tif_2ch(tif_fname);
            % load coordinates
            coords_fname = [sec_dir, sprintf('Coordinates %s.csv', uid)];
            [x1, y1, x2, y2] = load_csv_coords(coords_fname);
            % get line intensities
            this_line_ints = get_line_intensities(im, x1, y1, x2, y2);
            % store coordinates, and unique identifiers for each line
            cell_coords = mat2cell([x1, y1, x2, y2], ones(1, length(x1)), 4);
            cell_uids = mat2cell(repmat(uid, length(x1), 1), ones(1, length(x1)));
            all_line_ints = [all_line_ints; [cell_uids, cell_coords, this_line_ints]];
        end
    end
end

%% align line intensities
num_lines = size(all_line_ints, 1);
fprintf('*** Aligning %d lines ***\n', num_lines);
% return;
line_lengths = cellfun(@(x) length(x), all_line_ints(:, 3));
[max_length, max_idx] = max(line_lengths);
padded_lines_c1 = nan(num_lines, max_length * 3);
padded_lines_c2 = nan(num_lines, max_length * 3);
left_pad = zeros(num_lines, 1);

% get maxima coordinate
[~, max_x] = max(all_line_ints{max_idx, 3});
for k = 1:num_lines
    c1_ints = all_line_ints{k, 3};
    c2_ints = all_line_ints{k, 4};
    [~, max_this] = max(c1_ints);
    left_pad(k) = max_length + (max_x - max_this);
    padded_lines_c1(k, left_pad(k) : left_pad(k) + length(c1_ints) - 1) = c1_ints;
    padded_lines_c2(k, left_pad(k) : left_pad(k) + length(c2_ints) - 1) = c2_ints;
end

% check line orientations
% figure;
% for k = 1:num_lines
%     clf;
%     plot(padded_lines_c1(k, :));
%     title(k);
%     pause;
% end
% keyboard;

%% create pretty plots
% create x-coordinates, such that 0 micron is at the center with peak for c1
step_x = 0.086;
zero_at = max_length + max_x - 1;
plot_x = -step_x * (zero_at - 1) : step_x : step_x * (size(padded_lines_c1, 2) - zero_at);

% divide by value at point closest to -1.7
start_x_at = abs(plot_x - (-1.7));
start_x_at = find(start_x_at == min(abs(start_x_at)));
padded_lines_c1 = padded_lines_c1 / nanmean(padded_lines_c1(:, start_x_at));
padded_lines_c2 = padded_lines_c2 / nanmean(padded_lines_c2(:, start_x_at));

% compute mean lines to plot "line"
avg_c1 = nanmean(padded_lines_c1, 1);
avg_c2 = nanmean(padded_lines_c2, 1);

%%% draw the main lines thicker
figure; hold on;
plot(plot_x, avg_c1, 'r', 'LineWidth', 1);
plot(plot_x, avg_c2, 'g', 'LineWidth', 1);

% get left/right start/end points for plotting
mean_factor = 0.2;  % at least 20% of the lines should have a non-nan value here
left_right = sum(~isnan(padded_lines_c1), 1) >= mean_factor * num_lines;
xlim_left = find(left_right, 1, 'first');
xlim_right = find(left_right, 1, 'last');
xlim([plot_x(xlim_left), plot_x(xlim_right)]);
grid on;

% keyboard;
assert(-1.7 >= plot_x(xlim_left), 'XLim left failed');
assert(2 <= plot_x(xlim_right), 'XLim right failed');
xlim([-1.7, 2]);
ylim([0.8, 1.4]);


% create title
if use_raw_intensities
    fname = sprintf('Phase-%s-Raw', phase);
else
    fname = sprintf('Phase-%s-Matlab', phase);
end
title(fname);

% keyboard;
% save figure
H = gcf;
saveas(H, ['figures\linescans\', fname, '.fig']);
saveas(H, ['figures\linescans\', fname, '.png']);

end


function folders = dirnames(data_dir)
% get the list of folder names, convert to numbers
folders = dir(data_dir);
folders(1:2) = [];
folders = arrayfun(@(x) str2double(x.name), folders);
end        


function im = load_tif_2ch(fname)
% load a 2 channel tif image
ch1 = imread(fname, 1);
ch2 = imread(fname, 2);
im = cat(3, ch1, ch2);
end


function [x1, y1, x2, y2] = load_csv_coords(fname)
% load coordinates
data = readtable(fname);
% +1 because ImageJ outputs are 0-indexed pixel coordinates
x = reshape(data.X, 2, length(data.X)/2) + 1;
x1 = x(1, :)';
x2 = x(2, :)';
y = reshape(data.Y, 2, length(data.Y)/2) + 1;
y1 = y(1, :)';
y2 = y(2, :)';
end


function line_ints = get_line_intensities(im, x1, y1, x2, y2, interp)
% use improfile to get line intensities
if ~exist('interp', 'var'), interp = 'bicubic'; end

num_lines = length(x1);
line_ints = {num_lines, 2};
line_distances = sqrt((x1 - x2).^2 + (y1 - y2).^2);

for k = 1:length(x1)
    N = round(line_distances(k)) + 1;
    xk = [x1(k), x2(k)];
    yk = [y1(k), y2(k)];
    line_ints{k, 1} = improfile(im(:, :, 1), xk, yk, N, interp);
    line_ints{k, 2} = improfile(im(:, :, 2), xk, yk, N, interp);
end

end



function [x1, y1, x2, y2] = line_xy_from_rect_angle(xywh, angles)
% xywh: B x 4, rectangle coordinates
% angle: B, -180 to 180, indicates the direction in which line was drawn
% returns the correct diagonal of the ROI, asymmetric diagonals

B = size(xywh, 1);
assert(length(angles) == B, 'Invalid rectangles and angles');

% xy-coordinates of rectangles
topleft = xywh(:, 1:2);
topright = [xywh(:, 1) + xywh(:, 3), xywh(:, 2)];
botleft = [xywh(:, 1), xywh(:, 2) + xywh(:, 4)];
botright = [xywh(:, 1) + xywh(:, 3), xywh(:, 2) + xywh(:, 4)];

% get x1y1, x2y2 for the lines
x1 = zeros(1, B); y1 = zeros(1, B);
x2 = zeros(1, B); y2 = zeros(1, B);
for k = 1:length(angles)
    if angles(k) > 0 && angles(k) < 90  % first quadrant BL --> TR
        x1(k) = botleft(k, 1);  y1(k) = botleft(k, 2);
        x2(k) = topright(k, 1); y2(k) = topright(k, 2);

    elseif angles(k) > 90 && angles(k) < 180  % second quadrant BR --> TL
        x1(k) = botright(k, 1); y1(k) = botright(k, 2);
        x2(k) = topleft(k, 1);  y2(k) = topleft(k, 2);

    elseif angles(k) < -90 && angles(k) > -180  % third quadrant TR --> BL
        x1(k) = topright(k, 1); y1(k) = topright(k, 2);
        x2(k) = botleft(k, 1);  y2(k) = botleft(k, 2);

    elseif angles(k) < 0 && angles(k) > -90  % fourth quadrant TL --> BR
        x1(k) = topleft(k, 1);  y1(k) = topleft(k, 2);
        x2(k) = botright(k, 1); y2(k) = botright(k, 2);

    end
end

% 0 indexed to 1 indexed
% x1 = x1 + 1;
% x2 = x2 + 1;
% y1 = y1 + 1;
% y2 = y2 + 1;
end
