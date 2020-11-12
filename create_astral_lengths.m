function create_astral_lengths(dirname)
%CREATE_ASTRAL_LENGTHS Processes all comet tracks and poles
clc;
dirname = '../data/Exp03';

read_spindle_width = false;
all_folders = dir(dirname);
all_folders(1:2) = [];

%%% read spindle width file separately?
all_spindle_widths = read_all_spindle_widths(fullfile(dirname, 'all_spindle_widths.csv'));
% ext_spw = nan;  % set to nan to read from pole file

for k = 1:length(all_folders)
    % check that the thing is not a file
    if ~all_folders(k).isdir, continue; end

    % Found an actual data folder, process cell lines
    fprintf('=============================================\n');
    fprintf('%s\n', all_folders(k).name);
    cell_line_folders = dir(fullfile(dirname, all_folders(k).name));
    cell_line_folders(1:2) = [];

    % Get the list of actual cells
    cells = {};
    for c = 1:length(cell_line_folders)
        if ~cell_line_folders(c).isdir, continue; end
        cell = strrep(cell_line_folders(c).name, '_Statistics', '');
        cell = strrep(cell, '_comet', '');
        cells{end+1} = cell;
    end
    cells = unique(cells);
    fprintf('Found %d cells.\n', length(cells));
    
    % For each cell, create file templates and process
    for c = 1:length(cells)
        mat_fname = fullfile(dirname, all_folders(k).name, [cells{c}, '.mat']);
        % if mat file exists, continue
%         if exist(mat_fname, 'file')
%             fprintf('MAT exists. %s | %s\n', all_folders(k).name, cells{c});
%             continue;
%         end

        % generate pole and comet filenames
        pole_fname = fullfile(dirname, all_folders(k).name, [cells{c}, '_Statistics'], [cells{c}, '_Detailed.csv']);
        comet_fname = fullfile(dirname, all_folders(k).name, [cells{c}, '_comet_Statistics'], [cells{c}, '_comet_Detailed.csv']);
        % verify files exist
        if ~exist(pole_fname, 'file')
            fprintf(2, 'Pole file does not exist! %s %s\n', all_folders(k).name, cells{c});
            continue;
        end
        if ~exist(comet_fname, 'file')
            fprintf(2, 'Comet file does not exist! %s %s\n', all_folders(k).name, cells{c});
            continue;
        end

        % force feed spindle width, read from external csv file
        this_folder = all_folders(k).name;
        [~, mfn, ~] = fileparts(mat_fname);
        split_mfn = regexp(mfn, ' - ', 'split');
        this_type = split_mfn{1}(length(this_folder)+2:end);
        split_mfn2 = regexp(split_mfn{2}, '_', 'split');
        this_num = str2double(split_mfn2{1});
        idx = find(strcmp(all_spindle_widths(:, 1), this_folder) & ...
                   strcmp(all_spindle_widths(:, 2), this_type) & ...
                   [all_spindle_widths{:, 3}]' == this_num);
        if isempty(idx)
            fprintf(2, 'Failed to find spindle width!\n');
            keyboard;
        else
            ext_spw = all_spindle_widths{idx, 4};
        end

        % process!
        fprintf('Processing: %s | %s\n', all_folders(k).name, cells{c});
        process_one_pair(mat_fname, pole_fname, comet_fname, read_spindle_width, ext_spw);
    end
end

end


function [comet_lengths, astral_lengths] = ...
    process_one_pair(mat_fname, pole_fname, comet_fname, read_spindle_width, ext_spw)

if ~exist('read_spindle_width', 'var'), read_spindle_width = true; end

% read poles and comet files
[pole_a, pole_b, spindle_width] = read_poles_file(pole_fname, read_spindle_width);
if ~isnan(ext_spw) && isnan(spindle_width)
    spindle_width = ext_spw;
end
all_comets = read_comets_file(comet_fname);
if any(isnan(pole_a)), error('Pole A is nan!'); end
if any(isnan(pole_b)), error('Pole B is nan!'); end
% if any(isnan(spindle_width)), error('Spindle width is nan!'); end

% get spindle information
spindle_length = sqrt(sum((pole_a - pole_b).^2));
spindle_angle = atan((spindle_width/2)/(spindle_length/2)) * 180/pi;

% pole-comet distances
pole_comet_dists = [pdist2(pole_a, all_comets); pdist2(pole_b, all_comets)]';
comet_lengths = min(pole_comet_dists, [], 2);
comet_to_other_pole = max(pole_comet_dists, [], 2);

% cosine rule, find angle between comet and spindle axis
num = spindle_length.^2 + comet_lengths.^2 - comet_to_other_pole.^2;
den = 2 * spindle_length * comet_lengths;
comet_angle = acos(num./den) * 180/pi;

% find astrals and their length
which_astrals = comet_angle > spindle_angle;
astral_lengths = comet_lengths(which_astrals);

% save information as matfile
save(mat_fname, 'astral_lengths', 'comet_lengths');
fprintf('\b\t%4d astrals, %4d comets\n', size(astral_lengths, 1), size(comet_lengths, 1));

end


function collect = read_all_spindle_widths(fname)

collect = {};
fid = fopen(fname, 'r');
tline = fgetl(fid);
while ischar(tline)
    spl = regexp(tline, ',', 'split');
    collect = [collect; {spl{1}, spl{2}, str2double(spl{3}), str2double(spl{4})}];
    tline = fgetl(fid);
end

end





