function [pole_a, pole_b, spindle_width] = read_poles_file(pole_fname, read_spindle_width)
%READ_POLES_FILE Parses through pole coordinates file to read important
%information from it.
% Main things to read: pole A and B coordinates, and spindle width
%

if ~exist('read_spindle_width', 'var'), read_spindle_width = true; end

% open file, ignore first 3/4 lines
fid = fopen(pole_fname);
tline = fgetl(fid);
while ischar(tline)
    if startsWith(strtrim(tline), 'Position X')
        break
    end
    tline = fgetl(fid);
end
% fgetl(fid);  % added for prometaphase, skip 4 lines

% next 2 lines are pole A and B info, process and get coords
pole_a = process_pole_line(fgetl(fid));
pole_b = process_pole_line(fgetl(fid));

if read_spindle_width
    % ignore next 2 lines
    fgetl(fid); fgetl(fid);

    % next line is spindle width in the first number, extract
    info = regexp(fgetl(fid), ',', 'split');
    spindle_width = str2double(info{1});

else
    spindle_width = nan;
end
    
% close file
fclose(fid);

end


function pole = process_pole_line(line)

% split line based on commas
info = regexp(line, ',', 'split');
% get the first 3 elements and convert to numbers
pole = [str2double(info{1}), str2double(info{2}), str2double(info{3})];

assert(strcmp(info{5}, 'MeasurementPoint'), 'Invalid pole line\n');

end

