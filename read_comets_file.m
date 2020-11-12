function all_comets = read_comets_file(comet_fname)
%READ_COMETS_FILE Parses through comet file information. Includes several
%things, returns coordinates for now.
%

% open file, ignore first 4 lines
fid = fopen(comet_fname);
fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);

% collect comet coords
all_comets = [];
tline = fgetl(fid);
while ischar(tline)
    info = regexp(tline, ',', 'split');
    all_comets = [all_comets; ...
        [str2double(info{1}), str2double(info{2}), str2double(info{3})]];
    tline = fgetl(fid);
end

% print and bye
% fprintf('%s: Found %d comet coordinates\n', comet_fname, size(all_comets, 1));

fclose(fid);

end
