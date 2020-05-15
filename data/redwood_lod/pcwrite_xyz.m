function pcwrite_xyz(pc, filename)
% PCWRITE_XYZ
%
% See also: PCWRITE

% C. Sommer 05/2019

[fileID, err_msg] = fopen(filename, 'w');
disp(err_msg);
if size(pc.Normal, 2)<1
    fprintf(fileID, '%5.4f %5.4f %5.4f\n', pc.Location');
else
    fprintf(fileID, '%5.4f %5.4f %5.4f %5.4f %5.4f %5.4f\n', [pc.Location pc.Normal]');
end
fclose(fileID);