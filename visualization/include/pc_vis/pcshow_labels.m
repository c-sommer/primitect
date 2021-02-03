function pcshow_labels(pc, labels)
% PCSHOW_LABELS displays a point cloud where coloring is according to given labels
%
% PCSHOW_LABELS(pc, labels)
%
% Input:
% - pc: pointCloud object with normals
% - labels: set of labels for each point
%
% No output.

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer.

pcshow(pc.Location, labels)
num_colors = max(labels) - min(labels);
colormap([.75 .75 .75; distinguishable_colors(num_colors)])

% gray = .75 * [1 1 1];
% set(gcf,'color', gray)
% set(gca,'color', gray)
axis off