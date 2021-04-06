%% get_surfacedata()
% Reads surface data from CSV file.
% No input.
% Output:
%   - C : an (n x 3) matrix, where n is the number of surfaces, denoting
%         the location in body coordinates of each surface's centroid in m.
%   - N : an (n x 3) matrix denoting the unit outward-facing normal for
%         each surface.
%   - A : an (n x 1) matrix denoting the area for each surface, in m^2.

function [C, N, A] = get_surfacedata()
    centroid_file = 'CentroidData.csv';
    data = readtable(centroid_file, 'ReadVariableNames', false);
    C = [data.Var1 , data.Var2 , data.Var3];
    N = [data.Var4 , data.Var5 , data.Var6];
    A = data.Var7;
end