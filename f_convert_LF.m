% Converts the lead field or source structure of fieldtrip to the more
% convenient matrix and computes the projection
% 
%       [LFmatrix,LFproject] = f_convert_LF(data,specific)
% 
% Input
% grid:     structure containing lead field information after
%           ft_prepare_leadfield or source information after inverse model
%           ft_sourceanalysis (coming with fieldtrip)
% specific: ['LF' or 'source'] specification if the data input is the
%           leadfield or the already inverse modeled source
% 
% output
% LFmatrix: lead field (channels x sources or source matrix (sources x time
%           bins). Please keep in mind that x, y, and z directions are
%           concatinated (first 1/3 of the sources represent x, next 1/3 y...)
% LFproject: projection of leadfield or source: sqrt(x^2+y^2+z^2)
% 
function [LFmatrix,LFproject] = f_convert_LF(grid,specific)

if strcmp(specific,'LF')
    lead = grid.leadfield(grid.inside);
    LFMx = zeros(length(grid.label),numel(lead));
    LFMy = zeros(length(grid.label),numel(lead));
    LFMz = zeros(length(grid.label),numel(lead));
    for i = 1:numel(lead)
        LFMx(:,end+1) = lead{i}(:,1);
        LFMy(:,end+1) = lead{i}(:,2);
        LFMz(:,end+1) = lead{i}(:,3);
    end
    [LFmatrix] = [LFMx,LFMy,LFMz];
    [LFproject] = sqrt(LFMx.^2 + LFMy.^2 + LFMz.^2);
    
elseif strcmp(specific,'source')
    lead = grid.avg.mom(grid.inside);
    LFMx = zeros(numel(lead),length(grid.time));
    LFMy = zeros(numel(lead),length(grid.time));
    LFMz = zeros(numel(lead),length(grid.time));
    for i = 1:numel(lead)
        LFMx(end+1,:) = lead{i}(1,:);
        LFMy(end+1,:) = lead{i}(2,:);
        LFMz(end+1,:) = lead{i}(3,:);
    end
    [LFmatrix] = [LFMx;LFMy;LFMz];
    [LFproject] = sqrt(LFMx.^2 + LFMy.^2 + LFMz.^2);
else
    error('wrong specific for input of f_convert_LF');
end

end