function edges = plate_clamped_edges(bc)
%PLATE_CLAMPED_EDGES 解析课堂改参区中的固支边界条件。
% 例如：'x0_clamped'、'y1_clamped'、'x0_y1_clamped'、
% 'all_edges_clamped'。指定边上的 w、theta_x、theta_y 均固定。

if isstring(bc)
    assert(isscalar(bc),'cfg.bc must be one text value.');
    bc = char(bc);
end
assert(ischar(bc),'cfg.bc must be a character vector or string scalar.');
bc = lower(strtrim(bc));
if strcmp(bc,'all_edges_clamped')
    edges = {'x0','x1','y0','y1'};
    return;
end

suffix = '_clamped';
if numel(bc) <= numel(suffix) || ...
        ~strcmp(bc(end-numel(suffix)+1:end),suffix)
    error('Invalid cfg.bc: %s. Use e.g. x0_clamped or all_edges_clamped.',bc);
end
edges = strsplit(bc(1:end-numel(suffix)),'_');
allowed = {'x0','x1','y0','y1'};
if isempty(edges) || ~all(ismember(edges,allowed)) || ...
        numel(unique(edges)) ~= numel(edges)
    error('Invalid cfg.bc: %s. Edge names must be x0, x1, y0 or y1.',bc);
end
end
