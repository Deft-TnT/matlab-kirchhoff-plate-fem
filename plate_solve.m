function result = plate_solve(cfg, nx, ny)
%PLATE_SOLVE 组装整体方程，施加所选固支边界，并计算位移与支座反力。
[nodes,elements] = plate_mesh(cfg.Lx,cfg.Ly,nx,ny);
[K,F] = plate_assemble(nodes,elements,cfg);
ndof = numel(F);
tol = 1e-12*max(cfg.Lx,cfg.Ly);
edges = plate_clamped_edges(cfg.bc);
fixed_mask = false(size(nodes,1),1);
for k = 1:numel(edges)
    switch edges{k}
        case 'x0'
            fixed_mask = fixed_mask | abs(nodes(:,1)) < tol;
        case 'x1'
            fixed_mask = fixed_mask | abs(nodes(:,1)-cfg.Lx) < tol;
        case 'y0'
            fixed_mask = fixed_mask | abs(nodes(:,2)) < tol;
        case 'y1'
            fixed_mask = fixed_mask | abs(nodes(:,2)-cfg.Ly) < tol;
    end
end
fixed_nodes = find(fixed_mask);

fixed = reshape(bsxfun(@plus,(1:3)',3*(fixed_nodes(:)'-1)),[],1);
free = setdiff((1:ndof)',fixed);
d = zeros(ndof,1);
% 所有约束自由度的指定位移均为零；仅求解自由自由度对应的缩减方程组。
d(free) = K(free,free) \ F(free);
reaction = K*d-F; % 约束自由度上的残差即为支座反力。

w = d(1:3:end);
[wmax,imax] = max(abs(w));
reaction_z = sum(reaction(fixed(1:3:end)));
applied_z = sum(F(1:3:end));
W = reshape(w,nx+1,ny+1); % 行对应 x 方向，列对应 y 方向。

result.cfg = cfg;
result.nx = nx;
result.ny = ny;
result.nodes = nodes;
result.elements = elements;
result.K = K;
result.F = F;
result.d = d;
result.w = w;
result.reaction = reaction;
result.clamped_edges = edges;
result.fixed_dofs = fixed;
result.free_dofs = free;
result.wmax_abs = wmax;
result.wmax_signed = w(imax);
result.wmax_xy = nodes(imax,:);
result.wmax_over_h = wmax/cfg.h;
result.applied_z = applied_z;
result.reaction_z = reaction_z;
result.balance_rel = abs(reaction_z+applied_z)/max(abs(applied_z),1);
% 绕 y 轴的刚体虚转动对应 [w,theta_x,theta_y]=[-x,0,1]。
% 据此同时检查横向力及其共轭节点力矩的平衡。
result.applied_moment_y = sum(F(3:3:end)-nodes(:,1).*F(1:3:end));
result.reaction_moment_y = sum(reaction(3:3:end) ...
    -nodes(:,1).*reaction(1:3:end));
result.moment_y_balance_rel = abs(result.applied_moment_y ...
    +result.reaction_moment_y)/max(abs(result.applied_moment_y),1);
result.free_residual_rel = norm(reaction(free))/max(norm(F(free)),1);
% 仅当 y0、y1 两边同时固支或同时不固支时，才应关于 y=Ly/2 镜像对称。
result.ymid_symmetry_applicable = ...
    ismember('y0',edges) == ismember('y1',edges);
if result.ymid_symmetry_applicable
    result.symmetry_rel = max(abs(W(:)-reshape(W(:,end:-1:1),[],1))) ...
        /max(wmax,eps);
else
    result.symmetry_rel = nan;
end
end
