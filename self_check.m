function checks = self_check(cfg,verbose)
%SELF_CHECK 对局部单元和整体板模型进行独立校验。
% 在 MATLAB 中运行：self_check
% 校验依据已知多项式和物理平衡关系，避免用同一实现的不同副本相互验证。

if nargin < 1 || isempty(cfg)
    cfg = plate_default_parameters();
end
if nargin < 2
    verbose = true;
end
% 实现校验固定采用非零载荷、x=0 边固支的已知算例。
% 课堂修改后的载荷和边界条件由 run_assignment 检查。
cfg.q = 1000;
cfg.bc = 'x0_clamped';
dx = 0.4;
dy = 0.3;
a = dx/2;
b = dy/2;
[Ke,fe,T,D] = plate_element(dx,dy,cfg.E,cfg.nu,cfg.h,cfg.q);
xy = [-a,-b; a,-b; a,b; -a,b];

checks.interpolation = 0;
for k = 1:4
    [p,px,py] = plate_basis(xy(k,1)/a,xy(k,2)/b,a,b);
    P = [p;py;-px];
    target = eye(12);
    rows = 3*k-2:3*k;
    checks.interpolation = max(checks.interpolation, ...
        norm(P*T-target(rows,:),'fro'));
end

% 根据课程第 6 章第 53 页的解析式计算均布压力对应的节点载荷。
fa = zeros(12,1);
for k = 1:4
    fa(3*k-2:3*k) = cfg.q*dx*dy * ...
        [1/4; -xy(k,2)/12; xy(k,1)/12];
end
checks.load_rel = norm(fe-fa)/norm(fa);
checks.load_sum_rel = abs(sum(fe(1:3:end))-cfg.q*dx*dy)/(cfg.q*dx*dy);
% 用十二项基函数中各单项式的精确积分检验虚功关系。
% 此检验独立于计算 fe 时使用的数值积分。
powers = [0,0; 1,0; 0,1; 2,0; 1,1; 0,2; ...
          3,0; 2,1; 1,2; 0,3; 3,1; 1,3];
Cpoly = zeros(12,12);
intpoly = zeros(12,1);
for j = 1:12
    pxpow = powers(j,1);
    pypow = powers(j,2);
    for k = 1:4
        xi = xy(k,1)/a;
        eta = xy(k,2)/b;
        rows = 3*k-2:3*k;
        Cpoly(rows(1),j) = xi^pxpow * eta^pypow;
        if pypow > 0
            Cpoly(rows(2),j) = pypow/b * xi^pxpow * eta^(pypow-1);
        end
        if pxpow > 0
            Cpoly(rows(3),j) = -pxpow/a * xi^(pxpow-1) * eta^pypow;
        end
    end
    if mod(pxpow,2) == 0 && mod(pypow,2) == 0
        intpoly(j) = 4*a*b / ((pxpow+1)*(pypow+1));
    end
end
checks.load_work_rel = norm(Cpoly'*fe-cfg.q*intpoly) / ...
    norm(cfg.q*intpoly);
checks.symmetry_rel = norm(Ke-Ke','fro')/norm(Ke,'fro');

rigid = zeros(12,3);
dq = zeros(12,1);
for k = 1:4
    x = xy(k,1); y = xy(k,2);
    rigid(3*k-2:3*k,1) = [1;0;0];         % 常值刚体位移：w = 1。
    rigid(3*k-2:3*k,2) = [x;0;-1];        % x 方向刚体转动：w = x。
    rigid(3*k-2:3*k,3) = [y;1;0];         % y 方向刚体转动：w = y。
    % 精确二次弯曲场：w=x^2+2xy+3y^2。
    dq(3*k-2:3*k) = [x^2+2*x*y+3*y^2; 2*x+6*y; -2*x-2*y];
end
checks.rigid_rel = norm(Ke*rigid,'fro')/ ...
    (norm(Ke,'fro')*norm(rigid,'fro'));
kappa = [-2;-6;-4];
Eexact = dx*dy*(kappa'*D*kappa);
checks.quadratic_energy_rel = abs(dq'*Ke*dq-Eexact)/Eexact;
eigvals = sort(eig((Ke+Ke')/2));
checks.min_eigen_relative = eigvals(1)/eigvals(end);
checks.fourth_eigen_relative = eigvals(4)/eigvals(end);

% 构造 2×2 单元片，并用精确二次场指定所有外侧节点自由度。
% 任意二次场都满足齐次薄板方程，因此中心节点的三个未知自由度
% 应能准确重现该二次场。
[Kpatch, xpatch, ypatch] = quadratic_patch_matrix(cfg);
center = 13:15;                         % 中心节点坐标为 (x,y)=(1,1)。
outer = setdiff(1:27,center);
checks.patch_names = {'1','x','y','x^2','y^2','xy'};
checks.patch_residual_rel = zeros(1,6);
checks.patch_center_rel = zeros(1,6);
for mode = 1:6
    dpatch = zeros(27,1);
    for node = 1:9
        x = xpatch(node);
        y = ypatch(node);
        switch mode
            case 1, dof = [1; 0; 0];
            case 2, dof = [x; 0; -1];
            case 3, dof = [y; 1; 0];
            case 4, dof = [x^2; 0; -2*x];
            case 5, dof = [y^2; 2*y; 0];
            case 6, dof = [x*y; x; -y];
        end
        dpatch(3*node-2:3*node) = dof;
    end
    rcenter = Kpatch(center,:)*dpatch;
    scale = norm(Kpatch(center,:),inf)*norm(dpatch,inf);
    checks.patch_residual_rel(mode) = norm(rcenter,inf)/scale;
    center_solution = -Kpatch(center,center) \ ...
        (Kpatch(center,outer)*dpatch(outer));
    checks.patch_center_rel(mode) = ...
        norm(center_solution-dpatch(center),inf) / ...
        max(1,norm(dpatch(center),inf));
end
checks.patch_max_rel = max([checks.patch_residual_rel, ...
                            checks.patch_center_rel]);

base = plate_solve(cfg,4,4);
double_load = cfg;
double_load.q = 2*cfg.q;
twice = plate_solve(double_load,4,4);
checks.linearity_rel = norm(twice.d-2*base.d)/norm(base.d);
checks.balance_rel = base.balance_rel;
checks.moment_y_balance_rel = base.moment_y_balance_rel;
checks.free_residual_rel = base.free_residual_rel;
checks.midline_symmetry_rel = base.symmetry_rel;
allclamped = cfg;
allclamped.bc = 'all_edges_clamped';
other = plate_solve(allclamped,4,4);
checks.all_clamped_wmax = other.wmax_abs;
checks.one_edge_wmax = base.wmax_abs;

assert(checks.interpolation < 1e-11,'Nodal interpolation identity failed.');
assert(checks.load_rel < 1e-11,'Analytic consistent load failed.');
assert(checks.load_sum_rel < 1e-11,'Element load sum failed.');
assert(checks.load_work_rel < 1e-11, ...
    'Consistent load fails exact monomial virtual-work identities.');
assert(checks.symmetry_rel < 1e-11,'Element stiffness is not symmetric.');
assert(checks.rigid_rel < 1e-11,'Rigid plate motion produces energy.');
assert(checks.quadratic_energy_rel < 1e-11, ...
    'Constant-curvature energy identity failed.');
assert(checks.min_eigen_relative > -1e-10, ...
    'Element stiffness has a negative energy mode.');
assert(checks.fourth_eigen_relative > 1e-10, ...
    'Element stiffness has extra zero-energy modes.');
assert(checks.patch_max_rel < 1e-10, ...
    'Quadratic 2-by-2 patch test failed.');
assert(checks.linearity_rel < 1e-10,'Load scaling failed.');
assert(checks.balance_rel < 1e-10,'Global vertical equilibrium failed.');
assert(checks.moment_y_balance_rel < 1e-10, ...
    'Global bending moment equilibrium failed.');
assert(checks.free_residual_rel < 1e-10,'Global free DOF residual failed.');
assert(checks.midline_symmetry_rel < 1e-10,'Y-midline symmetry failed.');
assert(checks.all_clamped_wmax < checks.one_edge_wmax, ...
    'Four-edge clamping should reduce peak displacement.');

fprintf('All element and plate checks passed.\n');
if verbose
    disp(checks);
end
end

function [Kpatch, xnode, ynode] = quadratic_patch_matrix(cfg)
% 不调用正式计算使用的网格和组装函数，单独组装测试单元片，
% 避免由相同的索引实现自我验证。
Kpatch = zeros(27,27);
[xx,yy] = ndgrid(0:2,0:2);
xnode = xx(:);
ynode = yy(:);
[Klocal,~,~,~] = plate_element(1,1,cfg.E,cfg.nu,cfg.h,0);
for ey = 1:2
    for ex = 1:2
        lower_left = ex+3*(ey-1);
        nodes = [lower_left, lower_left+1, lower_left+4, lower_left+3];
        dofs = reshape([3*nodes-2; 3*nodes-1; 3*nodes],1,[]);
        Kpatch(dofs,dofs) = Kpatch(dofs,dofs)+Klocal;
    end
end
end
