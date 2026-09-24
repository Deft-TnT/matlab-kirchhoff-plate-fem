function summary = run_assignment(cfg,outdir)
%RUN_ASSIGNMENT 完成板弯曲有限元算例的全部计算与结果展示。
% 在 MATLAB R2020a 或更新版本中，切换到本文件夹并输入 run_assignment。
% 结果默认保存在源文件旁的 results 子文件夹中。

if nargin < 1 || isempty(cfg)
    cfg = plate_default_parameters();
end
if nargin < 2 || isempty(outdir)
    outdir = fullfile(fileparts(mfilename('fullpath')),'results');
end
validateattributes(cfg.q,{'numeric'},{'scalar','real','finite'}, ...
    mfilename,'cfg.q');
sizes = cfg.mesh_sizes(:);
assert(isnumeric(sizes) && isvector(cfg.mesh_sizes) && ~isempty(sizes) ...
    && all(isfinite(sizes)) && all(sizes >= 1) ...
    && all(sizes == floor(sizes)) && all(diff(sizes) > 0), ...
    'cfg.mesh_sizes must be an increasing vector of positive integers.');
edges = plate_clamped_edges(cfg.bc);
checks = self_check(cfg,false);
if ~exist(outdir,'dir')
    mkdir(outdir);
end

ns = numel(sizes);
solutions = cell(ns,1);
peak_m = zeros(ns,1);
peak_mm = zeros(ns,1);
peak_over_h = zeros(ns,1);
ndof = zeros(ns,1);
nelem = zeros(ns,1);
force_balance_rel = zeros(ns,1);
moment_balance_rel = zeros(ns,1);
symmetry_rel = zeros(ns,1);
free_residual_rel = zeros(ns,1);
peak_x_m = zeros(ns,1);
peak_y_m = zeros(ns,1);

fprintf('Kirchhoff plate, clamped edges: %s; q = %.6g N/m^2\n', ...
    strjoin(edges,', '),cfg.q);
fprintf('Mesh     DOF       peak w (mm)      w/h        relative change\n');
for i = 1:ns
    n = sizes(i);
    r = plate_solve(cfg,n,n);
    solutions{i} = r;
    peak_m(i) = r.wmax_abs;
    peak_mm(i) = 1e3*r.wmax_abs;
    peak_over_h(i) = r.wmax_over_h;
    ndof(i) = numel(r.d);
    nelem(i) = size(r.elements,1);
    force_balance_rel(i) = r.balance_rel;
    moment_balance_rel(i) = r.moment_y_balance_rel;
    symmetry_rel(i) = r.symmetry_rel;
    free_residual_rel(i) = r.free_residual_rel;
    peak_x_m(i) = r.wmax_xy(1);
    peak_y_m(i) = r.wmax_xy(2);
    if i == 1
        fprintf('%3d x %-3d %6d    %12.8f     %.6f       (first mesh)\n', ...
            n,n,ndof(i),peak_mm(i),peak_over_h(i));
    else
        change = relative_change(peak_m(i),peak_m(i-1));
        fprintf('%3d x %-3d %6d    %12.8f     %.6f       %.6f%%\n', ...
            n,n,ndof(i),peak_mm(i),peak_over_h(i),100*change);
    end
end

% 不同网格结果之间的变化仅反映数值收敛情况，并非相对于解析解的误差。
% 最细网格的计算结果也只是网格对照值。
change_from_previous_pct = nan(ns,1);
for i = 2:ns
    change_from_previous_pct(i) = 100*relative_change(peak_m(i),peak_m(i-1));
end
difference_from_finest_pct = nan(ns,1);
if peak_m(end) > 0
    difference_from_finest_pct = 100*abs(peak_m-peak_m(end))/peak_m(end);
end
% Chang Fo-van（1980）发表于 Applied Mathematics and Mechanics 1(3):371-383 的结果表明：
% 对均布载荷作用下、nu=0.3 的方形悬臂板，自由边中点挠度为 0.13102*q*L^4/D。
% Tian、Zhong 与 Li（2011）发表于 Archives of Civil and Mechanical Engineering 11(4)
% 的表 2 对同一点给出解析系数 0.129（三位有效数字）。
% 两篇文献的数值不一致；下文与 Chang 结果的差距不能视为精确误差。
chang1980_peak_m = nan;
chang1980_gap_pct = nan(ns,1);
if strcmp(cfg.bc,'x0_clamped') && cfg.q > 0 ...
        && abs(cfg.Lx-cfg.Ly) < 1e-12*cfg.Lx ...
        && abs(cfg.nu-0.3) < 1e-12
    bending_D = cfg.E*cfg.h^3/(12*(1-cfg.nu^2));
    chang1980_peak_m = 0.13102*cfg.q*cfg.Lx^4/bending_D;
    chang1980_gap_pct = 100*abs(peak_m-chang1980_peak_m) ...
        /chang1980_peak_m;
end
mesh_n = sizes;
convergence = table(mesh_n,mesh_n,nelem,ndof,peak_m,peak_mm, ...
    peak_over_h,peak_x_m,peak_y_m,change_from_previous_pct, ...
    difference_from_finest_pct,chang1980_gap_pct,force_balance_rel,moment_balance_rel, ...
    free_residual_rel,symmetry_rel,'VariableNames', ...
    {'nx','ny','elements','dofs','peak_w_m','peak_w_mm','peak_w_over_h', ...
    'peak_x_m','peak_y_m','change_from_previous_pct', ...
    'difference_from_finest_pct','chang1980_gap_pct','force_balance_rel', ...
    'moment_balance_rel','free_residual_rel','ymid_symmetry_rel'});
writetable(convergence,fullfile(outdir,'mesh_convergence.csv'));

finest = solutions{end};
% 边界条件与线性关系的辅助检查采用用户所选网格范围内的网格。
check_n = min(16,sizes(end));
idx_check = find(sizes==check_n,1);
if isempty(idx_check)
    active_check = plate_solve(cfg,check_n,check_n);
else
    active_check = solutions{idx_check};
end
one_edge_cfg = cfg;
one_edge_cfg.bc = 'x0_clamped';
if strcmp(cfg.bc,one_edge_cfg.bc)
    one_edge_check = active_check;
else
    one_edge_check = plate_solve(one_edge_cfg,check_n,check_n);
end
all_edges_cfg = cfg;
all_edges_cfg.bc = 'all_edges_clamped';
if strcmp(cfg.bc,all_edges_cfg.bc)
    all_edges_check = active_check;
else
    all_edges_check = plate_solve(all_edges_cfg,check_n,check_n);
end
if cfg.q ~= 0
    assert(all_edges_check.wmax_abs < one_edge_check.wmax_abs, ...
        'Four-edge clamping did not reduce the x=0 cantilever deflection.');
end

% 将载荷加倍后重新计算，以检查响应是否满足线性比例关系。
double_cfg = cfg;
double_cfg.q = 2*cfg.q;
double_load = plate_solve(double_cfg,check_n,check_n);
linearity_rel = norm(double_load.d-2*active_check.d) ...
    /max(norm(active_check.d),1);
assert(linearity_rel < 1e-10,'The doubled-load response is not linear.');

% 压力合力为 q*Lx*Ly；按自由度 [w, w_y, -w_x] 的符号约定，
% 绕 y 轴的载荷合力矩为 -q*Ly*Lx^2/2。
expected_force_n = cfg.q*cfg.Lx*cfg.Ly;
expected_moment_nm = -cfg.q*cfg.Ly*cfg.Lx^2/2;
assert(abs(finest.applied_z-expected_force_n) ...
    /max(abs(expected_force_n),1) < 1e-10);
assert(abs(finest.applied_moment_y-expected_moment_nm) ...
    /max(abs(expected_moment_nm),1) < 1e-10);

summary.cfg = cfg;
summary.convergence = convergence;
summary.checks = checks;
summary.finest = finest;
summary.check_mesh_n = check_n;
summary.active_check = active_check;
summary.one_edge_check = one_edge_check;
summary.all_edges_check = all_edges_check;
summary.linearity_rel = linearity_rel;
summary.expected_force_n = expected_force_n;
summary.expected_moment_nm = expected_moment_nm;
summary.chang1980_peak_m = chang1980_peak_m;
summary.tian2011_rounded_coefficient = nan;
if ~isnan(chang1980_peak_m)
    summary.tian2011_rounded_coefficient = 0.129;
end
save(fullfile(outdir,'study_results.mat'),'summary','-v7');
write_text_summary(summary,outdir);
write_free_edge_profile(summary,outdir);
make_figures(cfg,solutions,outdir);

fprintf('Finest peak: %.8f mm at (x,y) = (%.3f, %.3f) m\n', ...
    peak_mm(end),finest.wmax_xy(1),finest.wmax_xy(2));
if ~isnan(chang1980_peak_m)
    fprintf('Chang (1980): %.8f mm; gap: %.4f%%. Tian (2011): coeff. 0.129 (3 s.f.)\n', ...
        1000*chang1980_peak_m,chang1980_gap_pct(end));
end
fprintf('Vertical reaction: %.8f N; applied load: %.8f N\n', ...
    finest.reaction_z,finest.applied_z);
fprintf('Results saved in: %s\n',outdir);
end

function value = relative_change(current,previous)
if current == 0
    value = nan;
else
    value = abs(current-previous)/abs(current);
end
end

function write_text_summary(s,outdir)
fname = fullfile(outdir,'summary.txt');
fid = fopen(fname,'w');
if fid < 0
    error('Unable to open %s for writing.',fname);
end
cleanup = onCleanup(@() fclose(fid));
r = s.finest;
c = s.all_edges_check;
fprintf(fid,'STATIC KIRCHHOFF PLATE FEM ASSIGNMENT\n');
fprintf(fid,'Clamped edges: %s; all other edges free.\n', ...
    strjoin(r.clamped_edges,', '));
fprintf(fid,'cfg.bc = %s; mesh sizes = %s\n', ...
    char(s.cfg.bc),mat2str(s.cfg.mesh_sizes));
fprintf(fid,'Lx = %.6f m; Ly = %.6f m; h = %.6f m\n', ...
    s.cfg.Lx,s.cfg.Ly,s.cfg.h);
fprintf(fid,'E = %.6g Pa; nu = %.6f; q = %.6f N/m^2 (+w)\n', ...
    s.cfg.E,s.cfg.nu,s.cfg.q);
fprintf(fid,'Element: 4-node rectangular Adini/ACM Kirchhoff bending element.\n');
fprintf(fid,'Per-node DOFs: [w, dw/dy, -dw/dx]. 3x3 Gauss stiffness/load.\n');
fprintf(fid,'\nFINEST MESH\n');
fprintf(fid,'Elements = %d; nodes = %d; total DOFs = %d\n', ...
    size(r.elements,1),size(r.nodes,1),numel(r.d));
fprintf(fid,'Peak absolute nodal deflection = %.10f m = %.8f mm\n', ...
    r.wmax_abs,1000*r.wmax_abs);
fprintf(fid,'Peak coordinates: x = %.6f m, y = %.6f m\n',r.wmax_xy);
fprintf(fid,'Peak w/h = %.6f\n',r.wmax_over_h);
if ~isnan(s.chang1980_peak_m)
    D0 = s.cfg.E*s.cfg.h^3/(12*(1-s.cfg.nu^2));
    fe_coefficient = r.wmax_abs*D0/(s.cfg.q*s.cfg.Lx^4);
    fprintf(fid,'This FEM dimensionless coefficient = %.9f\n',fe_coefficient);
    fprintf(fid,'Chang (1980) analytic coefficient = 0.13102 (%.8f mm)\n', ...
        1000*s.chang1980_peak_m);
    fprintf(fid,'Gap from Chang (1980) = %.6f%%\n', ...
        s.convergence.chang1980_gap_pct(end));
    fprintf(fid,'Tian, Zhong & Li (2011) analytic coefficient = 0.129\n');
    fprintf(fid,'(three significant figures); this FEM rounds to 0.129.\n');
    fprintf(fid,'The two published analytic values disagree. Do not\n');
    fprintf(fid,'interpret the gap from Chang as verified exact error.\n');
    fprintf(fid,'Chang source: Applied Mathematics and Mechanics 1(3):371-383.\n');
    fprintf(fid,'https://www.amm.shu.edu.cn/EN/abstract/abstract14282.shtml\n');
    fprintf(fid,'Tian source: Archives of Civil and Mechanical Engineering\n');
    fprintf(fid,'11(4):1043-1052, DOI:10.1016/S1644-9665(12)60094-6.\n');
end
fprintf(fid,'Applied transverse force = %.10f N\n',r.applied_z);
fprintf(fid,'Support transverse reaction = %.10f N\n',r.reaction_z);
fprintf(fid,'Applied y-moment = %.10f N m\n',r.applied_moment_y);
fprintf(fid,'Support y-moment = %.10f N m\n',r.reaction_moment_y);
fprintf(fid,'Relative force equilibrium error = %.3e\n',r.balance_rel);
fprintf(fid,'Relative moment equilibrium error = %.3e\n',r.moment_y_balance_rel);
fprintf(fid,'Relative free-DOF residual = %.3e\n',r.free_residual_rel);
if r.ymid_symmetry_applicable
    fprintf(fid,'Relative y-midline symmetry error = %.3e\n',r.symmetry_rel);
else
    fprintf(fid,'Y-midline symmetry check: not applicable to these supports.\n');
end
fprintf(fid,'\nOTHER CHECKS\n');
fprintf(fid,'%dx%d all-edges-clamped peak = %.8f mm\n', ...
    s.check_mesh_n,s.check_mesh_n,1000*c.wmax_abs);
fprintf(fid,'%dx%d x=0 one-edge-clamped peak = %.8f mm\n', ...
    s.check_mesh_n,s.check_mesh_n,1000*s.one_edge_check.wmax_abs);
fprintf(fid,'%dx%d active-boundary peak = %.8f mm\n', ...
    s.check_mesh_n,s.check_mesh_n,1000*s.active_check.wmax_abs);
fprintf(fid,'%dx%d doubled-load linearity relative error = %.3e\n', ...
    s.check_mesh_n,s.check_mesh_n,s.linearity_rel);
fprintf(fid,'\nInterpretation: mesh-change percentages are convergence indicators,\n');
fprintf(fid,'not exact analytical errors. The element is nonconforming.\n');
if r.wmax_over_h >= 0.1
    fprintf(fid,'Peak w/h near %.2f makes the small-deflection assumption\n', ...
        r.wmax_over_h);
    fprintf(fid,'questionable for physical predictions; the result is the\n');
    fprintf(fid,'specified linear Kirchhoff model prediction only.\n');
else
    fprintf(fid,'Peak w/h = %.3g; this is a linear Kirchhoff model prediction.\n', ...
        r.wmax_over_h);
end
clear cleanup;
end

function write_free_edge_profile(s,outdir)
% 在 x=Lx 处取样；仅当 x1 未固支时，该边才是自由边。
r = s.finest;
xmax_edge_is_free = ~ismember('x1',r.clamped_edges);
y_fraction = [0;0.25;0.5;0.75;1];
y_m = y_fraction*s.cfg.Ly;
W = reshape(r.w,r.nx+1,r.ny+1);
w_m = interp1(linspace(0,s.cfg.Ly,r.ny+1),W(end,:),y_m);
w_mm = 1000*w_m;
coefficient = nan(size(w_m));
tian2011_rounded = nan(size(w_m));
chang1980 = nan(size(w_m));
if ~isnan(s.chang1980_peak_m)
    D0 = s.cfg.E*s.cfg.h^3/(12*(1-s.cfg.nu^2));
    coefficient = w_m*D0/(s.cfg.q*s.cfg.Lx^4);
    tian2011_rounded = [0.127;0.129;0.129;0.129;0.127];
    chang1980 = [0.12933;0.13056;0.13102;0.13056;0.12933];
end
profile = table(y_m,w_m,w_mm,repmat(xmax_edge_is_free,size(y_m)), ...
    coefficient,tian2011_rounded,chang1980, ...
    'VariableNames',{'y_m','fem_w_m','fem_w_mm','xmax_edge_is_free','fem_coefficient', ...
    'tian2011_coefficient_3sf','chang1980_coefficient'});
writetable(profile,fullfile(outdir,'free_edge_profile.csv'));
end

function make_figures(cfg,solutions,outdir)
r = solutions{end};
n = r.nx;
X = repmat(linspace(0,cfg.Lx,n+1)',1,n+1);
Y = repmat(linspace(0,cfg.Ly,n+1),n+1,1);
W = reshape(r.w,n+1,n+1);
show_plots = usejava('desktop');

f = case_figure('plate_fea_surface','Plate deflection', ...
    [100,100,920,650],show_plots);
surf(X,Y,1000*W,'EdgeColor','none');
colormap(parula(256));
shading interp;
colorbar;
view(38,30);
axis tight;
grid on;
xlabel('x (m)'); ylabel('y (m)'); zlabel('w (mm)');
title(sprintf('Clamped %s: deflection, %d x %d mesh', ...
    char(cfg.bc),n,n),'Interpreter','none');
print(f,fullfile(outdir,'deflection_surface.png'),'-dpng','-r180');
if ~show_plots, close(f); end

f = case_figure('plate_fea_centerline','Plate centerline', ...
    [100,100,900,550],show_plots);
hold on;
colors = lines(numel(solutions));
for k = 1:numel(solutions)
    t = solutions{k};
    T = reshape(t.w,t.nx+1,t.ny+1);
    centerline = interp1(linspace(0,cfg.Ly,t.ny+1),T',cfg.Ly/2)';
    plot(linspace(0,cfg.Lx,t.nx+1),1000*centerline, ...
        '-o','Color',colors(k,:),'LineWidth',1.5, ...
        'MarkerSize',4,'DisplayName',sprintf('%d x %d',t.nx,t.ny));
end
hold off;
grid on;
xlabel('x (m)'); ylabel(sprintf('w at y = %.3g m (mm)',cfg.Ly/2));
title(sprintf('Deflection along y-midline; clamped %s', ...
    char(cfg.bc)),'Interpreter','none');
legend('Location','northwest');
print(f,fullfile(outdir,'centerline_deflection.png'),'-dpng','-r180');
if ~show_plots, close(f); end

peaks = cellfun(@(s) 1000*s.wmax_abs,solutions);
mesh = cellfun(@(s) s.nx,solutions);
f = case_figure('plate_fea_convergence','Plate mesh convergence', ...
    [100,100,900,550],show_plots);
plot(mesh,peaks,'o-','LineWidth',1.8,'MarkerSize',7);
grid on;
xlabel('Elements per side (n)'); ylabel('Peak nodal deflection (mm)');
title(sprintf('Mesh convergence; clamped %s',char(cfg.bc)), ...
    'Interpreter','none');
xticks(mesh);
for k = 1:numel(mesh)
    text(mesh(k),peaks(k),sprintf('  %.5f',peaks(k)), ...
        'VerticalAlignment','bottom');
end
print(f,fullfile(outdir,'mesh_convergence.png'),'-dpng','-r180');
if ~show_plots, close(f); end
if show_plots, drawnow; end
end

function f = case_figure(tag,name,position,show_plots)
% 关闭本程序先前生成的同名图窗，便于课堂重新运行时查看新结果。
old = findall(0,'Type','figure','Tag',tag);
if ~isempty(old), close(old); end
if show_plots
    visibility = 'on';
else
    visibility = 'off';
end
f = figure('Visible',visibility,'Color','w','Position',position, ...
    'Tag',tag,'Name',name,'NumberTitle','off');
end
