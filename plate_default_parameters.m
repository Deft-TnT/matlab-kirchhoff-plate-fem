function cfg = plate_default_parameters()
%PLATE_DEFAULT_PARAMETERS 在此集中设置算例参数（国际单位制：N、m、Pa）。
% 题目要求一条边固支，但未指定具体位置。
% 默认固支 x=0 边；未列出的边均为自由边。
cfg.Lx = 1.0;
cfg.Ly = 1.0;
cfg.h = 0.010;
cfg.E = 210e9;
cfg.nu = 0.30;

% ===== 课堂改参区：只需修改下面三行，保存后运行 run_assignment =====
cfg.q = 1000.0;             % 单位为 N/m^2；正值表示沿 +z 方向加载。
cfg.mesh_sizes = [4, 8, 16, 32];
cfg.bc = 'x0_clamped';      % 也可选 'x1_clamped'、'y0_clamped'、
                           % 'x0_x1_clamped' 或 'all_edges_clamped'。
% 可将 x0、x1、y0、y1 中的任意边以 '_' 连接，再加上 '_clamped'。
% 例如，'x0_y1_clamped' 表示两条相邻边固支。
% ===== 课堂改参区结束 =====
end
