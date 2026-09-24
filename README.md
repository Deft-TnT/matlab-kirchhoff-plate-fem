# MATLAB 板弯曲有限元工程

本工程用 MATLAB 实现均布横向载荷下的矩形薄板静力弯曲计算。程序自行完成网格生成、四节点板单元插值与积分、稀疏总体组装、固支约束、线性求解和反力恢复；没有调用现成有限元求解器。源码中的说明注释为中文。

## 运行

已在 MATLAB R2020a 验证。将 MATLAB 的 **Current Folder** 切换到本仓库目录，在命令窗口输入：

```matlab
run_assignment
```

程序先运行 `self_check`，随后计算各级网格。MATLAB 桌面会显示挠度曲面、中线挠度和网格收敛三张图；命令窗口会列出峰值、位置及反力。`results/` 中的表格、图片、文字摘要和 MAT 文件也会重新生成。要单独查看收敛表，可以运行：

```matlab
s = run_assignment;
s.convergence
```

## 课堂改参

只需打开 [`plate_default_parameters.m`](plate_default_parameters.m)，修改“课堂改参区”对应的一行，保存后重新运行 `run_assignment`：

```matlab
cfg.q = 1000.0;                % 均布载荷，N/m²；正值沿 +z
cfg.mesh_sizes = [4, 8, 16, 32]; % 每个数 n 对应 n×n 网格
cfg.bc = 'x0_clamped';         % 固支边界
```

`cfg.bc` 支持 `x0_clamped`、`x1_clamped`、`y0_clamped`、`y1_clamped`、`all_edges_clamped`，也支持 `x0_y1_clamped` 这样的多边固支组合；未指定的边为自由边。目前没有简支选项。几何尺寸、板厚和材料参数也集中在同一个文件上方。更多例子见[课堂改参速查](课堂改参速查.md)。

## 默认模型与结果

默认板长宽均为 1 m、厚 0.010 m，弹性模量 210 GPa、泊松比 0.30；`x=0` 边固支，其余三边自由，均布载荷为 1000 N/m²。采用四节点、每节点三个自由度的 Adini/ACM 型 Kirchhoff 矩形板单元，单元刚度和一致节点载荷用 3×3 Gauss 积分计算。

| 网格 | 峰值挠度（mm） |
|---|---:|
| 4×4 | 6.73678971 |
| 8×8 | 6.71709678 |
| 16×16 | 6.71290147 |
| 32×32 | 6.71205957 |

32×32 网格的峰值位于 `(1, 0.5)` m；外载合力约为 +1000 N，固支边竖向反力约为 −1000 N。`self_check` 检查单元插值、一致载荷虚功、刚体零能量、二次场能量、二次多项式 patch、线性倍增、总体平衡及自由自由度残差。已有的 `results/` 是默认参数的示例输出；改参后以重新运行的结果为准。

`free_edge_profile.csv` 始终记录 `x=Lx` 边的取样值。若把 `x1` 设为固支边，该边不再是自由边，文件中的 `xmax_edge_is_free` 会标为 `0`。

默认算例的 `w_max/h≈0.67`，因此结果应理解为**线性 Kirchhoff 模型预测**，不宜直接当作真实板的大挠度响应。网格间变化反映数值收敛，不能直接称为相对于精确解的误差。

## 文件

- `run_assignment.m`：主程序，执行验证、计算并生成结果。
- `plate_default_parameters.m`：集中管理输入参数。
- `plate_mesh.m`、`plate_basis.m`、`plate_element.m`、`plate_assemble.m`、`plate_solve.m`、`plate_clamped_edges.m`：有限元核心实现。
- `self_check.m`：独立数值和物理校验。
- `课堂改参速查.md`：现场改参示例。
- `results/`：默认算例的可复现输出。
