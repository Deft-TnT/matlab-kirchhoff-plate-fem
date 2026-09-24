function [K, F] = plate_assemble(nodes, elements, cfg)
%PLATE_ASSEMBLE 将单元矩阵组装为整体稀疏刚度矩阵 K 和载荷向量 F。
% 显式完成自由度映射、三元组组装及节点载荷累加等有限元核心步骤。
nn = size(nodes,1);
ne = size(elements,1);
ndof = 3*nn;
I = zeros(144*ne,1);
J = zeros(144*ne,1);
V = zeros(144*ne,1);
F = zeros(ndof,1);

for e = 1:ne
    conn = elements(e,:);
    xy = nodes(conn,:);
    dx = xy(2,1)-xy(1,1);
    dy = xy(4,2)-xy(1,2);
    [Ke, fe] = plate_element(dx,dy,cfg.E,cfg.nu,cfg.h,cfg.q);

    dofs = reshape(bsxfun(@plus,(1:3)',3*(conn-1)),[],1);
    slice = (e-1)*144+(1:144);
    I(slice) = repmat(dofs,12,1);
    J(slice) = kron(dofs,ones(12,1));
    V(slice) = Ke(:);
    F(dofs) = F(dofs)+fe;
end
K = sparse(I,J,V,ndof,ndof); % sparse() 会自动累加相同位置的三元组数值。
end
