function [p, px, py, pxx, pyy, pxy] = plate_basis(xi, eta, a, b)
%PLATE_BASIS 计算矩形薄板单元十二项多项式及其导数。
% xi=x/a，eta=y/b；单元的物理坐标范围为 [-a,a] × [-b,b]。
% 将物理坐标无奇异地缩放为无量纲局部坐标后，各列依次对应
% [1,x,y,x^2,xy,y^2,x^3,x^2*y,x*y^2,y^3,x^3*y,x*y^3]。

p = [1, xi, eta, xi^2, xi*eta, eta^2, xi^3, xi^2*eta, ...
     xi*eta^2, eta^3, xi^3*eta, xi*eta^3];
p_xi = [0, 1, 0, 2*xi, eta, 0, 3*xi^2, 2*xi*eta, eta^2, 0, ...
        3*xi^2*eta, eta^3];
p_eta = [0, 0, 1, 0, xi, 2*eta, 0, xi^2, 2*xi*eta, 3*eta^2, ...
         xi^3, 3*xi*eta^2];
p_xixi = [0, 0, 0, 2, 0, 0, 6*xi, 2*eta, 0, 0, 6*xi*eta, 0];
p_etaeta = [0, 0, 0, 0, 0, 2, 0, 0, 2*xi, 6*eta, 0, 6*xi*eta];
p_xieta = [0, 0, 0, 0, 1, 0, 0, 2*xi, 2*eta, 0, 3*xi^2, 3*eta^2];

px = p_xi / a;
py = p_eta / b;
pxx = p_xixi / a^2;
pyy = p_etaeta / b^2;
pxy = p_xieta / (a*b);
end
