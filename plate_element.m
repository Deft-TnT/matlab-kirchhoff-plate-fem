function [Ke, fe, T, D] = plate_element(dx, dy, E, nu, h, q)
%PLATE_ELEMENT 四节点、十二自由度的 Kirchhoff 矩形板单元。
% 节点顺序为 (-a,-b)、(+a,-b)、(+a,+b)、(-a,+b)，其中 a=dx/2、b=dy/2。
% 每节点自由度：[w, theta_x, theta_y] = [w, dw/dy, -dw/dx]。
% 单元刚度矩阵 Ke = ∫B' D B dA，等效节点载荷 fe = ∫Nw' q dA。
% 对均匀矩形单元，每个方向采用三点高斯积分，可精确积分此处的多项式乘积。

validateattributes(dx, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(dy, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(E, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(h, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(nu, {'numeric'}, {'scalar','real','>=',0,'<',0.5});
validateattributes(q, {'numeric'}, {'scalar','real','finite'});

a = dx/2;
b = dy/2;
corners = [-1,-1; 1,-1; 1,1; -1,1];
C = zeros(12,12);
for n = 1:4
    [p, px, py] = plate_basis(corners(n,1), corners(n,2), a, b);
    rows = 3*n-2:3*n;
    C(rows,:) = [p; py; -px];
end
% T 将物理节点自由度映射到多项式系数。
% 求解此 12×12 线性方程组属于数值代数运算，并非调用外部有限元求解器。
T = C \ eye(12);

D0 = E*h^3/(12*(1-nu^2));
D = D0 * [1,nu,0; nu,1,0; 0,0,(1-nu)/2];

gp = [-sqrt(3/5), 0, sqrt(3/5)];
gw = [5/9, 8/9, 5/9];
Ke = zeros(12,12);
fe = zeros(12,1);
for iy = 1:3
    for ix = 1:3
        [p, ~, ~, pxx, pyy, pxy] = plate_basis(gp(ix), gp(iy), a, b);
        Nw = p*T;
        B = [-pxx; -pyy; -2*pxy]*T;
        dA = gw(ix)*gw(iy)*a*b;
        Ke = Ke + (B'*D*B)*dA;
        fe = fe + (Nw'*q)*dA;
    end
end
end
