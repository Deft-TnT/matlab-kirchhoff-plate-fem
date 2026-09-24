function [nodes, elements] = plate_mesh(Lx, Ly, nx, ny)
%PLATE_MESH 生成结构化矩形网格，单元节点按逆时针排列。
% 节点编号先沿 x 方向递增，再沿 y 方向递增；坐标为整体坐标（m）。
validateattributes(Lx, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(Ly, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(nx, {'numeric'}, {'scalar','integer','positive'});
validateattributes(ny, {'numeric'}, {'scalar','integer','positive'});

[X,Y] = meshgrid(linspace(0,Lx,nx+1), linspace(0,Ly,ny+1));
Xt = X.';
Yt = Y.';
nodes = [Xt(:), Yt(:)];

elements = zeros(nx*ny,4);
e = 0;
for j = 1:ny
    for i = 1:nx
        e = e+1;
        n1 = (j-1)*(nx+1)+i;
        n2 = n1+1;
        n4 = n1+(nx+1);
        n3 = n4+1;
        elements(e,:) = [n1,n2,n3,n4];
    end
end
end
