function [S,M,f,Kx,Ky] = AssemblyMatrices(mesh, b, alpha, g, beta, edge, dir)

% Equation:
% alpha * Δu = Q
% Robin boundary condition:
% alpha * du/dn + beta * u = g

dx1   = mesh.n2c(mesh.e2n(:,2),1) - mesh.n2c(mesh.e2n(:,1),1);
dx2   = mesh.n2c(mesh.e2n(:,3),1) - mesh.n2c(mesh.e2n(:,1),1);
dy1   = mesh.n2c(mesh.e2n(:,2),2) - mesh.n2c(mesh.e2n(:,1),2);
dy2   = mesh.n2c(mesh.e2n(:,3),2) - mesh.n2c(mesh.e2n(:,1),2);

dx = mesh.n2c(mesh.n2e(edge,1),1) - mesh.n2c(mesh.n2e(edge,2),1);
dy = mesh.n2c(mesh.n2e(edge,1),2) - mesh.n2c(mesh.n2e(edge,2),2);

T_area   = (1/2)*(dx1.*dy2 - dx2.*dy1);
E_len = sqrt(dx.^2+dy.^2);

aa_x = (dy1-dy2).^2;        aa_y = (dx1-dx2).^2;
bb_x = dy2.^2;              bb_y = dx2.^2;
cc_x = dy1.^2;              cc_y = dx1.^2;
ab_x = dy2.*(dy1-dy2);      ab_y = dx2.*(dx1-dx2);
ac_x = -dy1.*(dy1-dy2);     ac_y = -dx1.*(dx1-dx2);
bc_x = -dy1.*dy2;           bc_y = -dx1.*dx2;


A_row = (1/2).*alpha.*(1./T_area).*[aa_x+aa_y,ab_x+ab_y,ac_x+ac_y, ...
                                ab_x+ab_y,bb_x+bb_y,bc_x+bc_y, ...
                                ac_x+ac_y,bc_x+bc_y,cc_x+cc_y];    
S_row = (1/6).*E_len.*beta(edge).*[2,1,1,2];
M_row = 2.*T_area.*(1/24).*[2,1,1,1,2,1,1,1,2];
K_row_x = (1/6).*repmat([dy1-dy2,dy2,-dy1],1,3,1);
K_row_y = (1/6).*repmat([dx2-dx1,-dx2,dx1],1,3,1);


jj = repmat(mesh.e2n,1,3);
ii = jj(:,[1,1,1,2,2,2,3,3,3]);
kk = mesh.e2n;
ll = ones(mesh.Nt,3);
mm = repmat(mesh.n2e(edge,:),1,2);
nn = mm(:,[1,1,2,2]);
oo = mesh.n2e(edge,:);
pp = ones(length(mesh.n2e(edge)),2);

A  = sparse(ii(:),jj(:),A_row(:));
B  = sparse(mm(:),nn(:),S_row(:),mesh.Nn,mesh.Nn);
M  = sparse(ii(:),jj(:),M_row(:));
Kx = sparse(ii(:),jj(:),K_row_x(:));
Ky = sparse(ii(:),jj(:),K_row_y(:));

f_row = repmat((2/3)*T_area.*b,1,3); 
f = sparse(kk,ll,f_row(:));

r_row = repmat((1/2).*E_len,1,2);
r     = g.*sparse(oo,pp,r_row(:),mesh.Nn,1);


d_ = zeros(mesh.Nn,1);
d_(dir) = g(dir);
d = -A*d_;

S = A+2*B;
f = f+2*r+d;

end