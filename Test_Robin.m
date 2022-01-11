clear;

%% Test problem definiton:

% alpha * Δu = b
% Upper and lower edge - Robin boundary condition: alpha * du/dn + beta * u = g
% Left edge  - Dirichlet boundary condition u = g
% Right edge - Neumann boundary condition alpha * du/dn = g
% Center square region - Load b

%% Load mesh
% Mesh was build using MESH2D package
% https://www.mathworks.com/matlabcentral/fileexchange/25555-mesh2d-delaunay-based-unstructured-mesh-generation

load('mesh.mat','mesh');

% PlotMesh(mesh);

%% Edges for boundary conditions

outer_nodes = mesh.n2c(:,1) == -0.5 | mesh.n2c(:,1) == 0.5 | ...
              mesh.n2c(:,2) == -0.5 | mesh.n2c(:,2) == 0.5;

outer_edges = ismember(mesh.n2e(:,1),find(outer_nodes));

right_nodes = mesh.n2c(:,1) == 0.5; 
right_edges = ismember(mesh.n2e(:,1),find(right_nodes));

left_nodes = mesh.n2c(:,1) == -0.5;

%% Materials
% Defined on elements

alpha = zeros(mesh.Nt,1);
alpha(mesh.reg == 1) = 1;
alpha(mesh.reg == 2) = 3;
alpha(mesh.reg == 3) = 5;

% PlotElem(mesh,alpha)

%% Load
% Defined on elements

b = zeros(mesh.Nt,1);
b(mesh.reg == 3) = 200;

% PlotElem(mesh,b)

%% Robin BC
% Defined on edges and at nodes

beta = 0.5*ones(mesh.Ne,1);
beta(right_edges) = 0;

g = 2.5*ones(mesh.Nn,1);
g(right_nodes) = 5;
g(left_nodes) = 10;

nodes_Dir = left_nodes; % nodes for Dirichlet BC

%% Assembly matrices

[S,M,f,Kx,Ky] = AssemblyMatrices(mesh, b, alpha, g, beta, outer_edges, nodes_Dir);

%% Solve

u = zeros(mesh.Nn,1);
u(~left_nodes) = S(~left_nodes,~left_nodes)\f(~left_nodes);
u(left_nodes) = 10;

%% Compute gradient

dudx = M\(Kx*u);
dudy = M\(Ky*u);

%% Plot

PlotData(mesh,u);
PlotData(mesh,dudx);
PlotData(mesh,dudy);

