function [] = PlotMesh(mesh)

figure;
patch('Faces',mesh.e2n,'Vertices',mesh.n2c,'FaceColor','none')
axis equal
xlim([min(mesh.n2c(:,1)) max(mesh.n2c(:,1))]);
ylim([min(mesh.n2c(:,2)) max(mesh.n2c(:,2))]);

end