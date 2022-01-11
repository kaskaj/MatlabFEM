function [] = PlotElem(mesh,field)

figure;
patch('Faces',mesh.e2n,'Vertices',mesh.n2c,'FaceVertexCData',field,'FaceColor','flat','EdgeColor','none')
h = colorbar;
t = get(h,'Limits');
set(h,'Ticks',linspace(t(1),t(2),5));
colormap jet;
axis equal;
xlim([min(mesh.n2c(:,1)) max(mesh.n2c(:,1))]);
ylim([min(mesh.n2c(:,2)) max(mesh.n2c(:,2))]);

end