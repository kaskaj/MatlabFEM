function [] = PlotData(mesh,field,file_name,save)

if nargin < 4
    save = 0;
end

fig1 = figure;
patch('Faces',mesh.e2n,'Vertices',mesh.n2c,'FaceVertexCData',field,'FaceColor','interp','EdgeColor','none');
h = colorbar;
t = get(h,'Limits');
set(h,'Ticks',linspace(t(1),t(2),5));
colormap jet;
axis equal;
xlim([min(mesh.n2c(:,1)) max(mesh.n2c(:,1))]);
ylim([min(mesh.n2c(:,2)) max(mesh.n2c(:,2))]);

if save
    saveas(fig1, file_name);
end

end

