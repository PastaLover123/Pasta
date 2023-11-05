
isovalue = 0.5*max(max(max(fields_d.psi)));

isosurface(fields_d.psi,isovalue)
isocaps(fields_d.psi,isovalue)

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'ztick',[])
set(gca,'zticklabel',[])

xlabel('X (fm)', 'fontSize', 12)
ylabel('Y (fm)', 'fontSize', 12)
zlabel('Z (fm)', 'fontSize', 12)
title('Density Field (dimensionless)', 'fontSize', 12)
colorbar
%set(get(colorbar,'title'), 'fontSize', 14)

xlim([1 64])
ylim([1 64])
zlim([1 64])

set(gca, 'linewidth', 1.5, 'fontsize', 16)
%set(gcf, 'Position',  [100, 100, 800, 700])

view(45, 45)
axis equal