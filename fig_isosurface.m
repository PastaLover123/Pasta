
isovalue = 0.5*max(max(max(fields_d.psi)));

isosurface(fields_d.psi,isovalue)
isocaps(fields_d.psi,isovalue)

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'ztick',[])
set(gca,'zticklabel',[])

xlim([1 64])
ylim([1 64])
zlim([1 64])

set(gcf, 'Position',  [100, 100, 800, 700])