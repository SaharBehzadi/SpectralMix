function plot_cluster_grid (nx, ny, Nx, Ny)
% PLOT_CLUSTER_GRID
%
% $Id: plot_cluster_grid.m,v 1.1 2004/02/04 09:45:35 spapadim Exp $

% Plot horizontal lines
for x = cumsum(Nx)
  if x < nx
    line([0.5, ny+0.5], [x+0.5 x+0.5],'Color','r','LineStyle','--');
  end
end

% Plot vertical lines
for y = cumsum(Ny)
  if y < ny
    line([y+0.5 y+0.5], [0.5, nx+0.5],'Color','r','LineStyle','--');
  end
end

return

% $Log: plot_cluster_grid.m,v $
% Revision 1.1  2004/02/04 09:45:35  spapadim
% Initial checkin
%

