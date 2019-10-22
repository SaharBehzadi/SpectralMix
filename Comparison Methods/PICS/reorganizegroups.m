function [A F Nx] = reorganizegroups(Afile, Ffile, tracefile,p,p2)


load(Afile);
load(Ffile);
load(tracefile);

xlabels = {'conservative','liberal'}; %  'neutral'
%xlabels = {'liberal','conservative'};

tr = trace{end};
tr.Nx

[nx ny] = size(F);

[Px, Py] = cc_perm_cohesive(tr.k, tr.l, tr.Qx, tr.Qy, F);


 figure;
 plot_binary_duo(F(Px,Py),1,2);
 plot_cluster_grid(nx, ny, tr.Nx, tr.Ny);
set(gca,'XTick',[1 2 ],'XTicklabel',xlabels(Py))
 
  figure;
 plot_binary_duo(A(Px,Px),1,1);
 plot_cluster_grid(nx, nx, tr.Nx, tr.Nx);
 show_layout_shaded(tr.Qx,tr.Qx,tr.Dnz,1)

 pause
 
Nx = tr.Nx(p)
Dnz = tr.Dnz(p,:);
Dnz = Dnz(:,p);

Qxt = tr.Qx;

 Qx = zeros(1, nx);
 for i=1:tr.k
     ind = find(Qxt==i);
     Qx(ind) = p2(i);
 end

 [Px, Py] = cc_perm_cohesive(tr.k, tr.l, Qx, tr.Qy, F);

 F = F(Px,Py);
 figure;
 plot_binary_duo(F,1,2);
 plot_cluster_grid(nx, ny, Nx, tr.Ny);
set(gca,'XTick',[1 2 ],'XTicklabel',xlabels(Py))
 
A = A(Px,Px);
  figure;
 plot_binary_duo(A,1,1);
 plot_cluster_grid(nx, nx, Nx, Nx);
 
 
  Qx = zeros(1, nx);
 for i=1:tr.k
     ind = find(Qxt==i);
     Qx(ind) = p2(i);
 end
  show_layout_shaded(Qx,Qx,Dnz,1)
 

end