function [cA A2 F2] = test_reality(name, plotorg, mincost)

%load(strcat('D:/DATA/new_reality_mining_data/',name,'.mat'))
%load('D:/DATA/new_reality_mining_data/F.mat')
%A=device;
%A = triu(A,1)+tril(A,-1);

load('data/A_call.mat')
load('data/F_call.mat')
xlabels = {'prof','grad','grad-1','ugrad','ugrad-1','staff','sloan'};





if(plotorg)
% Plot the original connectivity and feature matrices
figure; plot_binary_duo(A,0,1);
figure; plot_binary_duo(F,0,2);
set(gca,'XTick',1:7,'XTicklabel',xlabels)
xticklabel_rotate

%[Px, Py] = cc_perm_cohesive(1, 1, ones(1,size(A,1)), ones(1,size(F,2)), F);
% figure; plot_binary_duo(A(Px,Px),0,1);
% figure; plot_binary_duo(F(Px,Py),0,2);
% set(gca,'XTick',1:7,'XTicklabel',xlabels)
% xticklabel_rotate
% pause
end
 

% k=3;
while (1)
 % permute randomly --due to local optimum
%xperm = randperm(size(A,1)); 
%A = A(xperm, xperm);
%F = F(xperm,:);
% figure; plot_binary_duo(A,1);
% figure; plot_binary_duo(F,0);
% set(gca,'XTick',[1 2 3],'XTicklabel',xlabels)
%  pause
 

% [k,l,Nx,Ny,Qx,Qy,Dnz, cA, c2] = cc_search(A,'hellscream',1,name);
% %end
% 
% Nx
% [Px, Py] = cc_perm(k, l, Qx, Qy);
%   
%  [nx ny] = size(A);
%  figure;
%  plot_binary_duo(F(Px,:),0);
%  set(gca,'XTicklabel',xlabels)
%  for x = cumsum(Nx)
%   if x < nx
%    line([0.5, ny+0.5], [x+0.5 x+0.5],'LineStyle','--');
%   end
%  end
%  
%  figure;
%  plot_binary_duo(A(Px,Py),1);
%  plot_cluster_grid(nx, ny, Nx, Ny);
%  pause

 %%

 [k,l,Nx,Ny,Qx,Qy,Dnz, DnzF, cA, c2] = cc_search_duo(A,F,name);
 
 if(cA<mincost)
[nx ny] = size(F);
% [Px, Py] = cc_perm(k, l, Qx, Qy);
[Px, Py] = cc_perm_cohesive(k, l, Qx, Qy, F);

F2 = F(Px,Py);
 figure;
 plot_binary_duo(F2,1,3);
 plot_cluster_grid(nx, ny, Nx, Ny);
 set(gca,'XTick',1:7,'XTicklabel',xlabels(Py))
 xticklabel_rotate
 
A2 = A(Px,Px);
 figure;
 plot_binary_duo(A2,1,1);
 plot_cluster_grid(nx, nx, Nx, Nx);
 
 %save(strcat('A_',name,'_rand.mat'),'A')
 %save(strcat('F_',name,'_rand.mat'),'F')
  
 break;
 end
end

  

end