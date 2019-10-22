function plot_trace(A, F)

%figure; plot_binary_duo(A,1); title(strcat('nnz = ',num2str(sum(sum(A)))))
%figure; plot_binary_duo(F,0); title(strcat('nnz = ',num2str(sum(sum(F)))))

disp 'plot done'
pause

[nx ny] = size(F)

trace={};
load('trace_youtube.mat')
length(tr)

for i=length(tr):-1:1
    %tr = tr{i};
    [tr.k tr.l]
    if(tr.k>1 && tr.l > 1)
        tr.Nx
        cumsum( tr.Nx)
        tr.Ny
    show_layout_shaded(tr.Qx,tr.Qx,tr.Dnz,1)
     %plot_cluster_grid (nx, nx, tr.Nx, tr.Nx)
     pause
    show_layout_shaded(tr.Qx,tr.Qy,tr.DnzF,0)
    end

    
%      [tr.k tr.l]
%      tr.Nx 
%      tr.Ny
%      pause
%       [nx ny] = size(F);
%      [Px, Py] =  cc_perm(tr.k, tr.l, tr.Qx, tr.Qy);
% 
% 
%      figure;
%      plot_binary_duo(F(Px,Py),1);
%      plot_cluster_grid(nx, ny, tr.Nx, tr.Ny);
%     
%      figure;
%      plot_binary_duo(A(Px,Px),1);
%      plot_cluster_grid(nx, nx, tr.Nx, tr.Nx);
    
    
    pause
end

end