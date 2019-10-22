function plot_trace_bin(A, F, tracename, isSelfGraph)

%figure; plot_binary_duo(A,1); title(strcat('nnz = ',num2str(sum(sum(A)))))
%figure; plot_binary_duo(F,0); title(strcat('nnz = ',num2str(sum(sum(F)))))
%xlabels = {'liberal','conservative'};
%xlabels = {'male','female','UK-gender','young','middle','older','UK-age'};



%[A F] = buildperfect(A,F);


[nx ny] = size(F)

trace={};
load(tracename)
length(trace)

for i=length(trace):-1:1
    tr = trace{i};
    [tr.k tr.l]
%     if(tr.k>1 && tr.l > 1)
%         tr.Nx
%         cumsum( tr.Nx)
%         tr.Ny
%     show_layout_shaded(tr.Qx,tr.Qx,tr.Dnz,1)
%      %plot_cluster_grid (nx, nx, tr.Nx, tr.Nx)
%      pause
%     show_layout_shaded(tr.Qx,tr.Qy,tr.DnzF,0)
%     end

    
     [tr.k tr.l]
     tr.Nx 
     tr.Ny
     pause
      [nx ny] = size(F);
     [Px, Py] =  cc_perm(tr.k, tr.l, tr.Qx, tr.Qy);

    if(isSelfGraph)
         figure;
         plot_binary_duo(F(Px,:),0);
          for x = cumsum(tr.Nx)
              if x < nx
               
               line([0.5, ny+0.5], [x+0.5 x+0.5],'LineStyle','--');
              end
          end
         %set(gca,'XTick',1:7,'XTicklabel',xlabels)
    else
         figure;
         plot_binary_duo(F(Px,Py),1,2);
         plot_cluster_grid(nx, ny, tr.Nx, tr.Ny);
         %set(gca,'XTick',1:7,'XTicklabel',xlabels(Py))
%          if(tr.k>1 && tr.l>1)
%          show_layout_shaded(tr.Qx,tr.Qy,tr.DnzF,0)
%          end
    end
    
    
     figure;
     plot_binary_duo(A(Px,Px),1,1);
     plot_cluster_grid(nx, nx, tr.Nx, tr.Nx);
    % if(tr.k>1)
     %show_layout_shaded(tr.Qx,tr.Qx,tr.Dnz,1)
    % end
     
    pause
    break;
end

end