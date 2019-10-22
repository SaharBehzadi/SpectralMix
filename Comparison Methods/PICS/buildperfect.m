function [Atemp Ftemp] = buildperfect(A,F)

xlabels = {'conservative','liberal'}; %  'neutral'
%%%%%%%%%%%%%%%%%%%%
% compute cost if all groups were perfectly homogeneous
perfind = [];
Nx = [];
for j=1:size(F,2)
    ind = find(F(:,j) == 1);
    perfind = [perfind; ind];
    Nx = [Nx length(ind)];
end
Atemp = A(perfind,perfind);
Ftemp = F(perfind,:);

[nx ny] = size(F);

 figure;
 plot_binary_duo(Ftemp,1);
 plot_cluster_grid(nx, ny, Nx, [1 1 1]);
 set(gca,'XTick',[1 2 ],'XTicklabel',xlabels)
  figure;
 plot_binary_duo(Atemp,1);
 plot_cluster_grid(nx, nx, Nx, Nx);
 
 pause
 cumNx = cumsum(Nx);
% compute nnz
Dnztemp = zeros(size(F,2),size(F,2));
for i=1:size(F,2)
    for j=1:size(F,2)
        if(i>1)
            range1 = cumNx(i-1)+1:cumNx(i);
        else
            range1 = 1:cumNx(i);
        end
        if(j>1)
            range2 = cumNx(j-1)+1:cumNx(j);
        else
            range2 = 1:cumNx(j);
        end
       
        Dnztemp(i,j) = sum(sum(Atemp( range1,range2 )));
    end
end

DnzFtemp = zeros(size(F,2),size(F,2));
for i=1:size(F,2)
    for j=1:size(F,2)
        if(i>1)
            range1 = cumNx(i-1)+1:cumNx(i);
        else
            range1 = 1:cumNx(i);
        end
        DnzFtemp(i,j) = sum(sum(Ftemp( range1,j )));
    end
end

Nx
Dnztemp
DnzFtemp
[costt, costt2] = cc_cost_duo(size(F,2),size(F,2), Nx, ones(1,size(F,2)), Dnztemp, DnzFtemp);
disp(sprintf('## Perfect homogeneous cost %f (C2: %f)', costt, costt2));
pause
%%%%%%%%%%%%%%%%%%%%

end