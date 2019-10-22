function density = computeDensity(A, tracefiles, tit)

A = triu(A,1)+tril(A,-1);

E = sum(sum(A))

density = zeros(2,1);

for file=1:2
trace={};
load(tracefiles{file})
length(trace)

for i=length(trace):-1:1
    tr = trace{i};
    [tr.k tr.l]
   
    tr.Nx
    
     [Px, Py] = cc_perm(tr.k, tr.l, tr.Qx, tr.Qy);
     figure; plot_binary_duo(A(Px,Px),1);
     plot_cluster_grid(size(A,1), size(A,1), tr.Nx, tr.Nx);
     
     
    EinC = 0;
    for c=1:tr.k
        ind = find(tr.Qx == c);
        assert(length(ind) == tr.Nx(c));
        sum(sum(A(ind, ind)))
        EinC = EinC + sum(sum(A(ind, ind)));
    end
    
    break;
end

density(file) = EinC/E
pause
end

figure; hold all;
bar(1, density(1),'b')
bar(2, density(2),'r')
set(gca, 'XTick',[1 2],'XTickLabel',{'CC','CC-DUAL'})
set(gca, 'FontSize',16)
title(tit)
ylabel('Density')

pause
close all

end