function entropy = computeEntropy(F, tracefiles, tit)


entropy = zeros(2,1);
[nx ny] = size(F);
for file=1:2
% trace={};
% load(tracefiles{file})
% length(trace)
load(tracefiles)
for i=length(tr):-1:1
    %tr = trace{i};
    [tr.k tr.l]
   
    Nx=tr.Nx
    
%      if(file==1)
%              [Px, Py] = cc_perm(tr.k, tr.l, tr.Qx, tr.Qy);
%              figure; 
%              plot_binary_duo(F(Px,:),0);
%              for x = cumsum(Nx)
%               if x < nx
%                line([0.5, ny+0.5], [x+0.5 x+0.5],'LineStyle','--');
%               end
%              end
%      else
%             [Px, Py] = cc_perm(tr.k, tr.l, tr.Qx, tr.Qy);
%             figure;
%             plot_binary_duo(F(Px,Py),1);
%             plot_cluster_grid(nx, ny, Nx, tr.Ny);         
%      end
     
    outersum = 0; 
    for f=1:size(F,2)
        innersum = 0;
        for c=1:tr.k
            ind = find(tr.Qx == c);
            assert(length(ind) == tr.Nx(c));
            % compute entropy
            fv = F(ind,f);
            fv1 = sum(fv);
            pfc1 = fv1/Nx(c);
            pfc0 = 1 - pfc1;
            if(pfc0 == 0 || pfc0 ==1)
                ent=0;
            else
                ent =  -1*(pfc1* log2(pfc1) + pfc0 * log2(pfc0));
            end
            assert(ent<=1);
            assert(ent>=0);
            innersum = innersum + ent*(Nx(c)/sum(Nx));
        end
        assert(innersum<=1);
        outersum = outersum+innersum;
    end
    
    break;
end

entropy(file) = outersum / size(F,2)
pause
end

entropy

figure; hold all;
bar(1, entropy(1),'b')
bar(2, entropy(2),'r')
set(gca, 'XTick',[1 2],'XTickLabel',{'CC','CC-DUAL'})
set(gca, 'FontSize',16)
title(tit)
ylabel('Entropy')

pause
close all

end