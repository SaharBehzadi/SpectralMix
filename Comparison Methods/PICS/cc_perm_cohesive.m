function [Px, Py] = cc_perm_cohesive (k, l, Qx, Qy, F)

Py = cc_label_perm(l, Qy);

f = size(F,2);


Px = [];
for i=1:k
    ind = find(Qx==i);
    ind = ind';
    for j=1:f
       ix = find(F(ind,j)==1);
        Px = [Px; ind(ix)];
    end
    ix = find(sum(F(ind,:),2)==0);
    Px = [Px; ind(ix)];
end



end