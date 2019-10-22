function [CA, C2] = contour(A,kmax,lmax,repetitions)
% Find the best clusters for each k and l from 1:kmax and 1:lmax.
% "repetitions" is the number of times the experiment is repeated for
% each setting.

CA = zeros(kmax,lmax);
C2 = zeros(kmax,lmax);
for k=1:kmax
  for l = 1:lmax
    cAmin = Inf;
    c2min = Inf;
    for p = 1:repetitions
      [Nx,Ny,Qx,Qy,Dnz,cA,c2] = cc(A,k,l);
      if cA < cAmin
        cAmin=cA;
      end;
      if c2 < c2min
        c2min=c2;
      end;
    end
    CA(k,l) = cA;
    C2(k,l) = c2;
  end
end
