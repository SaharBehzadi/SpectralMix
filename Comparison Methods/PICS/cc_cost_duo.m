function [c, c2] = cc_cost_duo(k, l, Nx, Ny, Dnz, DnzF)
% CC_COST  co-cluster code cost
%
% Returned values:
% c   total encoding cost
% c2  C_2 cost (per-block 0/1s only)
%
% $Id: cc_cost.m,v 1.5 2004/02/21 03:03:10 spapadim Exp $

% 0. If there are empty clusters, disregard them
real_k = nnz(Nx);
real_l = nnz(Ny);
if real_k ~= k || real_l ~= l
  non_empty_row_clusters = find(Nx);
  non_empty_col_clusters = find(Ny);
  k = real_k;
  l = real_l;
  Nx = Nx(non_empty_row_clusters);
  Ny = Ny(non_empty_col_clusters);
  Dnz = Dnz(non_empty_row_clusters,non_empty_col_clusters);
end

% 1. Encoding cost for k and l
c = logstar2(k) + logstar2(l);   % FIXME - Use int_bits instead?

% 2. Encoding cost for row cluster sizes
% Nx_bar = cumsum(Nx);
% Nx_bar = Nx_bar(end:-1:1) - k + (1:k);
% c = c + sum(int_bits(Nx_bar));
entrow = Nx .* entropy_bits(Nx/sum(Nx));
c = c + sum(entrow);

% 3. Encoding cost for column cluster sizes
% Ny_bar = cumsum(Ny);
% Ny_bar = Ny_bar(end:-1:1) - l + (1:l);
% c = c + sum(int_bits(Ny_bar));
entcol = Ny .* entropy_bits(Ny/sum(Ny));
c = c + sum(entcol);


% 4. Encoding cost for each block

% 4.0a First compute Nxy, the size of each cluster in A
Nxx = Nx' * Nx;  % outer-product, but these are row vectors
% 4.0b First compute Nxy, the size of each cluster in F
Nxy = Nx' * Ny;  % outer-product, but these are row vectors
% 4.0a Compute number of zeros in blocks of A
Dz = Nxx - Dnz;
% 4.0b Compute number of zeros in blocks of F
DzF = Nxy - DnzF;

% 4.1 Encoding cost for number of non-zeros
%c = c + eltsum(int_bits(Nxx + 1)) + eltsum(int_bits(Nxy + 1));
for i=1:k
    for j=1:k
        c = c + logstar2(Dnz(i,j));
    end
end
for i=1:k
    for j=1:l
        c = c + logstar2(DnzF(i,j));
    end
end

% 4.1a Encoding cost for data in each block of A
Pz = Dz./Nxx; Pz(~isfinite(Pz)) = 0;
Pnz = Dnz./Nxx; Pnz(~isfinite(Pnz)) = 0;
entropy_terms =  Dz .* entropy_bits(Pz) + Dnz .* entropy_bits(Pnz);
c2 = ceil(eltsum(entropy_terms));
fprintf('A-encoding-cost: %f\n', c2)
c = c + c2;
% 4.1b Encoding cost for data in each block of F
Pz = DzF./Nxy; Pz(~isfinite(Pz)) = 0;
Pnz = DnzF./Nxy; Pnz(~isfinite(Pnz)) = 0;
entropy_terms =  DzF .* entropy_bits(Pz) + DnzF .* entropy_bits(Pnz);
c3 = ceil(eltsum(entropy_terms));
fprintf('F-encoding-cost: %f\n', c3)
c = c + c3;

c2 = c2+c3; %cluster encoding cost

% $Log: cc_cost.m,v $
% Revision 1.5  2004/02/21 03:03:10  spapadim
% Check for empty clusters and disregard them.
% **NOTE** This silently changes the "effective" k,l; should not be a problem
% to cc() and cc_search(), but make sure anyway...
%
% Revision 1.4  2004/02/19 21:09:31  spapadim
% Fix what others are scared to touch for some reason :)
% (should not make much of a difference anyway, except for extreme k,l=M,N)
%
% Revision 1.3  2004/02/04 23:56:17  spapadim
% Misc touchups
%

