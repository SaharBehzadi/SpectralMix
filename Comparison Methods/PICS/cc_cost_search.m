function [CA, C2] = cc_cost_search(A, maxk, maxl, restarts, strategy)
% CC_COST_SEARCH
%
% Inputs:
% A           binary matrix
% maxk, maxl  values of k, l to estimate cost up to
% restarts    number of random restarts to try for each k,l (0 for none)
% strategy    (optional) cluster adding strategy
%
% $Id: cc_cost_search.m,v 1.2 2004/02/12 22:26:54 deepay Exp $

if nargin < 4, restarts = 0; end
if nargin < 5, strategy = 'aggresive'; end

add_cluster_func = str2func(strcat('add_cluster_', strategy));

[nx, ny] = size(A);

CA = repmat(+Inf, maxk, maxl);
C2 = repmat(+Inf, maxk, maxl);

At = A';

Dnz = nnz(A);
Nx = nx; Ny = ny;
Qx = ones(1, nx);
Qy = ones(1, ny);

% Get single-block encoding cost
[cost, cost2] = cc_cost(1, 1, Nx, Ny, Dnz);
disp(sprintf('## Starting cost %d (C2: %d)', cost, cost2));
CA(1,1) = cost;
C2(1,1) = cost2;

[CA, C2] = search_rec_l(CA, C2, maxk, maxl, restarts, ...
                        add_cluster_func, ...
                        A, At, 1, 1, Nx, Ny, Qx, Qy, Dnz);
[CA, C2] = search_rec_k(CA, C2, maxk, maxl, restarts, ...
                        add_cluster_func, ...
                        A, At, 1, 1, Nx, Ny, Qx, Qy, Dnz);

return


% -----  Utility functions  -----

function [CA, C2] = search_rec_l (CA0, C20, maxk, maxl, ...
                                  restarts, add_cluster_func, ...
                                  A, At, k, l, Nx0, Ny0, Qx0, Qy0, Dnz0)

CA = CA0;  C2 = C20;

if l >= maxl, return; end

% Do search for k, l+1
cost = +Inf;
for i = 1:(1+restarts)
  disp(sprintf('### Restart %d (k,l = %d,%d)', i, k, l));
  [l1, Ny1, Qy1] = feval(add_cluster_func, l, Ny0, Qy0, k, Nx0, Qx0, Dnz0, A);
  [Nx1, Ny1, Qx1, Qy1, Dnz1, new_cost, new_cost2] = ...
      cc(A, k, l1, Qx0, Qy1, At);
  if new_cost < cost
    Nx = Nx1; Ny = Ny1;
    Qx = Qx1; Qy = Qy1;
    cost = new_cost; cost2 = new_cost2;
  end
end

if cost < CA(k,l1)
  CA(k,l1) = cost;
  C2(k,l1) = cost2;
end

% Recurse
[CA, C2] = search_rec_l(CA, C2, maxk, maxl, restarts, ...
                        add_cluster_func, ...
                        A, At, k, l1, Nx, Ny, Qx, Qy, Dnz1);
[CA, C2] = search_rec_k(CA, C2, maxk, maxl, restarts, ...
                        add_cluster_func, ...
                        A, At, k, l1, Nx, Ny, Qx, Qy, Dnz1);

return



function [CA, C2] = search_rec_k (CA0, C20, maxk, maxl, ...
                                  restarts, add_cluster_func, ...
                                  A, At, k, l, Nx0, Ny0, Qx0, Qy0, Dnz0)

CA = CA0;  C2 = C20;

if k >= maxk, return; end

% Do search for k+1, l
cost = +Inf;
for i = 1:(1+restarts)
  disp(sprintf('### (Re)start %d (k,l = %d,%d)', i, k, l));
  [k1, Nx1, Qx1] = feval(add_cluster_func, k, Nx0, Qx0, l, Ny0, Qy0, Dnz0', At);
  [Nx1, Ny1, Qx1, Qy1, Dnz1, new_cost, new_cost2] = ...
      cc(A, k1, l, Qx1, Qy0, At);
  if new_cost < cost
    Nx = Nx1; Ny = Ny1;
    Qx = Qx1; Qy = Qy1;
    cost = new_cost; cost2 = new_cost2;
  end
end

if cost < CA(k1,l)
  CA(k1,l) = cost;
  C2(k1,l) = cost2;
end

% Recurse
[CA, C2] = search_rec_l(CA, C2, maxk, maxl, restarts, ...
                        add_cluster_func, ...
                        A, At, k1, l, Nx, Ny, Qx, Qy, Dnz1);
[CA, C2] = search_rec_k(CA, C2, maxk, maxl, restarts, ...
                        add_cluster_func, ...
                        A, At, k1, l, Nx, Ny, Qx, Qy, Dnz1);

return


% $Log: cc_cost_search.m,v $
% Revision 1.2  2004/02/12 22:26:54  deepay
% New heuristic to increase k and l
%
% Revision 1.1  2004/02/06 10:55:50  spapadim
% Preliminary commit; code needs checking\!
%
