function [Nx, Ny, Qx, Qy, Dnz, cA, c2, CAt, C2t] = cc(A, k, l, isSelfGraph, Qx0, Qy0, At, Dnz0)
% CC  co-cluster A into k row-clusters and l column-clusters
%
% Required arguments:
% A        the binary matrix
% k, l     number of row and column clusters
% isSelfGraph	obvious
%
% Optional arguments:
% Qx0, Qy0 starting row and column label maps (i.e., clusterings)
% At       the transpose of A
% Dnz0     pre-computed Dnz
%
% Returned values:
% Nx, Ny   size of each row and column cluster
% Qx, Qy   row and column label maps
% Dnz      number of non-zeros for each D_{ij} matrix
%
% Optional returned values:
% cA       total encoding cost
% c2       C_2 encoding cost (per-block 0/1 values only)
% CAt, C2t vector of costs (as above), per iteration


% $Id: cc.m,v 1.11 2004/04/28 18:59:11 deepay Exp $

COST_THRESHOLD = 10^-3;  % if cost change less than this, we stop

SHOW_RESULTS = 0;  % Set to zero to stop showin results of each iteration

[nx, ny] = size(A);
if(isSelfGraph && nx ~= ny)
  if(nx > ny)
    ny = nx;
  else
    nx = ny;
  end
  A(nx,ny) = 0;  %% Makes sure all elements of A are accessible
end

if (nargout > 7) && (nargout < 9)
  error('Must return both cost vectors');
end

if nargout > 7
  CAt = [];  C2t = [];
end

if nargin > 4
  if nargin < 6
    error('Must provide both Qx0 and Qy0 or neither');
  end
  
  % Assign label maps
  Qx = Qx0;
  Qy = Qy0;
  
  % Initialize cluster sizes
  Nx = histc(Qx, 1:k);
  Ny = histc(Qy, 1:l);
  
else
  % Initialize row cluster sizes
  Nx = zeros(1,k);
  Nx(1:k) = floor(nx/k);
  Nx(k) = Nx(k) + mod(nx,k);

  % Initialize column cluster sizes
  Ny = zeros(1,l);
  Ny(1:l) = floor(ny/l);
  Ny(l) = Ny(l) + mod(ny,l);

  % Compute cluster edges
  Xedges = cumsum([1 Nx]);
  Yedges = cumsum([1 Ny]);
  
  % Initialize row label map
  Qx = zeros(1,nx);
  for ik = 1:k
    Qx(Xedges(ik):(Xedges(ik+1)-1)) = ik;
  end
  
  % Initialize column label map
  Qy = zeros(1,ny);
  for jl = 1:l
    Qy(Yedges(jl):(Yedges(jl+1)-1)) = jl;
  end
end

  if(SHOW_RESULTS)
    [Px1, Py1] = cc_perm(k, l, Qx, Qy);
    figure; plot_binary(A(Px1,Py1),isSelfGraph);
    plot_cluster_grid(nx, ny, Nx, Ny);
    input('continue?','s')
  end

% Initialize non-zero counts
if nargin < 8
  Dnz = zeros(k,l);
  for jl = 1:l
    Aslice = A(:,find(Qy==jl));
    for ik = 1:k
      Dnz(ik,jl) = nnz(Aslice(find(Qx==ik),:));
    end
  end
else
  Dnz = Dnz0;
end


warning off MATLAB:divideByZero;

% Initialize cost
[cost, cost2] = cc_cost(k, l, Nx, Ny, Dnz);

% Pre-compute transpose
% (this is the costliest operation, esp. on sparse matrices)
if nargin < 7
  At = A';
end

% Iterate
iter = 1;
while 1
  %disp(sprintf('*** Iteration %d', iter));
  %disp(sprintf('Starting cost %d (C2: %d)', cost, cost2));
  % Iterate over columns
  [Dnz1, Nx1, Ny1, Qx1, Qy1] = cc_iter(A, k, l, Dnz, Nx, Ny, Qx, Qy, isSelfGraph, At);


  % BEGIN_DEBUG
  %Nx1, Ny1, Dnz1
  %cc_cost(k, l, Nx1, Ny1, Dnz1)

  %title(sprintf('Iter %d (cols) - After', iter));
  % END_DEBUG
  [intermediate_cost, intermediate_cost2] = cc_cost(k, l, Nx1, Ny1, Dnz1);
  disp(sprintf('Intermediate cost %d (C2: %d)', intermediate_cost, intermediate_cost2));
  if nargout > 7
    CAt = [CAt intermediate_cost];
    C2t = [C2t intermediate_cost2];
  end
  
  if(isSelfGraph)
    if(cost - intermediate_cost < COST_THRESHOLD)
      break
    else
      cost = intermediate_cost; cost2 = intermediate_cost2;
    end
  end

  if(SHOW_RESULTS)
    [Px1, Py1] = cc_perm(k, l, Qx1, Qy1);
    figure; plot_binary(A(Px1,Py1),isSelfGraph);
    plot_cluster_grid(nx, ny, Nx1, Ny1);
    input('continue?','s')
  end

  
  if(isSelfGraph == false)
    % Iterate over rows
    [Dnz1, Ny1, Nx1, Qy1, Qx1] = cc_iter(At, l, k, Dnz1', Ny1, Nx1, Qy1, Qx1, isSelfGraph);
    Dnz1 = Dnz1';
    % Estimate new cost
    [new_cost, new_cost2] = cc_cost(k, l, Nx1, Ny1, Dnz1);
    disp(sprintf('New cost %d (C2: %d)', new_cost, new_cost2));
    if cost - new_cost < COST_THRESHOLD, break, end
    cost = new_cost; cost2 = new_cost2;
    %if iter == 1, break, end % DEBUG
    if nargout > 7
      CAt = [CAt new_cost];
      C2t = [C2t new_cost2];
    end
    
    if(SHOW_RESULTS)
      [Px1, Py1] = cc_perm(k, l, Qx1, Qy1);
      figure; plot_binary(A(Px1,Py1),isSelfGraph);
      plot_cluster_grid(nx, ny, Nx1, Ny1);
      input('continue?','s')
    end
  end

  [Dnz, Nx, Ny, Qx, Qy] = deal(Dnz1, Nx1, Ny1, Qx1, Qy1);
  iter = iter + 1;
end

%disp(sprintf('cc(%d,%d) - %d iterations', k, l, iter));

warning on MATLAB:divideByZero;

if nargout > 5, cA = cost; end
if nargout > 6, c2 = cost2; end

% $Log: cc.m,v $
% Revision 1.11  2004/04/28 18:59:11  deepay
% working with new plot_binary
%
% Revision 1.10  2004/04/12 17:21:59  deepay
% Little stuff
%
% Revision 1.9  2004/04/02 15:14:21  deepay
% Should mostly work for self-graphs
%
% Revision 1.8  2004/02/25 06:32:50  spapadim
% Don't remember
%
% Revision 1.7  2004/02/23 04:30:18  deepay
% *** empty log message ***
%
% Revision 1.6  2004/02/14 03:55:13  deepay
% Whatever
%
% Revision 1.5  2004/02/11 03:00:51  deepay
% Removing assert and div statements.
%
% Revision 1.4  2004/02/06 10:54:48  spapadim
% Minor reorganization
%
% Revision 1.3  2004/02/04 23:56:17  spapadim
% Misc touchups
%
