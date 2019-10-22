function [k, l, Nx, Ny, Qx, Qy, Dnz, cA, c2, CAt, C2t] = cc_search(A, strategy, isSelfGraph,name)
% CC_SEARCH  a naive first implementation
%
% $Id: cc_search.m,v 1.13 2004/04/30 14:00:34 deepay Exp $

cost_threshold = 10^-3;  % if cost change less than this, we stop

SHOW_RESULTS = 0;  % Set this to 0 to stop showing the clusters at each stage

if (nargout > 9) && (nargout < 10)
  error('Must return both cost vectors');
end

if nargout > 9
  CAt = [];  C2t = [];
end

if nargin < 2, strategy = 'aggresive'; end

add_cluster_func = str2func(strcat('add_cluster_', strategy));

[nx, ny] = size(A);
if(isSelfGraph && nx ~= ny)
  if(nx > ny)
    ny = nx;
  else
    nx = ny;
  end
  A(nx,ny) = 0;  %% Makes sure all elements of A are accessible
end

At = A';

k = 1; l = 1;
Dnz = nnz(A);
Nx = nx; Ny = ny;
Qx = ones(1, nx);
Qy = ones(1, ny);

% Get single-block encoding cost
[cost, cost2] = cc_cost(k,l, Nx, Ny, Dnz);
disp(sprintf('## Starting cost %d (C2: %d)', cost, cost2));

trace = {};
trace = add_trace(trace, k, l, Nx, Ny, Qx, Qy, Dnz, cost, cost2);

no_improvements = 0;

if SHOW_RESULTS==1
  figure;
  plot_binary(A,isSelfGraph);
  %title('Original matrix');
end

while 1
  % Add column cluster
  [l1, Ny1, Qy1] = feval(add_cluster_func, l, Ny, Qy, k, Nx, Qx, Dnz, A, isSelfGraph);
  disp(sprintf('#####  Increasing l to %d (k=%d)', l1, k));
  if l1 > l
    if(isSelfGraph)
      [Nx1, Ny1, Qx1, Qy1, Dnz1, new_cost, new_cost2, new_CAt, new_C2t] = ...
        cc(A, l1, l1, isSelfGraph, Qy1, Qy1, At);
    else
      [Nx1, Ny1, Qx1, Qy1, Dnz1, new_cost, new_cost2, new_CAt, new_C2t] = ...
        cc(A, k, l1, isSelfGraph, Qx, Qy1, At);
    end
    
    if cost - new_cost < cost_threshold
      % I am sorry, brother cluster, but you are no good.
      % "Caesar wills that you die!"
      no_improvements = no_improvements + 1;
    else
      no_improvements = 0;
      [l, Dnz, Nx, Ny, Qx, Qy] = deal(l1, Dnz1, Nx1, Ny1, Qx1, Qy1);
      cost = new_cost; cost2 = new_cost2;
      if(isSelfGraph)
	      k = l1;
      end   
      disp(sprintf('## (k=%d) l=%d Current cost %d (C2: %d)',k,l1,new_cost,new_cost2));
      trace = add_trace(trace, k, l, Nx, Ny, Qx, Qy, Dnz, cost, cost2);
      if nargout > 9
        CAt = [CAt new_CAt];
        C2t = [C2t new_C2t];
      end
    end
  else
    no_improvements = no_improvements + 0.5;
  end

  if(isSelfGraph == false)
    % Add row cluster
    [k1, Nx1, Qx1] = feval(add_cluster_func, k, Nx, Qx, l, Ny, Qy, Dnz', At, isSelfGraph);
    disp(sprintf('#####  Increasing k to %d (l=%d)', k1, l));
    if k1 > k
      if(isSelfGraph)
        [Nx1, Ny1, Qx1, Qy1, Dnz1, new_cost, new_cost2, new_CAt, new_C2t] = ...
          cc(A, k1, k1, isSelfGraph, Qx1, Qx1, At);
      else   
        [Nx1, Ny1, Qx1, Qy1, Dnz1, new_cost, new_cost2, new_CAt, new_C2t] = ...
          cc(A, k1, l, isSelfGraph, Qx1, Qy, At);
      end

      if cost - new_cost < cost_threshold
        no_improvements = no_improvements + 0.5;
      else
        no_improvements = 0;
        [k, Dnz, Nx, Ny, Qx, Qy] = deal(k1, Dnz1, Nx1, Ny1, Qx1, Qy1);
        cost = new_cost; cost2 = new_cost2;
        if(isSelfGraph)
  	      l = k1;
        end   
        disp(sprintf('## k=%d (l=%d) Current cost %d (C2: %d)',k1,l, new_cost, new_cost2));
        trace = add_trace(trace, k, l, Nx, Ny, Qx, Qy, Dnz, cost, cost2);
        if nargout > 9
          CAt = [CAt new_CAt];
          C2t = [C2t new_C2t];
        end
      end 
    else
      no_improvements = no_improvements + 0.5;
    end
  end

  if (no_improvements < 1 && SHOW_RESULTS == 1)
    [Px, Py] = cc_perm(k, l, Qx, Qy);

    figure;
    plot_binary(A(Px,Py),isSelfGraph);
    plot_cluster_grid(nx, ny, Nx, Ny);
    %title('Clustered matrix');

    input('continue?','s')
  end

  if no_improvements >= 2, break, end

end
%disp(sprintf('#### Stopping with cost %d (C2: %d)', cost, cost2));

if 0
%if SHOW_RESULTS==0  % show only the final result
  %layout(A,k,l,Qx,Qy,Dnz);

  %show_layout_shaded(Qx,Qy,Dnz,isSelfGraph);

    [Px, Py] = cc_perm(k, l, Qx, Qy);
  
    figure;
    plot_binary(A(Px,Py),isSelfGraph);
    plot_cluster_grid(nx, ny, Nx, Ny);
    %title('Clustered matrix');

end

if nargout > 7, cA = cost; end
if nargout > 8, c2 = cost2; end

disp(sprintf('## Final cost %d (C2: %d)', cA, c2));
save(strcat('trace_cc_',name,'.mat'), 'trace');

return

% -----  Utility functions  -----


function trace = add_trace (trace0, ...
                            k, l, ...
                            Nx, Ny, Qx, Qy, ...
                            Dnz, cA, c2)

trace = trace0;
trace{end+1} = ...
    struct('k', k, 'l', l, ...
           'Nx', Nx, 'Ny', Ny, ...
           'Qx', Qx, 'Qy', Qy, ...
           'Dnz', Dnz, ...
           'cA', cA, 'c2', c2);

return


% $Log: cc_search.m,v $
% Revision 1.13  2004/04/30 14:00:34  deepay
% *** empty log message ***
%
% Revision 1.12  2004/04/28 19:03:23  deepay
% works with new show_layout_shaded
%
% Revision 1.11  2004/04/28 18:59:11  deepay
% working with new plot_binary
%
% Revision 1.10  2004/04/28 18:52:47  deepay
% sundry
%
% Revision 1.9  2004/04/12 17:21:59  deepay
% Little stuff
%
% Revision 1.8  2004/04/02 15:14:21  deepay
% Should mostly work for self-graphs
%
% Revision 1.7  2004/02/19 23:26:14  deepay
% *** empty log message ***
%
% Revision 1.6  2004/02/14 05:50:42  deepay
% Prettification :-(
%
% Revision 1.5  2004/02/14 03:55:13  deepay
% Whatever
%
% Revision 1.4  2004/02/12 22:26:54  deepay
% New heuristic to increase k and l
%
% Revision 1.3  2004/02/06 10:54:48  spapadim
% Minor reorganization
%
% Revision 1.2  2004/02/04 23:56:17  spapadim
% Misc touchups
%
% Revision 1.1  2004/02/04 09:44:34  spapadim
% Initial checkin
%
