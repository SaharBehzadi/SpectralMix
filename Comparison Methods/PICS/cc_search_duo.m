function [k, l, Nx, Ny, Qx, Qy, Dnz, DnzF, cA, c2, CAt, C2t] = cc_search_duo(A, F,name)
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

%if nargin < 2, strategy = 'aggresive'; end --do not have
%add_cluster_aggresive.m


% [nx, nx] = size(A);
% if(isSelfGraph && nx ~= ny) % self-graph but not symmetric
%   if(nx > ny)
%     ny = nx;
%   else
%     nx = ny;
%   end
%   A(nx,ny) = 0;  %% Makes sure all elements of A are accessible
% end

%At = A';

%%%%%%%%%%%%%
[nx, ny] = size(F);
%%%%%%%%%%%%%

k = 1; l = 1;
Dnz = nnz(A);
DnzF = nnz(F);
Nx = nx; Ny = ny;
Qx = ones(1, nx);
Qy = ones(1, ny);

% Get single-block encoding cost
[cost, cost2] = cc_cost_duo(k,l, Nx, Ny, Dnz, DnzF);
disp(sprintf('## Starting cost %f (C2: %f)', cost, cost2));


trace = {};
trace = add_trace(trace, k, l, Nx, Ny, Qx, Qy, Dnz, DnzF, cost, cost2);

time = {};

no_improvements = 0;

if SHOW_RESULTS==1
  figure; hold all;
  subplot(2,1,1);
  plot_binary(A,1);
  subplot(2,1,2);
  plot_binary(F,0);
  %title('Original matrix');
end

while 1
  % Add column cluster to F --start easy  
  add_cluster_func = str2func('add_cluster_hellscream');
  to = tic;
  [l1, Ny1, Qy1] = feval(add_cluster_func, l, Ny, Qy, k, Nx, Qx, DnzF, F, 0);   %% ----> OUTER LOOP
  timeoutl = toc(to);
  disp(sprintf('#####  Increasing l to %d (k=%d)', l1, k));
  
  
   %[size(Ny1,2) length(unique(Qy1))]  pause
  if l1 > l
    %if(isSelfGraph)
     % [Nx1, Ny1, Qx1, Qy1, Dnz1, new_cost, new_cost2, new_CAt, new_C2t] = ...
        %cc(A, l1, l1, isSelfGraph, Qy1, Qy1, At);
    %else
      ti = tic;
      [Nx1, Ny1, Qx1, Qy1, Dnz1, Dnz2, new_cost, new_cost2, new_CAt, new_C2t, numiterl] = ... %% ------> INNER LOOP
        cc_duo(A, F, k, l1, Qx, Qy1,0); 
      timeinl = toc(ti);
    %end
    
    if cost - new_cost < cost_threshold
        % I am sorry, brother cluster, but you are no good.
        % "Caesar wills that you die!"
         no_improvements = no_improvements + 0.5;
    else
        no_improvements = 0;
        [l, Dnz, DnzF, Nx, Ny, Qx, Qy] = deal(l1, Dnz1, Dnz2, Nx1, Ny1, Qx1, Qy1);
        cost = new_cost; cost2 = new_cost2;
        %if(isSelfGraph)
        %    k = l1;
        %end   
        disp(sprintf('## (k=%d) l=%d Current cost %d (C2: %d)',k,l1,new_cost,new_cost2));
        trace = add_trace(trace, k, l, Nx, Ny, Qx, Qy, Dnz, DnzF, cost, cost2);
        if nargout > 9
                CAt = [CAt new_CAt];
                C2t = [C2t new_C2t];
        end
    end
  else
    no_improvements = no_improvements + 0.5;
  end

  %if(isSelfGraph == false)
    % Add row cluster
    add_cluster_func = str2func('add_cluster_hellscream_duo');
    to=tic;
    [k1, Nx1, Qx1] = feval(add_cluster_func, k, Nx, Qx, l, Ny, Qy, Dnz', DnzF', A', F');  %%----> OUTER LOOP
    timeoutk = toc(to);
    disp(sprintf('#####  Increasing k to %d (l=%d)', k1, l));
    
    
    if k1 > k
      %if(isSelfGraph)
        %[Nx1, Ny1, Qx1, Qy1, Dnz1, new_cost, new_cost2, new_CAt, new_C2t] = ...
          %cc(A, k1, k1, isSelfGraph, Qx1, Qx1, At);
      %else
	ti=tic;   
        [Nx1, Ny1, Qx1, Qy1, Dnz1, Dnz2, new_cost, new_cost2, new_CAt, new_C2t, numiterk] = ... %% ------> INNER LOOP
            cc_duo(A, F, k1, l, Qx1, Qy,1);
	timeink = toc(ti);
      %end

      if cost - new_cost < cost_threshold
        no_improvements = no_improvements + 0.5;
      else
        no_improvements = 0;
        [k, Dnz, DnzF, Nx, Ny, Qx, Qy] = deal(k1, Dnz1, Dnz2, Nx1, Ny1, Qx1, Qy1);
       
        cost = new_cost; cost2 = new_cost2;
        %if(isSelfGraph)
  	     % l = k1;
        %end   
        disp(sprintf('## k=%d (l=%d) Current cost %d (C2: %d)',k1,l, new_cost, new_cost2));
        trace = add_trace(trace, k, l, Nx, Ny, Qx, Qy, Dnz, DnzF, cost, cost2);
        if nargout > 9
          CAt = [CAt new_CAt];
          C2t = [C2t new_C2t];
        end
      end 
    else
      no_improvements = no_improvements + 0.5;
    end
  %end

  if (no_improvements < 1 && SHOW_RESULTS == 1)
    [Px, Py] = cc_perm(k, l, Qx, Qy);

    figure;
    plot_binary(A(Px,Py),isSelfGraph);
    plot_cluster_grid(nx, ny, Nx, Ny);
    %title('Clustered matrix');

    input('continue?','s')
  end

  time = add_time(time, timeoutl, timeinl, numiterl, timeoutk, timeink, numiterk);

  if no_improvements >= 1, break, end

end % while 1
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
save(strcat(name,'_trace_cc_duo_t.mat'), 'trace');
save(strcat(name,'_time_cc_duo_t.mat'),'time');

return

% -----  Utility functions  -----


function trace = add_trace (trace0, ...
                            k, l, ...
                            Nx, Ny, Qx, Qy, ...
                            Dnz, DnzF, ctot, cclus)

trace = trace0;
trace{end+1} = ...
    struct('k', k, 'l', l, ...
           'Nx', Nx, 'Ny', Ny, ...
           'Qx', Qx, 'Qy', Qy, ...
           'Dnz', Dnz, ...
           'DnzF', DnzF, ...
           'ctot', ctot, 'cclus', cclus);

return
end

function time = add_time (time0, ...
                            timeoutl, timeinl, numiterl, timeoutk, timeink, numiterk)

time = time0;
time{end+1} = ...
    struct('timeoutl', timeoutl, 'timeinl', timeinl, 'numiterl', numiterl, ...
	'timeoutk', timeoutk, 'timeink', timeink, 'numiterk', numiterk	);

return
end


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
end
