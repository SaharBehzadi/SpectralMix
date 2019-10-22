function [Dnz, DnzF, Nx, Ny, Qx, Qy] = ...
    cc_iter_duo(A, F, k, l, Dnz0, DnzF0, Nx0, Ny0, Qx0, Qy0)

% call as (from cc_duo --shuffle rows): [Dnz1, DnzF1, Ny1, Nx1, Qy1, Qx1] = cc_iter_duo(At, Ft, l, k, Dnz', DnzF1', Ny1, Nx1, Qy1, Qx1);
% CC_ITER  perform one iteration (over columns) for co-clustering
%
% $Id: cc_iter.m,v 1.10 2004/04/30 14:00:34 deepay Exp $

% if(nargin < 9)
%   isSelfGraph = false;
% end
% 
% if(nargin < 10 && isSelfGraph)
At = A';
% end

[nx, ny] = size(F);

% compute number of zeros
Nyy0 = Ny0' * Ny0;  % outer-product, but these are row-vectors
%assert(size(Nxy0) == [k l], 'OOPS');
Dz0 = Nyy0 - Dnz0;

Nxy0 = Nx0' * Ny0;  % outer-product, but these are row-vectors
%assert(size(Nxy0) == [k l], 'OOPS');
DzF0 = Nxy0 - DnzF0;

Nx = Nx0;
Ny = Ny0;
%Nxy = Nxy0;
%%%%
%Nxx = Nxx0;
%%%%
%Qx = zeros(1,nx);
%Qy = zeros(1,ny);
%Dnz = Dnz0;
%%%%
%DnzF = DnzF0;
%%%%

% Pre-compute transpose for time s(h)avings
Nx_t = Nx';
Ny_t = Ny';

% Pre-compute the logarithm of probabilities
% for A
Pnz0 = Dnz0./Nyy0; Pnz0(~isfinite(Pnz0)) = 0;
Pz0 = Dz0./Nyy0; Pz0(~isfinite(Pz0)) = 0;
Lnz0 = entropy_bits(Pnz0);
Lz0 = entropy_bits(Pz0);
%if(isSelfGraph)
  Lnz0t = Lnz0';
  Lz0t = Lz0';
%end 

%for F
PnzF0 = DnzF0./Nxy0; PnzF0(~isfinite(PnzF0)) = 0;
PzF0 = DzF0./Nxy0; PzF0(~isfinite(PzF0)) = 0;
LnzF0 = entropy_bits(PnzF0);
LzF0 = entropy_bits(PzF0);

% BEGIN_DEBUG
%Nx0, Ny0
%Dnz0, Nxy0, Pnz0, Lnz0
%Dz0, Nxy0, Pz0, Lz0
% END_DEBUG


Qy = Qy0;
Qx = Qx0;
for j = 1:ny
      C = A(:,j); %splice A into k parts : BUT now l (order change)
      Cnz = cc_col_nz(C, l, Qy0);  % Number of non-zeros, per cluster
      Cz = Ny_t - Cnz;
      AllEntropyTerms = Cz' * Lz0 + Cnz' * Lnz0;
      %if (isSelfGraph)
        C1 = At(:,j);
        C1nz = cc_col_nz(C1, l, Qy0);
        C1z = Ny_t - C1nz;
        AllEntropyTerms1 = C1z' * Lz0t + C1nz' * Lnz0t;
        DiagEntry = A(j,j);
      %end

      CF = F(:,j); %splice F into l parts : BUT now k (order change)
      CnzF = cc_col_nz(CF, k, Qx0);  % Number of non-zeros, per cluster
      CzF = Nx_t - CnzF;
      AllEntropyTermsF = CzF' * LzF0 + CnzF' * LnzF0;


      jl_min = 0;
      jl_min_cost = +Inf;
      % For each candidate column cluster, compute "cost"
      for jl = 1:l
        % Also tried pre-extracting slices to cell array, but makes
        % things slower(!?) instead...
        %entropy_terms = Cz .* Lz0(:,jl) + Cnz .* Lnz0(:,jl);
        %entropy_terms = AllEntropyTerms(jl);
        cost_jl = AllEntropyTerms(jl) + AllEntropyTermsF(jl);
        %cost_jl = cc_total(entropy_terms);
        %notSelfGraphCost = cost_jl;

        %if (isSelfGraph)
          %entropy_terms = C1z .* Lz0t(:,jl) + C1nz .* Lnz0t(:,jl);
          %entropy_terms = AllEntropyTerms1(jl);
          tmp1 = AllEntropyTerms1(jl);
          %tmp1 = cc_total(entropy_terms);
          cost_jl = cost_jl + tmp1;


          if(DiagEntry == 1)
            cost_jl = cost_jl - Lnz0(Qy0(j),jl) - Lnz0(jl,Qy0(j)) + Lnz0(jl,jl);
            %cost_jl = cost_jl - Lnz0(Qx0(j),Qx0(j)) + Lnz0(jl,jl);
          else
            cost_jl = cost_jl - Lz0(Qy0(j),jl) - Lz0(jl,Qy0(j)) + Lz0(jl,jl);
            %cost_jl = cost_jl - Lz0(Qx0(j),Qx0(j)) + Lz0(jl,jl);
          end

        %end

        if (cost_jl < jl_min_cost) % || (cost_jl<jl_min_cost+0.001 && jl<jl_min)) 
          jl_min = jl;
          jl_min_cost = cost_jl;
        end


        %disp(sprintf('j %d / jl %d / notSelfGraphCost %f / tmp_cost_jl %f / cost %f', j, jl, notSelfGraphCost, tmp_cost_jl, cost_jl));
        %debug_disp(sprintf('  Cz [%s] Lz0 [%s]', num2str(Cz'), num2str(Lz0(:,jl)')));
        %debug_disp(sprintf('  Cnz [%s] Lnz0 [%s]', num2str(Cnz'), num2str(Lnz0(:,jl)')));
      end


      % BEGIN_DEBUG
      %if Qy(j) ~= Qy0(j)
      %  debug_disp(sprintf(' Reassigning %d: %d -> %d (cost: %f)', ...
      %                     j, Qy0(j), Qy(j), jl_min_cost));
      %end
      % END_DEBUG


      % Set new cluster label
      %disp(sprintf('==> Choose %d',jl_min));
      Qy(j) = jl_min;
      %if(isSelfGraph)
      %  Qx(j) = jl_min; --> Qx, (in real Qy), holds col shuffles 
      % Only Qx (here Qy, above) holds row shuffles
      %end


      % Update non-zero counts
    %   if(isSelfGraph == false)
    %     %% In this case, incremental updates work.
    %     Dnz(:,Qy0(j)) = Dnz(:,Qy0(j)) - Cnz;
    %     Dnz(:,Qy(j)) = Dnz(:,Qy(j)) + Cnz;
    %   end

      % Update size
      Ny(Qy0(j)) = Ny(Qy0(j)) - 1;
      Ny(Qy(j)) = Ny(Qy(j)) + 1;
    %   if(isSelfGraph)
    %     Nx(Qx0(j)) = Nx(Qx0(j)) - 1;
    %     Nx(Qx(j)) = Nx(Qx(j)) + 1;
    %   end  

end

%if(isSelfGraph)
  %% If this isn't a self-graph, the incremental updates have been applied
  %% already. If it *is*, we must do the whole recomputation.
%   Dnz = zeros(k,l);
%   for jl = 1:l
%     Aslice = A(:,find(Qy==jl));
%     for ik = 1:k
%       Dnz(ik,jl) = nnz(Aslice(find(Qx==ik),:));
%     end
%   end
%end
  Dnz = zeros(l,l);
  for jl = 1:l
    Aslice = A(:,find(Qy==jl));
    for il = 1:l
      Dnz(il,jl) = nnz(Aslice(find(Qy==il),:));
    end
  end
  
  DnzF = zeros(k,l);
  for jl = 1:l
    Fslice = F(:,find(Qy==jl));
    for ik = 1:k
      DnzF(ik,jl) = nnz(Fslice(find(Qx==ik),:));
    end
  end 



% $Log: cc_iter.m,v $
% Revision 1.10  2004/04/30 14:00:34  deepay
% *** empty log message ***
%
% Revision 1.9  2004/04/28 19:07:20  deepay
% *** empty log message ***
%
% Revision 1.8  2004/04/12 17:21:59  deepay
% Little stuff
%
% Revision 1.7  2004/04/10 23:49:42  deepay
% *** empty log message ***
%
% Revision 1.6  2004/04/02 15:14:21  deepay
% Should mostly work for self-graphs
%
% Revision 1.5  2004/02/21 03:03:56  spapadim
% Very minor tweak to shave of dt
%
% Revision 1.4  2004/02/20 10:27:57  spapadim
% Remove call to eltsum() and use sum(sum()) instead, for efficiency.  Also,
% use the MEX version of cc_col_nz for greatest improvement.
%
% Revision 1.3  2004/02/11 03:00:51  deepay
% Removing assert and div statements.
%
% Revision 1.2  2004/02/04 09:43:46  spapadim
% Misc fixes/checks
%
