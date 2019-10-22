function [l, Ny, Qy] = add_cluster_hellscream(l0, Ny0, Qy0, k0, Nx0, Qx0, Dnz0, A, isSelfGraph, At)
% To add a new cluster: 
%   find column cluster with max entropy per column
%   switch all columns whose removal lessens the entropy per column
%
% NOTE: This function requires heavy doses of black magic for proper
%       functioning. Run it only during "... those hours of darkness
%       when the powers of evil are exalted".

warning off MATLAB:divideByZero;

EPSILON = 1e-05;



if(nargin<10 && isSelfGraph)
	At = A';
end

% Compute the entropies, a la Spiros
Nxy0 = Nx0' * Ny0;
Dz0 = Nxy0 - Dnz0;
Pz = Dz0./Nxy0; Pz(~isfinite(Pz)) = 0;
Pnz = Dnz0./Nxy0; Pnz(~isfinite(Pnz)) = 0;
%entropy_terms =  Dz0 .* my_ent_bits(Pz) + Dnz0 .* my_ent_bits(Pnz);
entropy_terms =  Dz0 .* entropy_bits(Pz) + Dnz0 .* entropy_bits(Pnz);

% Entropy per col
entropy_for_clusters = (sum(entropy_terms,1))./Ny0;
if(isSelfGraph)
  % We should also consider the entropy of the corresponding row clusters	
  entropy_for_clusters = entropy_for_clusters + (sum(entropy_terms,2))'./Nx0;	
end  
[Ent, Ind] = sort(entropy_for_clusters);
max_entropy_per_col = Ent(length(Ent));
max_entropy_cluster = Ind(length(Ind));

% Find the columns of the most painful cluster
PC = find(Qy0==max_entropy_cluster);

Qy = Qy0;
l = l0;
Ny = Ny0;  % No new cluster added

DnzMax = Dnz0(:,max_entropy_cluster)';
DzMax = Dz0(:,max_entropy_cluster)';
if(isSelfGraph)
  Dnz2Max = Dnz0(max_entropy_cluster,:);
  Dz2Max = Dz0(max_entropy_cluster,:);
end

Nxy1 = Nx0 .* (Ny0(max_entropy_cluster)-1);
numCols = Ny0(max_entropy_cluster)-1;
if(isSelfGraph)
  Nxy2 = (Nx0(max_entropy_cluster)-1) .* Ny0;
  numRows = Nx0(max_entropy_cluster)-1;
end


for i=1:length(PC)


  % How much do I gain if I send column PC(i) to the new cluster
  currentCol = A(:,PC(i));
  nnzByRow = cc_col_nz(currentCol, k0, Qx0)';
  
  %tmp = Qx0(find(currentCol));
  %if isempty(tmp)
  %  nnzByRow = zeros(1,k0); %% <--- Think it should be (k0,1).
  %else
  %  nnzByRow = histc(Qx0(find(currentCol)),1:k0);
  %end

  if(isSelfGraph)
    currentRow = At(:,PC(i));
    nnzByCol = cc_col_nz(currentRow, k0, Qx0)'; % Rows and cols are equiv.
    %DiagEntry = A(j,j);

    %tmp1 = Qy0(find(currentRow));
    %if isempty(tmp1)
    %  nnzByCol = zeros(1,l0);
    %else
    %  nnzByCol = histc(Qy0(find(currentRow)),1:l0);
    %end
  end

	

    %%disp(sprintf('%d  %d', size(nnzByRow), size(DnzMax)));
  Dnz1 = DnzMax - nnzByRow;
  Dz1 = DzMax - (Nx0-nnzByRow);
  Pz1 = double(Dz1)./Nxy1; Pz1(~isfinite(Pz1)) = 0;
  Pnz1 = double(Dnz1)./Nxy1; Pnz1(~isfinite(Pnz1)) = 0;
%%  new_entropy = double(Dz1) .* my_ent_bits(Pz1) + double(Dnz1) .* my_ent_bits(Pnz1);
  new_entropy = double(Dz1) .* entropy_bits(Pz1) + double(Dnz1) .* entropy_bits(Pnz1);
  new_entropy_per_col = (sum(new_entropy)) ./ numCols;

  if(isSelfGraph)
    Dnz2 = Dnz2Max - nnzByCol;
    Dz2 = Dz2Max - (Ny0-nnzByCol);
    Pz2 = double(Dz2)./Nxy2; Pz2(~isfinite(Pz2)) = 0;
    Pnz2 = double(Dnz2)./Nxy2; Pnz2(~isfinite(Pnz2)) = 0;
    new_entropy_per_col = new_entropy_per_col + ...
      (double(Dz2) .* entropy_bits(Pz2) + double(Dnz2) .* entropy_bits(Pnz2))./ numRows;
  end


    
  if(double(new_entropy_per_col) < double(max_entropy_per_col) - EPSILON)

%%%disp(sprintf('i=%d    maxEnt = %7.7f, newEnt=%7.7f',i,max_entropy_per_col,new_entropy_per_col));

    % This column was perhaps a good choice
    if(l==l0)
         l = l0+1;
         Ny = [Ny0 0];
         if(isSelfGraph)
            k = k0+1;
            Nx = [Nx0 0];
        end	
    end
    Ny(l) = Ny(l)+1;
    Ny(max_entropy_cluster) = Ny(max_entropy_cluster)-1;
    Qy(PC(i)) = l;
    if(isSelfGraph)
      Nx(k) = Nx(k)+1;
      Nx(max_entropy_cluster) = Nx(max_entropy_cluster)-1;
      Qx(PC(i)) = k;
    end  

    % Changing DnzMax and DzMax and Nxy1 and max_entropy_per_col
    % and numCols
    DnzMax = Dnz1;
    DzMax = Dz1;
    Nxy1 = Nxy1 - Nx0;
    numCols = numCols - 1;
    if(isSelfGraph)
      Dnz2Max = Dnz2;
      Dz2Max = Dz2;
      Nxy2 = Nxy2 - Ny0;
      numRows = numRows - 1;
    end
    max_entropy_per_col = new_entropy_per_col;
  end



end


