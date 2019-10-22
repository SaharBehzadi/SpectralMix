function [Px,Py,Nx,Ny,Dnz1] = layout(A,k,l,Qx,Qy,Dnz)
% Organize clusters so that the biggest ones are towards the top-left.
% Organize rows and columns within clusters so that densest ones are
% topmost/leftmost.

Py = zeros(size(Qy));
Px = zeros(size(Qx));
Nx = histc(Qx,1:k);
Ny = histc(Qy,1:l);

Dnz1 = zeros(size(Dnz));
start = 1;


rowsdone = [];
colsdone = [];

while start <= k & start <= l
  [t1,t2] = sort(Dnz);
  [t3,t4] = sort(t1(k,:));
  col = t4(length(t4));
  row = t2(k,col);

disp(sprintf('%d %d: %d, %d:  %d',row,col,Nx(row),Ny(col),Dnz(row,col)));

  Qrows = find(Qx==row);
  Qcols = find(Qy==col);
  if ~isempty(Qrows)
    Px(Qrows) = start;
    Nx(start) = length(Qrows);
    rowsdone = [rowsdone row];
  end
  if ~isempty(Qcols)
    Py(Qcols) = start;
    Ny(start) = length(Qcols);
    colsdone = [colsdone col];
  end
  
  Dnz(row,:) = -1;
  Dnz(:,col) = -1;

  start = start+1;
end


rowstodo = setdiff([1:k],rowsdone);
if ~isempty(rowstodo)
  start1 = start;
  for tmp=1:length(rowstodo)
    Qrows = find(Qx==rowstodo(tmp));
    if ~isempty(Qrows)
      Px(Qrows) = start1;
      Nx(start1) = length(Qrows);
      start1 = start1+1;
    end
  end
end

colstodo = setdiff([1:l],colsdone);
if ~isempty(colstodo)
  start1 = start;
  for tmp=1:length(colstodo)
    Qcols = find(Qy==colstodo(tmp));
    if ~isempty(Qcols)
      Py(Qcols) = start1;
      Ny(start1) = length(Qcols);
      start1 = start1+1;
    end
  end
end

Dnz1 = zeros(k,l);
for jl = 1:l
  Aslice = A(:,find(Py==jl));
  for ik = 1:k
    Dnz1(ik,jl) = nnz(Aslice(find(Px==ik),:));
  end
end


%[PX,PY] = cc_perm(k,l,Px,Py);

%[nx,ny] = size(A);
%figure;
%plot_binary(A(PX,PY));
%plot_cluster_grid(nx, ny, Nx, Ny);
%title('Clustered matrix');

%show_layout_shaded(Px,Py,Dnz1,isSelfGraph);
