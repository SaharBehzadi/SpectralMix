function D = celleval(F, C)
% CELLEVAL  evaluate function F on every element of 2-D cell array C

if ndims(C) ~= 2, error('C must be 2-d'), end

[m,n] = size(C)
D = zeros(m,n)
for i = 1:m
  for j = 1:n
    D(i,j) = feval(F, C{i,j});
  end
end

% $Log: celleval.m,v $
% Revision 1.2  2004/02/04 23:56:17  spapadim
% Misc touchups
%

