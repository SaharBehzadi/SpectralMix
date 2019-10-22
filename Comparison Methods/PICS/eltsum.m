function s = eltsum(A)
% ELTSUM  return sum of all elements of A in one fell swoop
%
% $Id: eltsum.m,v 1.2 2004/02/04 23:56:17 spapadim Exp $

s = sum(reshape(A, 1, prod(size(A))));

% $Log: eltsum.m,v $
% Revision 1.2  2004/02/04 23:56:17  spapadim
% Misc touchups
%

