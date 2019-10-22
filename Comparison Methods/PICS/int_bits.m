function L = int_bits(A)
% INT_BITS  Return log2 of A, defining log2(0) = 0 and rounding up

A = A + (A==0);
L = ceil(log2(A));

% $Log: int_bits.m,v $
% Revision 1.2  2004/02/04 23:56:17  spapadim
% Misc touchups
%

