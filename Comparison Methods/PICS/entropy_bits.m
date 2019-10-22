function L = entropy_bits(A)
% ENTROPY_BITS  Return -log2 of A, defining -log2(0) \approx Inf

tiny = exp(-700);
A = A + tiny;
L = -log2(A);

% $Log: entropy_bits.m,v $
% Revision 1.2  2004/02/04 23:56:17  spapadim
% Misc touchups
%

