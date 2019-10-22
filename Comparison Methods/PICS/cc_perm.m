function [Px, Py] = cc_perm (k, l, Qx, Qy)
% CC_PERM  return the row and column permutations, according to the
% row and column label maps; this is just a shorthand utility function
%
% $Id: cc_perm.m,v 1.2 2004/02/04 23:56:17 spapadim Exp $

Px = cc_label_perm(k, Qx);
Py = cc_label_perm(l, Qy);

% $Log: cc_perm.m,v $
% Revision 1.2  2004/02/04 23:56:17  spapadim
% Misc touchups
%

