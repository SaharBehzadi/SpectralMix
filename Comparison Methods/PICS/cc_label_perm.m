function Px = cc_label_perm (k, Qx)
% CC_LABEL_PERM  return permutation vector, based on label map
%   so that the elements are sorted by class label
%
% $Id: cc_label_perm.m,v 1.2 2004/02/04 23:56:17 spapadim Exp $

Px = zeros(size(Qx));

start = 1;
for ik = 1:k
  idxs_k = find(Qx==ik);
  len_k = length(idxs_k);
  Px(start:(start+len_k-1)) = idxs_k;
  start = start + len_k;
end

% $Log: cc_label_perm.m,v $
% Revision 1.2  2004/02/04 23:56:17  spapadim
% Misc touchups
%

