function nz = cc_col_nz (C, k, Qx)
% CC_COL_NZ  return k-by-1 row vector of non-zeros per row-cluster
%
% $Id: cc_col_nz_supersededByMex.m,v 1.1 2004/04/28 19:10:18 deepay Exp $

%tmp = find(C);
%length(Qx)
%max(tmp)
idxs = Qx(find(C));
if isempty(idxs)
  nz = zeros(k,1);
else
  nz = histc(idxs, 1:k)';
end

% old, defunct code
% nz = zeros(k,1);
% for ik = 1:k
%  nz(ik) = nz(ik) + sum(idxs==ik);
% end

% $Log: cc_col_nz_supersededByMex.m,v $
% Revision 1.1  2004/04/28 19:10:18  deepay
% *** empty log message ***
%
% Revision 1.5  2004/04/02 15:14:21  deepay
% Should mostly work for self-graphs
%
% Revision 1.4  2004/02/14 03:55:13  deepay
% Whatever
%
% Revision 1.3  2004/02/04 23:56:17  spapadim
% Misc touchups
%

