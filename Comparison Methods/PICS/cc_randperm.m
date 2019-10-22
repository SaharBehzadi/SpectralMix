function Qx = cc_randperm (Qx0, f)
% CC_RANDPERM  randomly choose a fraction f of elements in Qx0 and
%    permute these randomly
%
% $Id: cc_randperm.m,v 1.2 2004/02/04 23:56:17 spapadim Exp $

if (f <= 0.0) || (f >= 1)
  error('Fraction not in (0.0, 1.0)');
end

n = length(Qx0);
Qx = Qx0;

part_n = ceil(f*n);

rnd_idxs = unique_random(part_n, n);
rnd_perm = randperm(part_n);

Qx(rnd_idxs) = Qx(rnd_idxs(rnd_perm));

% $Log: cc_randperm.m,v $
% Revision 1.2  2004/02/04 23:56:17  spapadim
% Misc touchups
%

