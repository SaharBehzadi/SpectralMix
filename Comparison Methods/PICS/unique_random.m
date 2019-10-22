function R = unique_random(N, maximum)
% UNIQUE_RANDOM  return a (possibly sorted) row N-vector of
%   integers between 1 and maximum (inclusive)
%
% $Id: unique_random.m,v 1.2 2004/02/04 23:56:17 spapadim Exp $

R = unique(ceil(maximum*rand(1,N)));
curN = length(R);
while curN < N
  R = unique([R ceil(maximum*rand(1,N-curN))]);
  curN = length(R);
end

% old, defunct code
%rp1 = randperm(maximum);
%rp2 = randperm(maximum);
%R = rp1(rp2(1:N));

% $Log: unique_random.m,v $
% Revision 1.2  2004/02/04 23:56:17  spapadim
% Misc touchups
%

