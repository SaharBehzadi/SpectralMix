function l = logstar2(n)
% LOGSTAR2  return log-star (universal integer code length) of n
%
% $Id: logstar2.m,v 1.2 2004/02/04 23:56:17 spapadim Exp $

l = 0;
while n > 1
  l = l + 1;
  n = log2(n);
end

% $Log: logstar2.m,v $
% Revision 1.2  2004/02/04 23:56:17  spapadim
% Misc touchups
%

