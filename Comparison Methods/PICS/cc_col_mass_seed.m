function Qy = cc_col_mass_seed (A, l)
% CC_COL_MASS_SEED  seed initial column label map (i.e.,
%   clustering) based on equal-mass distribution among clusters)
%
% A  binary matrix
% l  number of clusters
%
% $Id: cc_col_mass_seed.m,v 1.4 2004/02/11 03:00:51 deepay Exp $

[nx ny] = size(A);

counts = sum(A);  % Could be done much faster with access to internals...
[sorted_counts, sort_index] = sort(counts);

if 2*sorted_counts(1) <= l
  error('Too many clusters!');
end

Qy = zeros(1, ny);

total_count = nnz(A);    % equal to sum(counts)
thresh_step = total_count/l;

thresh = thresh_step; % start new cluster when cum_count > thresh
cum_count = 0; % total count so far
cur_jl = 1; % current cluster number
% Access in ascending-count order
for j = sort_index
  Qy(j) = cur_jl;
  cum_count = cum_count + counts(j);
  if cum_count >= thresh
    %cum_count, thresh, cur_jl  % DEBUG
    thresh = thresh + thresh_step;
    cur_jl = cur_jl + 1;
  end
  % BEGIN_DEBUG
  %if cur_jl > l
  %  disp(sprintf('OOPS: j %d (%d) / cur_jl %d / %f (%f) %d', j, ...
  %               counts(j), cur_jl, thresh-thresh_step, ...
  %               thresh_step, find(sort_index==j)));
  %end
  % END_DEBUG
end

%assert(cur_jl == l+1, 'Wrong number of clusters assigned!');
%assert(Qy(sort_index(end)) == l, 'Wrong number of clusters!');

% $Log: cc_col_mass_seed.m,v $
% Revision 1.4  2004/02/11 03:00:51  deepay
% Removing assert and div statements.
%
% Revision 1.3  2004/02/04 23:56:17  spapadim
% Misc touchups
%

