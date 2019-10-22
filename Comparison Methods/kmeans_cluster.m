function [centroids,clusters,err,qerrs]=kmeans_cluster(method,D,n,epochs,verbose)
% KMEANS The k-means clustering method. 
%
% [centroids,clusters,err,qerrs] = kmeans_cluster(method,D,n,[epochs],[verbose])
%
%   C = kmeans_cluster('batch',D,10,20);
%   [C,partition] = kmeans_cluster('seq',D,D(1:10,:),20,1);
%
%  Input and output arguments ([]'s are optional):
%   method    (string) 'batch' or 'seq'
%   D         (matrix) size dlen x dim, the data
%             (struct) map or data struct
%   n         (scalar) number of centroids
%             (matrix) size n x dim, initial values for the centroids
%   [epochs]  (scalar) number of training epochs, 100 by default
%   [verbose] (scalar) verbose level, 0 by default
%
%   centroids (matrix) size n x dim, the k-means centroids
%   clusters  (vector) size dlen x 1, cluster number for each sample
%   err       (scalar) total quantization error for the data set
%
% See also KMEANS_CLUSTERS.

% References: 
%   Jain, A.K., Dubes, R.C., "Algorithms for Clustering Data", 
%   Prentice Hall, 1988, pp. 96-101.
%
%   Moody, J., Darken C.J., "Fast Learning in Networks of
%   Locally-Tuned Processing Units", Neural Computation, vol 1.,
%   no. 2, 1989, pp. 281-294. 

% Contributed to SOM Toolbox vs2, February 2nd, 2000 by Esa Alhoniemi
% Copyright (c) by Esa Alhoniemi
% http://www.cis.hut.fi/projects/somtoolbox/

error(nargchk(3,5,nargin))
data = D; 

[dlen,dim] = size(data);
n = min(dlen,n);

if prod(size(n))==1, 
  temp = randperm(dlen);
  centroids = data(temp(1:n),:);
else
  centroids = n;
  n = size(centroids,1);
end

if (nargin<4) epochs = 100; end
if (nargin<5) verbose = 0; end

lr = 0.5; % initial learning rate for sequential k-means
clusters  = zeros(1,dlen);       

switch method
 case 'seq'
  len = epochs*dlen;
  l_rate = linspace(lr,0,len);
  order  = randperm(dlen);
  for iter = 1:len
    x = data(order(rem(iter,dlen)+1),:); % pick one sample vector
    dx = x(ones(n,1),:)-centroids; % difference 
    [dist,nearest] = min(sum(dx.^2,2)); % find nearest centroid
    centroids(nearest,:) = centroids(nearest,:)+l_rate(iter)*dx(nearest,:);
  end 
 case 'batch'
  for iter=1:epochs,    
    old_clusters = clusters;
    [dummy,clusters] = min((ones(dlen,1)*sum((centroids.^2)',1)-2.*(data*(centroids')))');
    for i = 1:n
      f = find(clusters==i);
      s = length(f);
      if (~isempty(f)) centroids(i,:) = sum(data(f,:),1)/s; end
    end    
    if ((iter>0)&&all(old_clusters==clusters))
      if (verbose) fprintf('Convergence in %d iterations\n',iter); end
      break
    end    
  end  
 otherwise,
  fprintf(2, 'Unknown method\n');
end

[qerrs,clusters] = min(((ones(n,1)*sum((data.^2)',1))'+ones(dlen,1)*sum((centroids.^2)',1)-2.*(data*(centroids')))');
err = sum(qerrs);

if (size(clusters,1)==1) clusters=clusters'; end

return
