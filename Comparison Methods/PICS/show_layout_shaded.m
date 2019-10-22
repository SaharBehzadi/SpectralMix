function show_layout_shaded(Qx,Qy,Dnz,isSelfGraph)
% Show cluster density using some form of shading

if nargin < 4, isSelfGraph = false; end

[k,l] = size(Dnz);

Nx = histc(Qx, 1:k);
Ny = histc(Qy, 1:l);

Nxy = Nx' * Ny; %% This way, or the other way round?!!! :-(

Pxy = Dnz ./ Nxy;

[Y,I] = sort(Pxy(1:k*l));
%[Y,I] = sort(Dnz(1:k*l));
min = Y(1); max = Y(length(Y));

x = 1;
y = 1;

nx = length(Qx); ny = length(Qy);

figure;
axis([1,length(Qy),1,length(Qx)]);
%set(gca, 'YTick',length(Qy):-1:1)
if(isSelfGraph)
	ylabel('Node Groups','FontSize',14)
        xlabel('Node Groups','FontSize',14)
else
	 ylabel('Node Groups','FontSize',14)
        xlabel('Feature Groups','FontSize',14)
end
%title('Clustered and shaded matrix');

for i=1:k
  y = 1;
  for j=1:l
     fillFraction = Pxy(i,j);
    %fillFraction = Dnz(i,j);
    shading = 1 - ((fillFraction-min)/(max-min));
    patch([y (y+Ny(j)) (y+Ny(j)) y],[(nx-x+1) (nx-x+1) (nx-x-Nx(i)+1) (nx-x-Nx(i)+1)],[shading shading shading]);
    %disp(sprintf('x:%d->%d, y:%d->%d',x,x+Nx(i),y,y+Ny(j)));
    y = y+Ny(j);
  end
  x = x+Nx(i);
end  
