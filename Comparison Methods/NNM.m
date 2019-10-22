function f=NNM(A,F,k,filepath)

% A := weight matrix
% F := Attribute matrix
% k := number of clusters

disp('NNM is running ..... ');

N=size(F,1);
for i=1:N
    X(i,:)=F(i,:)/norm(F(i,:));
end

Y=X*X';
degs = sum(A, 2);
Dd    = sparse(1:size(A, 1), 1:size(A, 2), degs);

d=diag(Dd);
D=d*d';

M=D\A;

w=0.5;
degs(degs == 0) = eps;
Dd = spdiags(1./(degs.^0.5), 0, size(Dd, 1), size(Dd, 2));
L=sum(sum(A))*2;

Mw=Dd*(w*N/L^2*D-w*N/L*A-(1-w)/(2*N)*Y)*Dd;
Mw=sparse(Mw);
diff   = eps;
[U, ~] = eigs(Mw, k, diff);
[~,C,~,~]=kmeans_cluster('seq',U,k);
f=C;

% filename=[strcat(filepath,'\NNM_dfb_label.txt')];
% fileID = fopen(filename,'w');
% fprintf(fileID, '%d\n', C');
% fclose(fileID);
disp('NNM is done ..... ');