function [C,S,NScut]=SSCG(A,F,k,filepath)

% A := weight matrix
% F := Attribute matrix
% k := number of clusters

disp('SSCG is running ..... ');

W=A;
[n,d]=size(F);
s=ones(1,d)/d;
for i=1:n-1
    for j=i+1:n
        if A(i,j)
            tmp= F(i,:)-F(j,:);
            dis=(tmp*diag(s)*tmp')^0.5;
            W(i,j)=1/dis;
            W(j,i)=W(i,j);
        end
    end
end

iter=0;

for runtime=1:1
    iter=iter+1;
    degs = sum(W, 2);
    D    = sparse(1:size(W, 1), 1:size(W, 2), degs);
    
    % compute unnormalized Laplacian
    L = D - W;
    % avoid dividing by zero
    degs(degs == 0) = eps;
    % calculate inverse of D
    D = spdiags(1./degs, 0, size(D, 1), size(D, 2));
    
    % calculate normalized Laplacian
    L = D * L;
    
    diff   = eps;
    [U, lam] = eigs(L, k, diff);
    
    [~,C,~,~]=kmeans_cluster('seq',U,k);
    label(:,runtime)=C;
    
    B=zeros(n,k);
    
    for i=1:k
        index=find(C(:,1)==i);
        B(index,i)=1;
    end
    
    one=ones(n,1);
    nscut=0;
    for i=1:k
        nscut=nscut+B(:,i)'*W*(one-B(:,i))/(B(:,i)'*W*one+eps);
    end
    
    NScut(iter)=nscut;
    
    %     if iter>=2 && NScut(iter)>=NScut(iter-1)
    %         break;
    %     end
    
    index=1:d;
    
    run=1;
    
    
    for k1=1:k
        
        s=diag(ones(1,d));
        counter=1;
        while(counter<=2)
            counter=counter+1;
            %         while(1)
            for num=1:size(s,1)
                for i=1:n-1
                    for j=i+1:n
                        if A(i,j)
                            tmp= F(i,:)-F(j,:);
                            dis=(tmp*diag(s(num,:)/run)*tmp')^0.5;
                            W(i,j)=1/dis;
                            W(j,i)=W(i,j);
                        end
                    end
                end
                g(num)=B(:,k1)'*W*(one-B(:,k1))/(B(:,k1)'*W*one);
            end
            [gsorted,sub]=sort(g);
            result{run,1}=s(sub(1),:)';
            result1(run,1)=gsorted(1);
            
            tt=1;
            ind=find(s(sub(1),:)==1);
            compleSet=setdiff(index,ind);
            if isempty(compleSet)
                break;
            end
            for j=1:length(compleSet)
                s1=s(sub(1),:);
                s1(compleSet(j))=1;
                ss(tt,:)=s1;
                tt=tt+1;
            end
            
            s=[];
            g=[];
            s=ss;
            ss=[];
            run=run+1;
        end
        [~,best]=sort(result1);
        S(:,k1)=result{best(1),1};
        S(:,k1)=S(:,k1)/sum(S(:,k1));
    end
    
    %update W
    for i=1:n-1
        for j=i+1:n
            if A(i,j)
                if C(i)==C(j)
                    tmp= F(i,:)-F(j,:);
                    dis=(tmp*diag(S(:,C(i)))*tmp')^0.5;
                    W(i,j)=1/dis;
                    W(j,i)=W(i,j);
                else
                    tmp= F(i,:)-F(j,:);
                    dis=(tmp*diag(S(:,C(i)))*tmp')^0.5;
                    W(i,j)=1/dis;
                    tmp=-tmp;
                    dis=(tmp*diag(S(:,C(j)))*tmp')^0.5;
                    W(j,i)=1/dis;
                end
            end
        end
    end
end

[~,minindex]=sort(NScut);
C=label(:,minindex(1));
% filename=[strcat(filepath,'\SSCG_label.txt')];
% fileID = fopen(filename,'w');
% fprintf(fileID, '%d\n', C');
% fclose(fileID);
disp('SSCG is done ..... ');
