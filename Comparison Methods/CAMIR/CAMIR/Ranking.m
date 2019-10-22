function Ranking_matrix = Ranking(numOfClusters,Laplacians,Lambda,numOfNodes)

properties=size(Laplacians);
Ranking_matrix=zeros(properties,1);
U=cell(properties,1);
D=cell(properties,1);

for i=1:properties
    Laplac_p=Laplacians{i,1};
    [U{i,1},D{i,1}]=eig(Laplac_p);
end

for i=1:properties
    rest_properties=size(Laplacians);
    tr=zeros(rest_properties,1);
    for j=1:rest_properties
        temp=cell(rest_properties,1);
        for k=1:rest_properties
            if(k~=j)
                temp=temp+U{k,1}'*U{k,1};
            end
            temp=temp*Lambda;
        end
        temp_matrix=U{j,1}'*(Laplacians{j,1}+temp{j,1})*U{j,1};
        tr(j)=trace(temp_matrix);
    end
    
    [max_amount,max_index]=max(tr);
    Ranking_matrix(max_index)=rest_properties;
    Laplacians{max_index,1}=[];
end

end

