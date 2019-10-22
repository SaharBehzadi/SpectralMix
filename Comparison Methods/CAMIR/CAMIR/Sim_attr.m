function [sim_matrix,Laplac]  = Sim_attr( attr,All_edges )

numofNodes=length(attr);
sim_matrix=zeros(numofNodes,numofNodes);
sigma_values=zeros(1,numofNodes);

for i=1:numofNodes
    sigma_values=Sigma(i,attr,All_edges);
end

for i=1:numofNodes
    for j=1:numofNodes
        if(i==j)
            
            sim_matrix(i,j)=0;
        else
            sim_matrix(i,j) = exp(-(attr(i)-attr(j))*(attr(i)-attr(j))/(2*sigma_values(i)*sigma_values(j)));
        end
    end
end

D=zeros(numOfNodes,numOfNodes);
for i=1:numOfNodes
    for j=1:numOfNodes
        D(i,i)=D(i,i)+sim_matrix(i,j);
    end
end

Laplac=eye(numOfNodes)-(power(D,-1/2)*sim_matrix*power(D,-1/2));

end


function sigma_value = Sigma(node,attr, All_edges)

sigma_value=0;
connections=0;
edge_typ=size(All_edges);

for i=1:edge_typ
    edges=All_edges{i,1};
    for j=1:length(edges)
        v_i=edges(j,1);
        v_j=edges(j,2);
        if (v_i==node)
            sigma_value= sigma_value + attr(v_j);
            connections=connections+1;
        end
    end
end

sigma_value=sigma_value/connections;

end

