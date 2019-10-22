function [sim_matrix,Laplac]  = Sim_edge( Filename )

    data_edges=load(Filename);
    numberOfEdges=size(data_edges);
    numOfNodes=max(max(max(data_edges)),max(max(data_edges)));
    sim_matrix=zeros(numOfNodes,numOfNodes);
    for e=1:numberOfEdges
        v_i=data_edges(e,1)+1;
        v_j=data_edges(e,2)+1;
        sim_matrix(v_i,v_j)=data_edges(e,3);
    end

    D=zeros(numOfNodes,numOfNodes);
    for i=1:numOfNodes
        for j=1:numOfNodes
            D(i,i)=D(i,i)+sim_matrix(i,j);
        end      
    end
    
    Laplac=eye(numOfNodes)-(power(D,-1/2)*sim_matrix*power(D,-1/2));
    
end

