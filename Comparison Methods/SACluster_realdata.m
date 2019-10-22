function C=SACluster_realdata(data,edgeWeight,k,filepath)

% dataname='4area';
% k=50;

% filename=['D:\Wei\graph clustering\graph data\' dataname '\Data.txt'];
% data= importdata(filename);
% filename=['D:\Wei\graph clustering\graph data\' dataname '\edgeWeight.txt'];
% edgeWeight= importdata(filename);

% A=transfertoAffinity(edgeWeight,data);
A=transfertoAffinity(data,filepath);

[n,d]=size(data);

% dir=[strcat( filepath,'\',dataname ,'\')];
dir=[strcat( filepath,'\')];
C=SACluster_Wei(k,n,d,dir);
% [C,S,NScut]=SSCG(A,data,k);


function [DataAttribute]=transfertoAffinity(F,filepath)

    NumOfAttr=size(F,2);
    NumOfObj=size(F,1);

    count=1;
    for i=1:NumOfObj
        for j=1:NumOfAttr
            DataAttribute(count,1)=i;
            DataAttribute(count,2)=j;
            DataAttribute(count,3)=F(i,j);
            count=count+1;
        end
    end
    
fileName = strcat(filepath,'\dataAttribute.txt');
fid = fopen(fileName,'wt');
for ii = 1:size(DataAttribute,1)
    fprintf(fid,'%g\t',DataAttribute(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);