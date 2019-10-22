clear all;
close all;
clc;

filename='facebook';
clusters=17; %number of clusters

filepath=strcat('C:\Users\saharb63cs\Desktop\SpectralMix\Data\',filename);
addpath(filepath);

% .mat data type
% a = load('4area.mat');

% .txt data type
a.F=load('Data.txt'); % data
a.A=load('edgeWeight.txt'); % adjacency
a.GT=load('clusterlabel.txt'); % Ground truth

Nodes = size(a.F,1);
Attributes = size(a.F,2);
Adj=zeros(Nodes,Nodes);
Adj_weight=zeros(Nodes,Nodes);
edges=size(a.A,1);

for i=1:edges
        Adj(a.A(i,1)+1,a.A(i,2)+1)=1;
        Adj_weight(a.A(i,1)+1,a.A(i,2)+1)=a.A(i,3);
end

PrintPath=strcat('C:\Users\saharb63cs\workspace\SpectralMixed\real data\',filename);

[SSCG_result,~,~]=SSCG(Adj,a.F,clusters,filepath);
Print(PrintPath,'sscg',SSCG_result');

NNM_result=NNM(Adj,a.F,clusters,filepath);
Print(PrintPath,'NNM',NNM_result');

SA_result=SACluster_realdata(a.F,a.A,clusters,filepath);
Print(PrintPath,'SA',SA_result');

PICS_result=PICS(Adj,a.F);
Print(PrintPath,'PICS',PICS_result);

[BAGC_result,~,~,~]=bagc(a.A,a.F,clusters,5);
Print(PrintPath,'BAGC',BAGC_result');

% NMI_Ben(result,a.GT);

disp('done');
       
%name = strcat('HMGC_',num2str(s),'Strength_3Poisson_1Gamma_');
% xlswrite(strcat(name,'Crack.xlsx'),ComparingTable,1);
