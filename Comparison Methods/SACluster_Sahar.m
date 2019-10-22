function SACluster(A,F,k,filepath)

% A := weight matrix
% F := Attribute matrix
% k := number of clusters

% k = 400;
PRUNE_SIZE = 0.001;
Restart = 0.15;
Sigma = 0.005;
Gama = 0.005;
% VertexSize = 84170;
% AttributeSize = 2;
[VertexSize,AttributeSize]=size(F);

% AttributeValueSize1 = 3;
% AttributeValueSize2 = 99;

for i=1:AttributeSize
    AttributeValueSize_vector(i)=size(unique(F(:,i)),1);
end

% AttributeValueSize = AttributeValueSize1 + AttributeValueSize2;
AttributeValueSize=sum(AttributeValueSize_vector);
MatrixDimension = VertexSize + AttributeValueSize;
AttributeStartDimension = VertexSize + 1;
WeightOld= 1/2.4*ones(1,AttributeSize);
% WeightOld1 = 1/2.4;
% WeightOld2 = 1/2.4;
% WeightSum = WeightOld1 + WeightOld2;
WeightSum = sum(WeightOld);
IncrementFactor = 1.75;
OuterIterationTime = 1;
InnerIterationTime = 10;



tic;
Dataset = [A(:,1:2)+1,A(:,3)];
DataAttribute=matrixTransformation(F);
Attribute = spconvert(DataAttribute);
Transition = spconvert(Dataset);

Factor = (1.0 - Restart) * Transition;
clear Transition;
PruneSize = PRUNE_SIZE;
Value = logical(Factor > 0.0 & Factor <= PruneSize / Restart);
Factor(Value) = 0.0;
clear Value;
Factor = sparse(Factor);
Iteration = Restart * Factor;
RandomWalk = Iteration;
PruneSize = PruneSize / IncrementFactor;

for i = 2 : InnerIterationTime
    Iteration = Iteration * Factor;
	Value = logical(Iteration > 0.0 & Iteration <= PruneSize);
	Iteration(Value) = 0.0;
	clear Value;
    Iteration = sparse(Iteration);
    PruneSize = PruneSize / IncrementFactor;
	RandomWalk = RandomWalk + Iteration;
end;
clear Iteration;



% load Dataset.txt;
Transition = spconvert(Dataset);
% clear Dataset;

Density = zeros(VertexSize, 1);
for i = 1 : VertexSize
    Density(i, 1) = nnz(Transition(i, 1 : VertexSize));   
end;

CenterLabel = zeros(k, 1);
[Density, Index] = sortrows(Density);
for i = 1 : k
	CenterLabel(i, 1) = Index(VertexSize + 1 - i);
end;
clear Density;
clear Index;



AssignmentLast = zeros(VertexSize, 1);
CenterMatrix = zeros(VertexSize, k);
AttributeValueDistribution={};
for itr=1:AttributeSize
    AttributeValueDistribution(itr,1)=mat2cell(zeros(AttributeValueSize_vector(itr), 1),AttributeValueSize_vector(itr),1);
end
% AttributeValueDistribution1 = zeros(AttributeValueSize1, 1);
% AttributeValueDistribution2 = zeros(AttributeValueSize2, 1);

while (true)
    CenterMatrix = RandomWalk(1 : VertexSize, CenterLabel);
	[MaxSimilarity, AssignmentCurrent] = max(CenterMatrix, [], 2);
    clear MaxSimilarity;
    for i = 1 : k
        AssignmentCurrent(CenterLabel(i, 1)) = i;
    end;
       
    if((nnz(AssignmentCurrent - AssignmentLast) / VertexSize) < Sigma)
        TimeElapsed(OuterIterationTime) = toc;
        FileName = 'Runtime.txt';
        if(OuterIterationTime == 1)
            Fid = fopen(FileName, 'wt');
        else
            Fid = fopen(FileName, 'at');
        end;
        fprintf(Fid, '%d    %8.6f\n', OuterIterationTime, TimeElapsed(OuterIterationTime));
        fclose(Fid);
        break;
    end;
 	AssignmentLast = AssignmentCurrent;

    AttributeContribution=zeros(AttributeSize,1);
% 	AttributeContribution1 = 0;
% 	AttributeContribution2 = 0;   
 	for i = 1 : k
		Cluster = RandomWalk(AssignmentCurrent == i, AssignmentCurrent == i);
		RandomWalkIndex = find(AssignmentCurrent == i);
		ClusterSize = size(RandomWalkIndex, 1);
        
        Similarity = zeros(ClusterSize, 1);
        for j = 1 : ClusterSize
            Similarity = Similarity + Cluster(:, j).^2;
        end;

		[MaxSimilarity, ClusterIndex] = max(Similarity);
        clear MaxSimilarity;
		CenterLabel(i, 1) = RandomWalkIndex(ClusterIndex(1));

        for itr=1:AttributeSize
            AttributeContribution(itr)=AttributeContribution(itr) + size(find(Attribute(RandomWalkIndex, itr) == Attribute(CenterLabel(i, 1), itr)), 1);
        end
%         AttributeContribution1 = AttributeContribution1 + size(find(Attribute(RandomWalkIndex, 1) == Attribute(CenterLabel(i, 1), 1)), 1);
%         AttributeContribution2 = AttributeContribution2 + size(find(Attribute(RandomWalkIndex, 2) == Attribute(CenterLabel(i, 1), 2)), 1);
    end;

    
%     AttributeContributionSum = AttributeContribution1 + AttributeContribution2;
AttributeContributionSum = sum(AttributeContribution);
for itr=1:AttributeSize
    WeightNew(itr) = (WeightOld(itr) + WeightSum * AttributeContribution(itr) / AttributeContributionSum) / 2;
    WeightProportion(itr) = WeightNew(itr) / WeightOld(itr);
    WeightIncrement(itr) = WeightProportion(itr) - 1;
end
% 	WeightNew1 = (WeightOld1 + WeightSum * AttributeContribution1 / AttributeContributionSum) / 2;
% 	WeightNew2 = (WeightOld2 + WeightSum * AttributeContribution2 / AttributeContributionSum) / 2;
%     WeightProportion1 = WeightNew1 / WeightOld1;
%     WeightProportion2 = WeightNew2 / WeightOld2;
%     WeightIncrement1 = WeightProportion1 - 1;
%     WeightIncrement2 = WeightProportion2 - 1;
    
Thr = 0;
for itr=1:AttributeSize
    Thr=Thr+abs(WeightNew(itr) - WeightOld(itr));
end
Thr=Thr/AttributeSize;

    if (Thr< Gama)
        TimeElapsed(OuterIterationTime) = toc;
        FileName = 'Runtime.txt';
        if(OuterIterationTime == 1)
            Fid = fopen(FileName, 'wt');
        else
            Fid = fopen(FileName, 'at');
        end;
        fprintf(Fid, '%d    %8.6f\n', OuterIterationTime, TimeElapsed(OuterIterationTime));
        fclose(Fid);
        break;
    end;
    for itr=1:AttributeSize
        WeightOld(itr) = WeightNew(itr);
    end
%  	WeightOld1 = WeightNew1;
% 	WeightOld2 = WeightNew2;
    
    
    
    TimeElapsed(OuterIterationTime) = toc;
	FileName = 'Runtime.txt';
	if(OuterIterationTime == 1)
        Fid = fopen(FileName, 'wt');
    else
        Fid = fopen(FileName, 'at');
	end;
	fprintf(Fid, '%d    %8.6f\n', OuterIterationTime, TimeElapsed(OuterIterationTime));
	fclose(Fid);  
    OuterIterationTime = OuterIterationTime + 1;
    
    
    
    tic;
    PruneSize = PRUNE_SIZE;
    count=AttributeStartDimension;
    for itr=1:AttributeSize
        if itr==AttributeSize
            Factor(1 : VertexSize, count : MatrixDimension) = WeightProportion(itr) * Factor(1 : VertexSize, count : MatrixDimension);           
        else
            Factor(1 : VertexSize, count : count-1 + AttributeValueSize_vector(itr)) = WeightProportion(itr) * Factor(1 : VertexSize, count : count-1 + AttributeValueSize_vector(itr));           
        end
           count=count +AttributeValueSize_vector(itr);
    end

%     Factor(1 : VertexSize, AttributeStartDimension : VertexSize + AttributeValueSize1) = WeightProportion1 * Factor(1 : VertexSize, AttributeStartDimension : VertexSize + AttributeValueSize1);
%     Factor(1 : VertexSize, VertexSize + AttributeValueSize1 + 1 : MatrixDimension) = WeightProportion2 * Factor(1 : VertexSize, VertexSize + AttributeValueSize1 + 1 : MatrixDimension);
    Iteration = Restart * Factor;
    RandomWalk = Iteration;
    PruneSize = PruneSize / IncrementFactor;


    for i = 2 : InnerIterationTime
        Iteration = Iteration * Factor;
    	Value = logical(Iteration > 0.0 & Iteration <= PruneSize);
    	Iteration(Value) = 0.0;
    	clear Value;
        Iteration = sparse(Iteration);
        PruneSize = PruneSize / IncrementFactor;
        RandomWalk = RandomWalk + Iteration;
    end;
    clear Iteration;


%     load Dataset.txt;
    Transition = spconvert(Dataset);
%     clear Dataset;
end;


	FileName = 'Runtime.txt';
	Fid = fopen(FileName, 'at');
	fprintf(Fid, '\nTotalTime    %8.6f', sum(TimeElapsed));
	fclose(Fid);
clear Factor;
clear RandomWalk;


EntropySum = 0;

for itr=1:AttributeSize
    AttributeValueDistribution(itr,1)=mat2cell(zeros(AttributeValueSize_vector(itr), 1),AttributeValueSize_vector(itr),1);
end

% AttributeValueDistribution1 = zeros(AttributeValueSize1, 1);
% AttributeValueDistribution2 = zeros(AttributeValueSize2, 1);

for i = 1 : k
	TransitionIndex = find(AssignmentCurrent == i);
	ClusterSize = size(TransitionIndex, 1);
    
    for itr=1:AttributeSize
        Entropy = 0;
        for j = 1 : AttributeValueSize_vector(itr)
            AttributeValueDistribution(itr,j, 1) = size(find(Attribute(TransitionIndex, 1) == (j-1)), 1) / ClusterSize;
            if (AttributeValueDistribution(itr,j, 1) > 0)
                Entropy = Entropy - WeightOld(itr) * AttributeValueDistribution(itr,j, 1) * log(AttributeValueDistribution(itr,j, 1)) / log(2);
            end
        end
        EntropySum = EntropySum + Entropy * ClusterSize / VertexSize / WeightSum;
    end

% 	Entropy = 0;
% 	for j = 1 : AttributeValueSize2
%         AttributeValueDistribution2(j, 1) = size(find(Attribute(TransitionIndex, 2) == j), 1) / ClusterSize;
%         if (AttributeValueDistribution2(j, 1) > 0)            
%             Entropy = Entropy - WeightOld2 * AttributeValueDistribution2(j, 1) * log(AttributeValueDistribution2(j, 1)) / log(2);
%         end; 
% 	end;
% 	EntropySum = EntropySum + Entropy * ClusterSize / VertexSize / WeightSum;
end

fileName = 'Entropy.txt';
fid = fopen(fileName, 'wt');
fprintf(fid, 'Entropy: %8.6f\n', EntropySum);
fclose(fid);  


% load Dataset.txt;
Transition = spconvert(Dataset);
% clear Dataset;
Cohensive = 0;
for i = 1 : k
	Cohensive = Cohensive + nnz(Transition(AssignmentCurrent == i, AssignmentCurrent == i));
end;
Cohensive = Cohensive / nnz(Transition(1 : VertexSize, 1 : VertexSize));
clear Transition;

fileName = 'Cohensive.txt';
fid = fopen(fileName, 'wt');
fprintf(fid, 'Density: %8.6f\n', Cohensive);
fclose(fid);

function [DataAttribute]=matrixTransformation(F)

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
    
%     count=0;
%     for i=1:NumOfObj
%         for j=1:NumOfObj
%             if(A(i,j)~=0)
%                 Dataset(count,1)=i;
%                 Dataset(count,2)=j;
%                 Dataset(count,3)=A(i,j);
%             end
%         end
%     end


