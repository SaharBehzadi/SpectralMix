function f=SACluster_Wei(k,VertexSize,AttributeSize,dir)

VertexSize=VertexSize-AttributeSize;
AttributeValueSize1 = floor(AttributeSize/2);
AttributeValueSize2 = AttributeSize-AttributeValueSize1;
AttributeValueSize = AttributeValueSize1 + AttributeValueSize2;
MatrixDimension = VertexSize + AttributeValueSize;
AttributeStartDimension = VertexSize + 1;


PRUNE_SIZE = 0.001;
Restart = 0.15;
Sigma = 0.005;
Gama = 0.005;
WeightOld1 = 1/2.4;
WeightOld2 = 1/2.4;
WeightSum = WeightOld1 + WeightOld2;
IncrementFactor = 1.75;
OuterIterationTime = 1;
InnerIterationTime = 10;


tic;
filename=strcat(dir,'Dataset.txt');
Dataset=load(filename);
Dataset = [Dataset(:,1:2)+1,Dataset(:,3)];
Transition = spconvert(Dataset);
clear Dataset;
clear filename;

filename=strcat(dir,'dataAttribute.txt');
Attribute = load(filename);
clear filename;


Factor = (1.0 - Restart) * Transition;
clear Transition;
PruneSize = PRUNE_SIZE;
Value(1,:) = logical(Factor(1,:) > 0.0 & Factor(1,:) <= PruneSize / Restart);
Value=sparse(Value);
for i=2:size(Factor,1)
    Value(i,:) = logical(Factor(i,:) > 0.0 & Factor(i,:) <= PruneSize / Restart);
end
Factor(Value) = 0.0;
clear Value;
Factor = sparse(Factor);
Iteration = Restart * Factor;
RandomWalk = Iteration;
PruneSize = PruneSize / IncrementFactor;

for i = 2 : InnerIterationTime
    Iteration = Iteration * Factor;
    for j=2:size(Iteration,1)
    Value(j,:) = logical(Iteration(j,:) > 0.0 & Iteration(j,:) <= PruneSize);
    end
    Iteration(Value) = 0.0;
    clear Value;
    Iteration = sparse(Iteration);
    PruneSize = PruneSize / IncrementFactor;
    RandomWalk = RandomWalk + Iteration;
end;
clear Iteration;



filename=strcat(dir,'Dataset.txt');
Dataset=load(filename);
Dataset = [Dataset(:,1:2)+1,Dataset(:,3)];
Transition = spconvert(Dataset);
clear Dataset;
clear filename;

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

filename1=strcat(dir,'SAcluster');
mkdir(filename1);
FileName=strcat(filename1,'\Runtime.txt');


AssignmentLast = zeros(VertexSize, 1);
CenterMatrix = zeros(VertexSize, k);
AttributeValueDistribution1 = zeros(AttributeValueSize1, 1);
AttributeValueDistribution2 = zeros(AttributeValueSize2, 1);
while (true)
    CenterMatrix = RandomWalk(1 : VertexSize, CenterLabel);
    [MaxSimilarity, AssignmentCurrent] = max(CenterMatrix, [], 2);
    clear MaxSimilarity;
    for i = 1 : k
        AssignmentCurrent(CenterLabel(i, 1)) = i;
    end;
    
    if((nnz(AssignmentCurrent - AssignmentLast) / VertexSize) < Sigma)
        TimeElapsed(OuterIterationTime) = toc;
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
    
    
    AttributeContribution1 = 0;
    AttributeContribution2 = 0;
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
        
        AttributeContribution1 = AttributeContribution1 + size(find(Attribute(RandomWalkIndex, 1) == Attribute(CenterLabel(i, 1), 1)), 1);
        AttributeContribution2 = AttributeContribution2 + size(find(Attribute(RandomWalkIndex, 2) == Attribute(CenterLabel(i, 1), 2)), 1);
    end;
    
    
    AttributeContributionSum = AttributeContribution1 + AttributeContribution2;
    WeightNew1 = (WeightOld1 + WeightSum * AttributeContribution1 / AttributeContributionSum) / 2;
    WeightNew2 = (WeightOld2 + WeightSum * AttributeContribution2 / AttributeContributionSum) / 2;
    WeightProportion1 = WeightNew1 / WeightOld1;
    WeightProportion2 = WeightNew2 / WeightOld2;
    WeightIncrement1 = WeightProportion1 - 1;
    WeightIncrement2 = WeightProportion2 - 1;
    
    if (((abs(WeightNew1 - WeightOld1) + abs(WeightNew2 - WeightOld2)) / AttributeSize) < Gama)
        TimeElapsed(OuterIterationTime) = toc;
        if(OuterIterationTime == 1)
            Fid = fopen(FileName, 'wt');
        else
            Fid = fopen(FileName, 'at');
        end;
        fprintf(Fid, '%d    %8.6f\n', OuterIterationTime, TimeElapsed(OuterIterationTime));
        fclose(Fid);
        break;
    end;
    WeightOld1 = WeightNew1;
    WeightOld2 = WeightNew2;
    
    
    
    TimeElapsed(OuterIterationTime) = toc;
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
    Factor(1 : VertexSize, AttributeStartDimension : VertexSize + AttributeValueSize1) = WeightProportion1 * Factor(1 : VertexSize, AttributeStartDimension : VertexSize + AttributeValueSize1);
    Factor(1 : VertexSize, VertexSize + AttributeValueSize1 + 1 : MatrixDimension) = WeightProportion2 * Factor(1 : VertexSize, VertexSize + AttributeValueSize1 + 1 : MatrixDimension);
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
    
    
    filename=strcat(dir,'Dataset.txt');
    Dataset=load(filename);
    Dataset = [Dataset(:,1:2)+1,Dataset(:,3)];
    Transition = spconvert(Dataset);
    clear Dataset;
    clear filename;
end;


Fid = fopen(FileName, 'at');
fprintf(Fid, '\nTotalTime    %8.6f', sum(TimeElapsed));
fclose(Fid);
clear Factor;
clear RandomWalk;


EntropySum = 0;
AttributeValueDistribution1 = zeros(AttributeValueSize1, 1);
AttributeValueDistribution2 = zeros(AttributeValueSize2, 1);
for i = 1 : k
    TransitionIndex = find(AssignmentCurrent == i);
    ClusterSize = size(TransitionIndex, 1);
    Entropy = 0;
    for j = 1 : AttributeValueSize1
        AttributeValueDistribution1(j, 1) = size(find(Attribute(TransitionIndex, 1) == (j-1)), 1) / ClusterSize;
        if (AttributeValueDistribution1(j, 1) > 0)
            Entropy = Entropy - WeightOld1 * AttributeValueDistribution1(j, 1) * log(AttributeValueDistribution1(j, 1)) / log(2);
        end;
    end;
    EntropySum = EntropySum + Entropy * ClusterSize / VertexSize / WeightSum;
    Entropy = 0;
    for j = 1 : AttributeValueSize2
        AttributeValueDistribution2(j, 1) = size(find(Attribute(TransitionIndex, 2) == j), 1) / ClusterSize;
        if (AttributeValueDistribution2(j, 1) > 0)
            Entropy = Entropy - WeightOld2 * AttributeValueDistribution2(j, 1) * log(AttributeValueDistribution2(j, 1)) / log(2);
        end;
    end;
    EntropySum = EntropySum + Entropy * ClusterSize / VertexSize / WeightSum;
end;

fileName = strcat(filename1,'\Entropy.txt');
fid = fopen(fileName, 'wt');
fprintf(fid, 'Entropy: %8.6f\n', EntropySum);
fclose(fid);
clear fileName;


filename=strcat(dir,'Dataset.txt');
Dataset=load(filename);
Dataset = [Dataset(:,1:2)+1,Dataset(:,3)];
Transition = spconvert(Dataset);
clear Dataset;
clear filename;

Cohensive = 0;
for i = 1 : k
    Cohensive = Cohensive + nnz(Transition(AssignmentCurrent == i, AssignmentCurrent == i));
end;
Cohensive = Cohensive / nnz(Transition(1 : VertexSize, 1 : VertexSize));
clear Transition;

fileName = strcat(filename1,'\Density.txt');
fid = fopen(fileName, 'wt');
fprintf(fid, 'Density: %8.6f\n', Cohensive);
fclose(fid);
clear fileName;

% for ii=1:AttributeSize
%     fd=randperm(k);
%     onelast(ii,1)=fd(1);
% end
AssignmentCurrent=[AssignmentCurrent;ones(AttributeSize,1)];

fileName = strcat(filename1,'\asignment.txt');
fileID = fopen(fileName,'w');
fprintf(fileID, '%d\n', AssignmentCurrent');
fclose(fileID);


f=AssignmentCurrent;