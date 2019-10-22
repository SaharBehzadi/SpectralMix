k = 400;
PRUNE_SIZE = 0.001;
Restart = 0.15;
Sigma = 0.005;
Gama = 0.005;
VertexSize = 84170;
AttributeSize = 2;
AttributeValueSize1 = 3;
AttributeValueSize2 = 99;
AttributeValueSize = AttributeValueSize1 + AttributeValueSize2;
MatrixDimension = VertexSize + AttributeValueSize;
AttributeStartDimension = VertexSize + 1;
WeightOld1 = 1/2.4;
WeightOld2 = 1/2.4;
WeightSum = WeightOld1 + WeightOld2;
IncrementFactor = 1.75;
OuterIterationTime = 1;
InnerIterationTime = 10;


tic;
load Dataset.txt;
Transition = spconvert(Dataset);
clear Dataset;

load DataAttribute.txt;
Attribute = spconvert(DataAttribute);
clear DataAttribute;


AD = cell(1, InnerIterationTime);
AD{1} = zeros(VertexSize, AttributeValueSize);
for i = 2 : InnerIterationTime
	AD{i} = zeros(MatrixDimension, AttributeValueSize);
end;

ADIncrement = cell(1, InnerIterationTime);
for i = 1 : InnerIterationTime
	ADIncrement{i} = zeros(MatrixDimension, AttributeValueSize);
end;

Factor = (1.0 - Restart) * Transition;
clear Transition;
Iteration = Restart * Factor;
AD{1} = Iteration(1 : VertexSize, AttributeStartDimension : MatrixDimension);

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
    AD{i} = Iteration(:, AttributeStartDimension : MatrixDimension);
	Value = logical(Iteration > 0.0 & Iteration <= PruneSize);
	Iteration(Value) = 0.0;
	clear Value;
    Iteration = sparse(Iteration);
    PruneSize = PruneSize / IncrementFactor;
	RandomWalk = RandomWalk + Iteration;
end;
clear Iteration;


load Dataset.txt;
Transition = spconvert(Dataset);
clear Dataset;

Density = zeros(VertexSize, 1);
for i = 1 : VertexSize
    Density(i, 1) = nnz(Transition(i, 1 : VertexSize));   
end;
clear Transition;

CenterLabel = zeros(k, 1);
[Density, Index] = sortrows(Density);
for i = 1 : k
	CenterLabel(i, 1) = Index(VertexSize + 1 - i);
end;
clear Density;
clear Index;


TimeElapsed1(OuterIterationTime) = toc;
FileName = 'Runtime.txt';
Fid = fopen(FileName, 'wt');
fprintf(Fid, '%d Randomwalk %8.6f\n', OuterIterationTime, TimeElapsed1(OuterIterationTime));
fclose(Fid);
tic;



AssignmentLast = zeros(VertexSize, 1);
CenterMatrix = zeros(VertexSize, k);
while (true)
    CenterMatrix = RandomWalk(1 : VertexSize, CenterLabel);
	[MaxSimilarity, AssignmentCurrent] = max(CenterMatrix, [], 2);
    clear MaxSimilarity;
    for i = 1 : k
        AssignmentCurrent(CenterLabel(i, 1)) = i;
    end;
    
    if((nnz(AssignmentCurrent - AssignmentLast) / VertexSize) < Sigma)
        TimeElapsed2(OuterIterationTime) = toc;
        FileName = 'Runtime.txt';
        Fid = fopen(FileName, 'at');
        fprintf(Fid, '%d Clustering %8.6f\n', OuterIterationTime, TimeElapsed2(OuterIterationTime));
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
        TimeElapsed2(OuterIterationTime) = toc;
        FileName = 'Runtime.txt';
        Fid = fopen(FileName, 'at');
        fprintf(Fid, '%d Clustering %8.6f\n', OuterIterationTime, TimeElapsed2(OuterIterationTime));
        fclose(Fid);
        break;
    end;
 	WeightOld1 = WeightNew1;
	WeightOld2 = WeightNew2;
        
    
    TimeElapsed2(OuterIterationTime) = toc;
	FileName = 'Runtime.txt';
	Fid = fopen(FileName, 'at');
	fprintf(Fid, '%d Clustering %8.6f\n', OuterIterationTime, TimeElapsed2(OuterIterationTime));
	fclose(Fid);   
    OuterIterationTime = OuterIterationTime + 1;
    
    
    
    tic;
    
    PruneSize = PRUNE_SIZE;
    for i = 1 : InnerIterationTime
        ADIncrement{i} = [WeightIncrement1 * AD{i}(:, 1 : AttributeValueSize1) WeightIncrement2 * AD{i}(:, AttributeValueSize1 + 1 : AttributeValueSize)];
        AD{i} = [WeightProportion1 * AD{i}(:, 1 : AttributeValueSize1) WeightProportion2 * AD{i}(:, AttributeValueSize1 + 1 : AttributeValueSize)];
    end;

    Increment = zeros(MatrixDimension, MatrixDimension);
    Increment(1 : VertexSize, AttributeStartDimension : MatrixDimension) = ADIncrement{1};
	Value = logical(Increment > 0.0 & Increment <= PruneSize);
	Increment(Value) = 0.0;
	clear Value;
    Increment = sparse(Increment);
    PruneSize = PruneSize / IncrementFactor;
    RandomWalk = RandomWalk + Increment;


    Factor(1 : VertexSize, AttributeStartDimension : VertexSize + AttributeValueSize1) = WeightProportion1 * Factor(1 : VertexSize, AttributeStartDimension : VertexSize + AttributeValueSize1);
    Factor(1 : VertexSize, VertexSize + AttributeValueSize1 + 1 : MatrixDimension) = WeightProportion2 * Factor(1 : VertexSize, VertexSize + AttributeValueSize1 + 1 : MatrixDimension);

    for i = 2 : InnerIterationTime
        Increment = Increment * Factor;        
        AD{i} = AD{i} + Increment(1 : MatrixDimension, AttributeStartDimension : MatrixDimension);
        Increment(:, AttributeStartDimension : MatrixDimension) = Increment(:, AttributeStartDimension : MatrixDimension) + ADIncrement{i};
        Value = logical(Increment > 0.0 & Increment <= PruneSize);
        Increment(Value) = 0.0;
        clear Value;
        Increment = sparse(Increment);
        PruneSize = PruneSize / IncrementFactor;
        RandomWalk = RandomWalk + Increment;
    end;
    clear Increment;
    
    TimeElapsed1(OuterIterationTime) = toc;
	FileName = 'Runtime.txt';
	Fid = fopen(FileName, 'at');
	fprintf(Fid, '%d Increment %8.6f\n', OuterIterationTime, TimeElapsed1(OuterIterationTime));
	fclose(Fid);   
    tic;
end;
clear Factor;
clear RandomWalk;


FileName = 'Runtime.txt';
Fid = fopen(FileName, 'at');
fprintf(Fid, '\nTotal    %8.6f', sum(TimeElapsed1) + sum(TimeElapsed2));
fclose(Fid);


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

fileName = 'Entropy.txt';
fid = fopen(fileName, 'wt');
fprintf(fid, 'Entropy: %8.6f\n', EntropySum);
fclose(fid);   


load Dataset.txt;
Transition = spconvert(Dataset);
clear Dataset;
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