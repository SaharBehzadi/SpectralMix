function reto=NMI_Ben(one,two)

if(size(one)~= size(two))
    Disp('Sizes do not match');
end

maxone = max(one);
maxtwo = max(two);

count = zeros(maxone+1,maxtwo+1);
for i=1:size(one)
    count(one(i),two(i))= count(one(i),two(i))+1;
end

bj = zeros(maxtwo+1,1);
ai = zeros(maxone+1,1);

for m=1:maxtwo+1
    for l=1:maxone+1
        bj(m)=bj(m)+count(l,m);
    end
end

for m=1:(maxone+1)
    for l=1:(maxtwo+1)
        ai(m)=ai(m)+count(m,l);
    end
end

N=0;

for i=1:length(ai)
    N=N+ai(i);
end

HU = 0;
for l=1:length(ai)
    c=(ai(l)/N);
    if(c>0)
        HU=HU-c*log(c);
    end
end

HV = 0;
for l=1:length(bj)
    c=(bj(l)/N);
    if(c>0)
        HV=HV-c*log(c);
    end
end

HUstrichV=0;
for i=1:maxone+1
    for j=1:(maxtwo+1)
        if(count(i,j)>0)
            HUstrichV=HUstrichV-count(i,j)/N*log(((count(i,j))/(bj(j))));
        end
    end
end

IUV = HU-HUstrichV;
reto = IUV/(max(HU, HV));

disp(strcat('NMI: ', num2str(reto)));