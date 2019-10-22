
a=load("1684.edges");

circles=17;

fid = fopen("1684_circles.txt");

for i=1:circles
    line_ex = fgetl(fid) ;
    
end


u=size(unique(a));

for j=1:u[1]
    m=min(min(a));
    s=size(a);
    for i=1:s[1]
        if(a(i,1)==m)
            a(i,1)=10000+j;
        end
        if(a(i,2)==m)
            a(i,2)=10000+j;
        end
    end
end