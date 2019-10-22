function plottimes(data, n, err,tit)

sizes = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

avgwholetimes = [];
avgparttimes = [];
avgperiters = [];
avgiters = [];

load(strcat('runs/',data,'/','totaledges',data,'_run',num2str(1),'.mat'))
totaledges = totaledges + round(n.*sizes)';


for r=1:10

load(strcat('runs/',data,'/','totaltimes',data,'_run',num2str(r),'.mat'))
avgwholetimes = [avgwholetimes; totaltimes'];


tots = zeros(1,length(sizes));
iters = zeros(1,length(sizes));
for t=1:length(sizes)
	file = strcat('runs/',data,'/',data,'_run',num2str(r),'_sizes',num2str(t),'_time_cc_duo_t.mat');
	load(file);
	tot = 0;
	iter = 0;
	for i=1:length(time)
		tot = tot + time{i}.timeink+time{i}.timeinl+time{i}.timeoutk+time{i}.timeoutl;
		iter = iter + time{i}.numiterk + time{i}.numiterl; 
	end
	tots(t) = tot;
	iters(t) = iter;
end
avgparttimes = [avgparttimes; tots];
avgperiters = [avgperiters; (tots./iters)];
avgiters = [avgiters; iters];

%[totaltimes tots']
%pause
end %r

% figure;
% boxplot(avgperiters,'labels',num2str(totaledges),'position',totaledges,'widths',10*ones(length(totaledges),1))
% pause


% q25avgwholetimes = quantile(avgwholetimes, 0.25, 1);
% q25avgparttimes = quantile(avgparttimes, 0.25, 1);
% q25avgperiters = quantile(avgperiters, 0.25, 1);
% q25avgiters = quantile(avgiters, 0.25, 1);
% 
% q75avgwholetimes = quantile(avgwholetimes, 0.75, 1);
% q75avgparttimes = quantile(avgparttimes, 0.75, 1);
% q75avgperiters = quantile(avgperiters, 0.75, 1);
% q75avgiters = quantile(avgiters, 0.75, 1);

avgwholetimes2 = sum(avgwholetimes,1) / r;
avgparttimes2 = sum(avgparttimes,1) / r;
avgperiters2 = sum(avgperiters,1) / r;
avgiters2 = sum(avgiters,1) / r;


q25avgwholetimes = std(avgwholetimes, 0, 1);
q25avgparttimes = std(avgparttimes, 0, 1);
q25avgperiters = std(avgperiters, 0, 1);
q25avgiters = std(avgiters, 0, 1);






% figure; hold all;
% plot(totaledges,avgwholetimes,'-^');
% %errorbar(totaledges,avgwholetimes,q25avgwholetimes,q75avgwholetimes)
% xlabel('Number of Edges','FontSize',14)
% ylabel('Time (sec)','FontSize',14)
% set(gca,'FontSize',14)


figure; hold all;
plot(totaledges,avgparttimes2,'-o','MarkerSize',8,'LineWidth',2);
xlabel('Number of Edges','FontSize',14)
ylabel('Time (sec)','FontSize',14)
set(gca,'FontSize',14)

figure; hold all;
plot(totaledges,avgiters2,'-*','MarkerSize',8,'LineWidth',2);
xlabel('Number of Edges','FontSize',14)
ylabel('Number of Iterations','FontSize',14)
set(gca,'FontSize',14)

figure; hold all;
plot(totaledges,avgperiters2,'square','Color','r','MarkerSize',8,'LineWidth',2);
if(err), errorbar(totaledges,avgperiters2,3*q25avgperiters,3*q25avgperiters,'LineStyle','none','Color','b','LineWidth',2); end
ss= ylim;
ylim([0 ss(2)]);
xlabel('Number of Edges |A|+|F|','FontSize',14)
ylabel('Time per Iteration (sec)','FontSize',14)
title(tit,'FontSize',14);
set(gca,'FontSize',14)


end

