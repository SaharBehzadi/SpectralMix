function [DistGroup, DistPerNode, DistRelative] = cc_distances(k,Nx,Dnz)
% Find the distance of each cluster from every other cluster.
% Distance is the increase in encoding cost when the two groups
% are mashed into one.
% ONLY FOR SELFGRAPHS!

DistGroup = zeros(k,k);
DistPerNode = zeros(k,k);
DistRelative = zeros(k,k);

Nxy = Nx' * Nx;
Dz = Nxy - Dnz;
Pz = Dz./Nxy; Pz(~isfinite(Pz)) = 0;
Pnz = Dnz./Nxy; Pnz(~isfinite(Pnz)) = 0;
entropy_terms =  Dz .* entropy_bits(Pz) + Dnz .* entropy_bits(Pnz);

for i=1:k-1
	sum1 = sum(entropy_terms(:,i)) + sum(entropy_terms(i,:));
	for j=i+1:k
		%% Combine groups i and j
		oldsum = sum1 + sum(entropy_terms(:,j)) + sum(entropy_terms(j,:));
		oldsum = oldsum - entropy_terms(i,i) - entropy_terms(i,j) - entropy_terms(j,i) - entropy_terms(j,j);

		newnz1 = Dnz(i,:) + Dnz(j,:);
		newnz2 = Dnz(:,i) + Dnz(:,j);
		newxy  = Nx * (Nx(i)+Nx(j));
		newz1 = newxy - newnz1;
		newz2 = newxy' - newnz2;
		newpnz1 = newnz1 ./ newxy;
		newpz1  = newz1 ./ newxy;
		newpnz2 = newnz2 ./ newxy';
		newpz2  = newz2 ./ newxy';
		e1 = newnz1 .* entropy_bits(newpnz1) + newz1 .* entropy_bits(newpz1);
		e2 = newnz2 .* entropy_bits(newpnz2) + newz2 .* entropy_bits(newpz2);
		newsum = sum(e1) + sum(e2) - e1(i) - e1(j) - e2(i) - e2(j);

		midblknz = Dnz(i,i) + Dnz(i,j) + Dnz(j,i) + Dnz(j,j);
		midblkxy = (Nx(i)+Nx(j)) * (Nx(i)+Nx(j));
		newsum = newsum + midblknz*entropy_bits(midblknz/midblkxy) + (midblkxy-midblknz)*entropy_bits(1-midblknz/midblkxy);

		DistGroup(i,j) = newsum - oldsum;
		DistPerNode(i,j) = DistGroup(i,j) / (Nx(i)+Nx(j));
		DistRelative(i,j) = DistGroup(i,j) / oldsum;
	end
end

