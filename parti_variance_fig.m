load_pbmcdata
%run PARTI code (Y. Hart, et al., Inferring biological tasks using Pareto analysis of high-dimensional data. Nat. Methods 2015)
   [~, arcOrig] = ParTI_lite(pc_data, 1, 10);
   [arc_near,~] = knnsearch(arcOrig,pc_data);
   [arc_rep,~] = knnsearch(pc_data,arcOrig);
   close all

%% Projection of PBMC data that helps visualize archetypes

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,[1 2])
xvec = [1;1;0;1;-1;0;0;1;0;0]; % Weight of PCs to visualize archetypes
yvec = [0;1;2;-1;1;0;0;0;0;0]; % Weight of PCs to visualize archetypes
scatter(pc_data*xvec,pc_data*yvec,4,clusters6,'filled'); 
hold on; scatter(arcOrig*xvec,arcOrig*yvec,100,'m','pentagram','filled')
scatter(pc_data(arc_rep,:)*xvec,pc_data(arc_rep,:)*yvec,50,'m','filled')
text(pc_data(arc_rep,:)*xvec,pc_data(arc_rep,:)*yvec,num2str([1:6]'),'FontSize',16,'HorizontalAlignment','right','VerticalAlignment','bottom')
colormap(jet); axis square
xlabel('PC 1/2/4/5/8'); ylabel('PC 2/3/4/5')
legend('PBMC data','Projected archetypes','Nearest data pt to archetype')

%% Computing the variance of each Cell-type cluster along the edges between archetypes

k = 6;
allpairs = nchoosek(1:k,2);
n_k = size(allpairs,1);
percentvar = zeros(k,n_k);

unitvecs = zeros(n_k,10);
arcrep = pc_data(arc_rep,:);
for i = 1:n_k
    unitvecs(i,:) = arcrep(allpairs(i,1),:) - arcrep(allpairs(i,2),:);
    unitvecs(i,:) = unitvecs(i,:)/sqrt(sum(unitvecs(i,:).^2,2));
end

for j = 1:k
    clusterdata = pc_data(clusters6==j,:);
    var_tot = sum(trace(cov(clusterdata)));
    percentvar(j,:) = std(clusterdata*unitvecs').^2/var_tot;
end

subplot(1,3,3)
imagesc(percentvar')
colormap(jet); colorbar()
set(gca,'ytick',1:size(allpairs,1))
set(gca,'yticklabels',num2str(allpairs))
set(gca,'xtick',1:k)
set(gca,'xticklabels',{'T-cells','Monocytes','NK','B','Eryth','Platelet'})
ylabel('1-simplex edges between archetypes')
xlabel('K-means clusters')
title('Percent total variance along 1-simplex edges')


