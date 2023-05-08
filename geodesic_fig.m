load_pbmcdata
%%
kvalue = 30; %k-nearest neighbor parameter to use for defining topology
n_size = size(pc_data,1);

%constructing the undirected graph object that defines topology
source_edge = repmat(1:n_size,kvalue,1);
[target_edge,distances] = knnsearch(pc_data,pc_data,'k',kvalue,'NSMethod','kdtree');
target_edge = target_edge';
distances = distances';
k_graph = graph(source_edge(:),target_edge(:),distances(:),'omitselfloops');

%choosing pairs of points for which to compute geodesics
colorvec = {'r','g','b','c','m','y','k'};
%endnodes = randi(5108, [7,2]);
endnodes = [4004,1855; 3019,149; 2284,2073; 4338,4271; 4163,261; 828,3714; 3923,118]; %geodesics from paper
    for i = 1:size(endnodes,1)
        P_list_knn{i} = shortestpath(k_graph,endnodes(i,1),endnodes(i,2)); %finding the geodesic for a given pair
    end

figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    scatter(tsne_data(:,1),tsne_data(:,2),4,0.5*[1 1 1],'filled'); axis square
    hold on
    for i = 1:size(endnodes,1)
        plot(tsne_data(P_list_knn{i},1),tsne_data(P_list_knn{i},2),[colorvec{i} '-o'],'LineWidth',2)
    end
    title ('tSNE')
    
    subplot(1,2,2)
    scatter(umap_data(:,1),umap_data(:,2),4,0.5*[1 1 1],'filled'); axis square
    hold on
    for i = 1:size(endnodes,1)
        plot(umap_data(P_list_knn{i},1),umap_data(P_list_knn{i},2),[colorvec{i} '-o'],'LineWidth',2)
    end
    title('UMAP'); legend(['data'; 'Geodesic' + string(num2str([1:7]'))], 'location','southwest')