load_pbmcdata
%% compute the topology of the data
    kvalue = 30; %k-nearest neighbor parameter to use for defining topology
    n_size = size(pc_data,1);
    
    %constructing the undirected graph object that defines topology
    source_edge = repmat(1:n_size,kvalue,1);
    [target_edge,distances] = knnsearch(pc_data,pc_data,'k',kvalue,'NSMethod','kdtree');
    target_edge = target_edge';
    distances = distances';
    k_graph = graph(source_edge(:),target_edge(:),distances(:),'omitselfloops');
    
%% compute pairwise distances of cluster centroids for different defn's of distance
        
        % 5-clusters
        sequence5_full = pdist(grpstats(pc_data,clusters5)); %10-PC Euclidean Distance
        sequence5_geo = pdist_geodesic(k_graph, knnsearch(pc_data,grpstats(pc_data,clusters5))); %kNN-graph geodesic distance
        sequence5_tsne = pdist(grpstats(tsne_data,clusters5)); %tSNE distance
        sequence5_umap = pdist(grpstats(umap_data,clusters5)); %UMAP distance
        
        % 6 clusters
        sequence6_full = pdist(grpstats(pc_data,clusters6));
        sequence6_geo = pdist_geodesic(k_graph, knnsearch(pc_data,grpstats(pc_data,clusters6)));
        sequence6_tsne = pdist(grpstats(tsne_data,clusters6));
        sequence6_umap = pdist(grpstats(umap_data,clusters6));
        
        % 7 clusters
        sequence7_full = pdist(grpstats(pc_data,clusters7));
        sequence7_geo = pdist_geodesic(k_graph, knnsearch(pc_data,grpstats(pc_data,clusters7)));
        sequence7_tsne = pdist(grpstats(tsne_data,clusters7));
        sequence7_umap = pdist(grpstats(umap_data,clusters7));
    
    % figures use tSNE projection, but replacing tsne_data with umap_data
    % will show UMAP projection
    figure('units','normalized','outerposition',[0 0 1 1]) 
    subplot(2,3,1)
        scatter(tsne_data(:,1),tsne_data(:,2),4,clusters5,'filled'); 
        daspect([1 1 1]); colormap(jet)
    subplot(2,3,2)
        scatter(tsne_data(:,1),tsne_data(:,2),4,clusters6,'filled'); 
        daspect([1 1 1]); colormap(jet)
    subplot(2,3,3)
        scatter(tsne_data(:,1),tsne_data(:,2),4,clusters7,'filled'); 
        daspect([1 1 1]); colormap(jet)

    %example here compares correlation of 10-PC Euclidean Distance with tSNE distance
    %replace the sequences with whichever ones are of interest
    subplot(2,3,4)
        scatter_pdist(sequence5_full,sequence5_tsne) 
    subplot(2,3,5)
        scatter_pdist(sequence6_full,sequence6_tsne) 
    subplot(2,3,6)
        scatter_pdist(sequence7_full,sequence7_tsne) 

%% Auxilary functions
function scatter_pdist(seq_hiD,seq_2D)  
    [rho, p] = corr(seq_hiD', seq_2D','Type','Spearman');
    scatter(seq_hiD', seq_2D', 20, 'filled');
    title(['Spearman Corr=' num2str(rho) '; p-val=' num2str(p)])
    xlabel('Hi-Dim Distance'); ylabel('Visual Distance')
    axis square
end
    

function pdist_geo = pdist_geodesic(graph, nodes)
    endnodes = nchoosek(nodes,2);
    for i = 1:size(endnodes,1)
        [~,pdist_geo(i)] = shortestpath(graph,endnodes(i,1),endnodes(i,2));
    end
end
     