%files found in https://www.10xgenomics.com/resources/datasets/pbmcs-3p_acda_sepmate-3-1-standard
%specifically the unzipped file 3p_ACDA_SepMate_analysis.tar.gz under
%"Clustering Analysis"

temp = readtable('analysis/pca/10_components/projection.csv');
pc_data = temp{:,2:11};

temp = readtable('analysis/tsne/2_components/projection.csv');
tsne_data = temp{:,2:3};

temp = readtable('analysis/umap/2_components/projection.csv');
umap_data = temp{:,2:3};

temp = readtable('analysis/clustering/kmeans_5_clusters/clusters.csv');
clusters5 = temp{:,2};
temp = readtable('analysis/clustering/kmeans_6_clusters/clusters.csv');
clusters6 = temp{:,2};
temp = readtable('analysis/clustering/kmeans_7_clusters/clusters.csv');
clusters7 = temp{:,2};

clear temp