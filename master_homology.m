    %% install JavaPlex package (instructions at https://github.com/appliedtopology/javaplex/releases/)
    % NOTE: increase MATLAB Java Heap Memory in Settings to several GB

    addpath lib
    addpath utility
    load_javaplex
    import edu.stanford.math.plex4.*
    
    %% Compute homology groups for PBMC data

    load_pbmcdata
    num_landmark_points = 100; %Number of landmark points to use for simplicial complex
    max_dim = 3; %sets maximum homology group H_k to compute
    max_filtration_value = 4;
    
    perm_sampsize = 5; %number of landmark-point-sets to use for statistical; paper uses 100
    
    
    %store summary of each barcode plot in max_lengths_sample
    max_lengths_sample = zeros(perm_sampsize,max_dim);
    for i = 1:perm_sampsize
        landmark_selector = api.Plex4.createRandomSelector(pc_data, num_landmark_points);
        max_bar_lengths = summarize_homology(landmark_selector, max_dim, max_filtration_value, 0);
        max_lengths_sample(i,:) = max_bar_lengths';
    end
    
    %perform same computation for null scenario of gaussian distribution
    max_lengths_null = zeros(perm_sampsize,max_dim);
    for i = 1:perm_sampsize
        nulldata = normrnd(0,1,size(pc_data)).*std(pc_data);
        landmark_selector = api.Plex4.createRandomSelector(nulldata, num_landmark_points);
        max_bar_lengths = summarize_homology(landmark_selector, max_dim, max_filtration_value, 0);
        max_lengths_null(i,:) = max_bar_lengths';
    end
 %% Plots
    %make one barcode plot for one random set of landmark pts as an example
    landmark_selector = api.Plex4.createRandomSelector(pc_data, num_landmark_points);
    summarize_homology(landmark_selector, max_dim, max_filtration_value, 1);
       
    %summary plot
    figure()
        subplot(1,2,1)
        histogram(max_lengths_sample(:,2),linspace(0,30,20)); hold on
        histogram(max_lengths_null(:,2),linspace(0,30,20))
        xlabel('Normalized bar-length')
        title('H_1')

        subplot(1,2,2)
        histogram(max_lengths_sample(:,3),linspace(0,30,20)); hold on
        histogram(max_lengths_null(:,3),linspace(0,30,20))
        legend('PBMC data','Normal Distribution')
        xlabel('Normalized bar-length')
        title('H_2')


    %%
    function max_barlengths = summarize_homology(landmark_selector, max_dim, max_filtration_value,plotflag)
        import edu.stanford.math.plex4.*
        tic
        R = landmark_selector.getMaxDistanceFromPointsToLandmarks();
        stream = api.Plex4.createWitnessStream(landmark_selector, max_dim, max_filtration_value, 200);
        stream.getSize() 
        persistence = api.Plex4.getModularSimplicialAlgorithm(max_dim, 2);
        intervals = persistence.computeIntervals(stream);
        
        for k = 1:max_dim
            interval_endpts = homology.barcodes.BarcodeUtility.getEndpoints(intervals,k-1,0);
            interval_endpts(interval_endpts==Inf) = max_filtration_value;
            interval_lngths = interval_endpts(:,2) - interval_endpts(:,1);
            max_barlengths(k) = max(interval_lngths)/median(interval_lngths);
        end
        
        if plotflag==1
        options.max_filtration_value = max_filtration_value;
        options.side_by_side = false;
        plot_barcodes(intervals,options);
        end
        
        toc
        
    end
   