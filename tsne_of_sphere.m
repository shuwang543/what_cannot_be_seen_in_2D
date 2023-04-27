n_size = 10000; %number of data points to sample from the sphere

spheredata = normrnd(0,1,[n_size,3]);
spheredata = spheredata./sqrt(sum(spheredata.^2,2));

Y = tsne(spheredata);

ptsz = 3;
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
    scatter(Y(:,1),Y(:,2),ptsz,spheredata(:,1),'filled')
    daspect([1 1 1]); title('X-value')
subplot(1,3,2)
    scatter(Y(:,1),Y(:,2),ptsz,spheredata(:,2),'filled')
    daspect([1 1 1]); title('Y-value')
subplot(1,3,3)
    scatter(Y(:,1),Y(:,2),ptsz,spheredata(:,3),'filled')
    daspect([1 1 1]); title('Z-value')
