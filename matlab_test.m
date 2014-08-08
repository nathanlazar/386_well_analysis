%% Clear
clear all;

%% Read and adjust images
pref = 'C:\Users\lazar\Desktop\Thesis\CSMA-MEMA\386images\HCC1143.KRTVIM.NVS1.dmso\C5--W00053--P00001--Z00000--T00000--'
dapi = imadjust(imread(strcat(pref, 'DAPI.tif')));
alexa647 = imadjust(imread(strcat(pref, 'Alexa 647.tif')));
alexa488 = imadjust(imread(strcat(pref, 'Alexa 488.tif')));
alexa568 = imadjust(imread(strcat(pref, 'Alexa 568.tif')));

%% Display adjusted images
figure;
name = {'DAPI' 'Alexa 647' 'Alexa 488' 'Alexa 568'};
image = {dapi alexa647 alexa488 alexa568};
for i=1:4,
    h(i) = subplot(2,2,i); 
    imshow(image{i}); 
    title(name{i});
end
linkaxes(h, 'xy') %allows to zoom images together

%% Naive segmentation
figure;
for i=1:4,
    h(i) = subplot(2,2,i);
    hist(double(reshape(image{i}, 1, 1024*1344)),200);
    title(name{i});
end
linkaxes(h, 'xy')

for i=1:4,
    figure;
    h(1) = subplot(2,2,1); 
    imshow(image{i}); title(name{i});
    h(2) = subplot(2,2,2);
    seg2 = image{i} >= 1500;
    imshow(seg2); title('pixels above 1500');
    h(3) = subplot(2,2,3);
    seg2 = image{i} >= 15000;
    imshow(seg2); title('pixels above 15000');
    h(4) = subplot(2,2,4); 
    seg1 = dapi==65535;
    imshow(seg1); title('pixels at max');
    linkaxes(h, 'xy');
end
%% Naive segmentation of cells (Alexa 568)
figure;
hist(double(reshape(alexa568, 1, 1024*1344)),200);
figure;
h(1) = subplot(2,2,1); 
imshow(alexa568); title('Alexa 568')
h(2) = subplot(2,2,2);
seg2 = dapi>=1500;
imshow(seg2); title('pixels above 1500')
h(3) = subplot(2,2,3);
seg2 = dapi>=15000;
imshow(seg2); title('pixels above 15000')
h(4) = subplot(2,2,4); 
seg1 = dapi==65535;
imshow(seg1); title('pixels at max')

linkaxes(h, 'xy')
%%
quantile(single(reshape(alexa568, 1, 1024*1344)), 10)
%%
q = quantile(single(reshape(dapi, 1, 1024*1344)), 200)
%% Denoise (using median filter)
for j=1:4,    
    figure;
    h(1) = subplot(2,3,1);
    imshow(image{j});
    title(name{j})
    for i=1:5,
        h(i+1) = subplot(2,3,i+1);
        imshow(medfilt2(image{j}, [i,i])); 
        title(strcat(num2str(i), 'X', num2str(i), ' filter'));
    end
    linkaxes(h, 'xy')
end

%% Edge detection and fill
figure; title('canny')
for i=1:4,
    h(i) = subplot(3,4,i); 
    imshow(image{i}); title(name{i})
    h(i+4) = subplot(3,4,i+4); 
    imedge{i} = edge(image{i}, 'canny');
    imshow(imedge{i}); title(strcat(name{i}, ' canny'));
    h(i+8) = subplot(3,4,i+8);
    imshow(imfill(imedge{i}, 'holes')); title(strcat(name{i}, ' filled'));
end
linkaxes(h, 'xy')

%% Dapi edge detection, filling
figure;
h(1) = subplot(4,4,1); 
imshow(dapi); title('DAPI');
e_types = {'canny', 'sobel', 'prewitt', 'roberts', 'log', 'zerocross'};
placement = {1,2,3,4,10,11,12};
for i=1:6,
    h(i+1) = subplot(4,4,placement{i+1});
    e{i} = edge(dapi, e_types{i});
    imshow(e{i}); title(e_types{i});
    h(i+7) = subplot(4,4,placement{i+1}+4);
    imshow(imfill(e{i}, 'holes')); title(strcat(e_types{i}, ' filled'));
end 
linkaxes(h, 'xy')

%% Filter out objects by size
size_min = 200;
size_max = 2000;
[L, n] = bwlabel(imfill(e{6}, 'holes'), 4);
j = 0;
for i=1:n,
    s = sum(sum(L==i)); %number of pixels in object
    if s < size_min | s > size_max
        L(find(L==i)) = 0; %set mask values to zero
    else
        j = j+ 1;
    end
end
figure;
imshow(L)

%% Crop out individual objects
objs = unique(L(L~=0))
i = 1;
for obj=objs',
    col = max(L==obj);
    x = find(col);
    row = max(L==obj, [], 2);
    y = find(row);
    a = dapi;
    a(find(L~=obj)) = 0; %sets non-object pixels to black
    b = imcrop(a, [min(x)-5 min(y)-5 max(x)-min(x)+10 max(y)-min(y)+10]);
    imwrite(b, strcat('DAPI_objects/DAPI_object', num2str(i), '.tiff'));
    i = i + 1;
end

%% Play with adding images
figure;
dapi_canny = imadd(dapi, uint16(e{6})*2^16-1); %e{6} are zero-cross edges
h(1) = subplot(1,3,1); imshow(dapi); title('dapi')
h(2) = subplot(1,3,2); imshow(e{6}); title('edges')
h(3) = subplot(1,3,3); imshow(dapi_canny); title('sum')
linkaxes(h, 'xy')

%% Convert to color 
dapi_blue = repmat(dapi,[1 1 3]);
%alternatively: dapi_blue = cat(3,dapi2,dapi2,dapi2);
%zero out the red and green channels
dapi_blue(:,:,1) = 0;
dapi_blue(:,:,1) = 0;
figure;
imshow(dapi_blue)

%% Combine blue dapi with green and red channels
figure;
comb1 = cat(3, alexa647, alexa488, dapi);
h(1) = subplot(1,2,1); 
imshow(comb1); title('647 red, 488 green, DAPI blue')
comb2 = cat(3, alexa647, alexa568, dapi);
h(2) = subplot(1,2,2); 
imshow(comb2); title('647 red, 568 green, DAPI blue')
linkaxes(h, 'xy')

%% Phase Congruency (edge detection method by Peter Kovesi)
figure;
for i=1:4,
    h(i) = subplot(2,2,i); 
    imshow(phasecong(image{i})); 
    title(name{i});
end
linkaxes(h, 'xy') %allows to zoom images together

%% Use DBScan to find objects
% we first need to translate pixels above a certain threshold into points
% w/ x,y as their coordinates
filt = dapi > 2500;
filt = filt(200:400,200:400)
figure;
imshow(filt);
[x,y] = ind2sub(size(filt), find(filt));
%%
[C, ptsC, centres] = dbscan([x,y]', 3, 10);
%%
figure;
scatter(centres(2,:), 200-centres(1,:),'blue','filled');
hold on;
%imshow(filt(C{1}))
%scatter(centres(1,:), centres(2,:), 'red'); %just testing, changes these to the points from ptsC
hold off;
