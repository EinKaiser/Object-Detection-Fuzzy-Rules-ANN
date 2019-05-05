clc;
clear all;
close all;
warning('off','all');
% % % Get The Input Image
% [filename,pathname] = uigetfile('BAN-DATASET\*.jpg');
% img = imread([pathname,filename]);

a1 = dir('ROI_images\*.jpg');
% f1 = {};
for k = 1 : length(a1)
    a = a1(k).name;
% figure,imshow(img);title('Input Image');
img = imread(fullfile('ROI_images\',a));
% %  Scaling

%% SURF Feature Extraction 

imx = rgb2gray(img);
p1 = detectSURFFeatures(imx);
% figure,
% imshow(imx);
% hold on;
% plot(selectStrongest(p1,100));
% % Extract the features.
f2 = extractFeatures(imx,p1);
f1{k,1} = mean(f2);
end
final_feature = cell2mat(f1);
%  feature1 = final_feature;
% save('final_feature','final_feature');