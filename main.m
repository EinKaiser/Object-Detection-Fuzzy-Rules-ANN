clc;
clear;
close all;
%% get video
delete('Frame\*.jpg');
delete('Output Frame\*.jpg');
label=ones(43,1);
label(14:24)=2;
label(25:32)=2;
load feat.mat
model=fitcknn(f,label);
[filename,pathname] = uigetfile({'test video\*.avi'},'Select A Video File');
I = VideoReader([pathname,filename]);
nFrames = I.numberofFrames;
vidHeight =  I.Height;
vidWidth =  I.Width;
mov(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
    'colormap', []);
for k = 1:nFrames
    mov(k).cdata = read( I, k);
    mov(k).cdata = imresize(mov(k).cdata,[256,320]);
    imwrite(mov(k).cdata,['Frame\',sprintf('%06d.jpg',k)]);
end
%% get image
WWW=dir('Frame\*.jpg');
count=1;
SS=[];
ccf=0;
% figure;
outp=0;
for j=1:length(WWW)-2
    j
    Y1= fullfile('Frame',WWW(j).name);
    a1=imresize(imread(Y1),[256 320]);
    Y2 = fullfile('Frame',WWW(j+1).name);
    if length(WWW)>1000
        outp=1;
    end
    b1=imresize(imread(Y2),[256 320]);
    Y3 = fullfile('Frame',WWW(j+2).name);
    c3=imresize(imread(Y3),[256 320]);
    a=(rgb2gray(a1));
    b=(rgb2gray(b1));
    %% Preprocessing
    
    [r1,c1]=size(a);
    blockSizeR=r1/4;
    blockSizeC=c1/5;
    
    coordinating=forward(a1,b1);
    coord=forward(b1,c3);
    %% Moving object segmentation
    if ~isempty(coordinating) && ~isempty(coord)
        % Object region tracking
        coordinatefor=coordinating(:,1:2);
        coordinatebac=coord(:,1:2);
        forwa=zeros(256,320);
        backwa=zeros(256,320);
        for i=1:size(coordinatefor,1)
            c=coordinatefor(i,1);
            r=coordinatefor(i,2);
            forwa(r:r+7,c:c+7)=1;
        end
        for i=1:size(coordinatebac,1)
            c=coordinatebac(i,1);
            r=coordinatebac(i,2);
            backwa(r:r+8,c:c+8)=1;
        end
        coordination=[];
        for rowi = 1 : blockSizeR : r1
            for coli = 1 : blockSizeC : c1
                row1 = rowi;
                row2 = row1 + blockSizeR - 1;
                col1 = coli;
                col2 = col1 + blockSizeC - 1;
                % Extract out the block into a single subimage.
                mainimage = forwa(row1:row2, col1:col2);
                referenceimage = backwa(row1:row2, col1:col2);
                intersection=mainimage.*referenceimage;
                intersection=sum(intersection(:));
                mainimage=sum(mainimage(:));
                if (intersection/mainimage)>0.2
                    coordination=[coordination;coli rowi 64 64];
                end
            end
        end
        if ~isempty(coordination)
            coordinate2=[];
            for i=1:size(coordination,1)
                R=coordination(i,2);
                C=coordination(i,1);
                for ii=1:size(coordinatefor,1);
                    r=coordinatefor(ii,2);
                    c=coordinatefor(ii,1);
                    if r>R && r<R+65 && c>C && c<C+65
                        coordinate2=[coordinate2; c r 7 7];
                    end
                end
            end
            % object Boundary Refinement
            backwa=zeros(256,320);
            for i=1:size(coordinate2,1)
                c=coordinate2(i,1);
                r=coordinate2(i,2);
                backwa(r:r+7,c:c+7)=1;
            end
            backwa = bwareaopen(backwa, 350);
            [r,c]=find(backwa==1);
            cod=zeros(256,320);
            for i=1:size(r,1)
                R=r(i);
                C=c(i);
                cod(R,C)=1;
            end
            cod1=zeros(256,320);
            for ii=1:size(coordinate2,1);
                r=coordinate2(ii,2);
                c=coordinate2(ii,1);
                cod1(r,c)=1;
            end
            coordinate=cod.*cod1;
            [r,c]=find(coordinate==1);
            coordinate=[];
            refine=zeros(256,320);
            for i=1:length(r)
                R=r(i);
                C=c(i);
                refine(R:R+7,C:C+7)=1;
                coordinate=[coordinate;C R 7 7];
            end
            L = bwlabel(refine,8);
            numb=max(max(L));
            S1=zeros(256,320);
            for i=1:numb
                [r2,c2] = find(L == i);
                S=zeros(256,320);
                for ii=1:length(r2)
                    R=r2(ii);
                    C=c2(ii);
                    S(R,C)=1;
                end
                co=cod.*cod1;
                cord=S.*co;
                [r,c]=find(cord==1);
                [Avgdepth,CUdepth]=CUDEPTH(a,b,r,c);
                
                j0=min(min(r));
                j1=max(max(r));
                CUdepth2=0;CUdepth3=0;CUdepth4=0;CUdepth5=0;p1=0;
                while CUdepth2<Avgdepth || max(max(cord))==0
                    I1=find(r==j0);
                    L0=[];L1=[];H0=[];H1=[];
                    for ii=1:length(I1)
                        H0(ii)=c(I1(ii));
                        L0=[L0;j0];
                    end
                    [Avgdepth2,CUdepth2]=CUDEPTH(a,b,L0,H0);
                    if p1==1
                        cord(j0-8,:)=0;
                    end
                    j0=j0+8;
                    p1=1;
                end
                p1=0;
                while CUdepth3<Avgdepth || max(max(cord))==0
                    I2=find(r==j1);
                    L0=[];L1=[];H0=[];H1=[];
                    for iii=1:length(I2)
                        H1(iii)=c(I2(iii));
                        L1=[L1;j1];
                    end
                    [Avgdepth3,CUdepth3]=CUDEPTH(a,b,L1,H1);
                    if p1==1
                        cord(j1+8,:)=0;
                    end
                    j1=j1-8;
                    p1=1;
                end
                
                j0=min(min(c));
                j1=max(max(c));
                p1=0;
                while CUdepth4<Avgdepth || max(max(cord))==0
                    I1=find(c==j0);
                    L0=[];L1=[];H0=[];H1=[];
                    for ii=1:length(I1)
                        H0(ii)=r(I1(ii));
                        L0=[L0;j0];
                    end
                    [Avgdepth4,CUdepth4]=CUDEPTH(a,b,H0,L0);
                    if p1==1
                        cord(:,j0-8)=0;
                    end
                    j0=j0+8;
                    p1=1;
                end
                
                p1=0;
                while CUdepth5<Avgdepth || max(max(cord))==0
                    I2=find(c==j1);
                    L0=[];L1=[];H0=[];H1=[];
                    for iii=1:length(I2)
                        H1(iii)=r(I2(iii));
                        L1=[L1;j1];
                    end
                    [Avgdepth5,CUdepth5]=CUDEPTH(a,b,H1,L1);
                    if p1==1
                        cord(:,j1+8)=0;
                    end
                    j1=j1-8;
                    p1=1;
                end
                
                [r2,c2] = find(cord == 1);
                
                for ii=1:length(r2)
                    R=r2(ii);
                    C=c2(ii);
                    S1(R:R+8,C:C+8)=1;
                end
            end
            s2=regionprops(im2bw(S1),'BoundingBox');
            s=[];
            for i=1:length(s2)
                s1=s2(i).BoundingBox;
                s=[s;s1];
            end
            %% Moving object classification in HEVC compressed domain
            if ~isempty(s)
                feature1=[];
                s3=[];
                for i=1:size(s,1)
                    S=s(i,:);
                    blockR=4;
                    clockC=4;
                    a=imcrop(rgb2gray(a1),S);
                    b=imcrop(rgb2gray(b1),S);
                    c=imcrop(rgb2gray(c3),S);
                    feat = extractHOGFeatures (b);
                    BW = edge(b);
                    L = bwlabel(BW,8);
                    numb=max(max(L));
                    S1=zeros(256,320);
                    JJ=zeros(1,numb);
                    for iii=1:numb
                        [r2,c2] = find(L == iii);
                        S9=zeros(256,320);
                        for iiii=1:length(r2)
                            R=r2(iiii);
                            C=c2(iiii);
                            S9(R,C)=1;
                        end
                        JJ(iii)=nnz(S9);
                    end
                    if isempty(JJ)
                        JJ=0;
                    end
                    
                    if outp==1 && max(JJ)>20
                        b1 = insertText(b1,[S(1) S(2)],'Person','FontSize',8);
                        s3=[s3;S];
                    elseif outp==0 && max(JJ)>20
                        b1 = insertText(b1,[S(1) S(2)],'Vehicle','FontSize',8);
                        s3=[s3;S];
                        %                             b1 = insertMarker(b1, Centroids, 'o', 'Size', 6, 'Color', 'r');
                        %                             RGB = insertObjectAnnotation(b1,'rectangle',S,'Vehicle');
                    end
                    
                end
                b1 = insertShape(b1,'rectangle',s3,'LineWidth',1);
                imwrite(b1,['Output Frame\',sprintf('%06d.jpg',j)]);
                %                 drawnow;
                %                 imshow(b1);
            else
                imwrite(b1,['Output Frame\',sprintf('%06d.jpg',j)]);
                %                 drawnow;
                %                 imshow(b1);
            end
        else
            imwrite(b1,['Output Frame\',sprintf('%06d.jpg',j)]);
            %             drawnow;
            %             imshow(b1);
        end
    else
        imwrite(b1,['Output Frame\',sprintf('%06d.jpg',j)]);
        %         drawnow;
        %         imshow(b1);
    end
    
end
% load feature.mat
% feature2=[feature2;feat];
% save feature.mat feature2
files = dir('Output Frame\*.jpg');
aviobj = VideoWriter('Output1.avi'); %creating a movie object
open(aviobj);
for i=1:numel(files) %number of images to be read
    b = fullfile('Output Frame',files(i).name);
    a = imread(b);
    a = uint8(a);%convert the images into unit8 type
    M = im2frame(a);%convert the images into frames
    writeVideo(aviobj,M);
    %     fprintf('adding frame = %i\n', i);
end
disp('Closing movie file...')
close(aviobj);
%% accuracy predictinon
final_accuracy();
%disp('accuracy of the video....')