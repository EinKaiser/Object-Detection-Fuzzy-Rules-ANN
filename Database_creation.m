clc;
clear;
close all;
%% get image
[f,d]=uigetfile('Frame\*.jpg');
a1=imresize(imread([d f]),[256 320]);
[f,d]=uigetfile('Frame\*.jpg');
b1=imresize(imread([d f]),[256 320]);
[f,d]=uigetfile('Frame\*.jpg');
c3=imresize(imread([d f]),[256 320]);
a=(rgb2gray(a1));
b=((rgb2gray(b1)));
% c=a1-b1;
% figure,imshow(c);
%% Preprocessing
% BW = edge(b);

[r1,c1]=size(a);
blockSizeR=r1/4;
blockSizeC=c1/5;
sliceNumber=1;
coordinates=[];
coor=zeros(256,320);
figure;
hold on
for rowi = 1 : blockSizeR : r1
    for coli = 1 : blockSizeC : c1
        row1 = rowi;
        row2 = row1 + blockSizeR - 1;
        col1 = coli;
        col2 = col1 + blockSizeC - 1;
        % Extract out the block into a single subimage.
        mainimage = a(row1:row2, col1:col2);
        referenceimage = b(row1:row2, col1:col2);
        for k=1:8:length(mainimage)
            for l=1:8:length(referenceimage)
                ro1=k;
                ro2=ro1+7;
                co1=l;
                co2=co1+7;
                mainpatch=mainimage(ro1:ro2,co1:co2);
                referencepatch=referenceimage(ro1:ro2,co1:co2);
                Mv1=mainpatch(3:4,3:4);
                Mv2=mainpatch(4,5);
                Mv3=mainpatch(4,6);
                Mv4=mainpatch(4,7:8);
                Mv5=mainpatch(5:8,3:4);
                Mv6=mainpatch(5:6,7:8);
                Mv7=mainpatch(7:8,7:8);
                MV1=referencepatch(3:4,3:4);
                MV2=referencepatch(4,5);
                MV3=referencepatch(4,6);
                MV4=referencepatch(4,7:8);
                MV5=referencepatch(5:8,3:4);
                MV6=referencepatch(5:6,7:8);
                MV7=referencepatch(7:8,7:8);
                M1=floor((sum(sum(Mv1-MV1)))/9);
                M2=floor((sum(sum(Mv2-MV2)))/9);
                M3=floor((sum(sum(Mv3-MV3)))/9);
                M4=floor((sum(sum(Mv4-MV4)))/9);
                M5=floor((sum(sum(Mv5-MV5)))/9);
                M6=floor((sum(sum(Mv6-MV6)))/9);
                M7=floor((sum(sum(Mv7-MV7)))/9);
                Motionvector=[M1 M2 M3 M4 M5 M6 M7];
                num=nnz(Motionvector);
                M=Motionvector(Motionvector~=0);
                R = spones(Motionvector);
                if num~=0
                    ind=find(R==1);
                else
                    ind=[];
                end
                zz=[];
                if ~isempty(ind)
                    for i=1:length(ind)
                        in1=ind(i);
%                         in2=0;
%                         in3=0;
                        if i~=length(ind)
                            in2=ind(i+1);
                        else
                            in2=ind(1);
                            if in2==1
                                in2=2;
                            end
                        end
                        if i~=1
                            in3=ind(i-1);
                        else
                            in3=ind(length(ind));
                            if in3==7
                                in3=6;
                            end
                        end
                        if (in1==(in2-1))||(in1==(in3+1))
                            zz=[zz in1];
                        end
                    end
                end
                if length(zz)>1 && (M1==0 || M2==0 || M3==0 || M4==0 || M5==0 || M6==0 || M7==0)
                    coordinates=[coli-1+l rowi-1+k 8 8; coordinates];
                    coor(rowi-1+k,coli-1+l)=1;
                end
            end
        end
        mainimage3D(:, :, sliceNumber) = mainimage;
        referenceimage3D(:, :, sliceNumber) = referenceimage;
        sliceNumber = sliceNumber + 1;
    end
end

[ssx,ssy]=find(coor==1);
coordinating=[];
for i=1:length(ssx)
    x=ssx(i);
    y=ssy(i);
    z=0;
    if x<246
        x1x=x+8;
        x1y=y;
        if coor(x,y)==coor(x1x,x1y)
            z=z+1;
        end
    else
        z=z+0;
    end
    if x>9
        x3x=x-8;
        x3y=y;
        if coor(x,y)==coor(x3x,x3y)
            z=z+1;
        end
    else
        z=z+0;
    end
    if y<310
        x2x=x;
        x2y=y+8;
        if coor(x,y)==coor(x2x,x2y)
            z=z+1;
        end
    else
        z=z+0;
    end
    if y>9
        x4x=x;
        x4y=y-8;
        if coor(x,y)==coor(x4x,x4y)
            z=z+1;
        end
    else
        z=z+0;
    end
    if z>1
        coordinating=[coordinating;y x 8 8];
    end
end
RGB = insertShape(b1,'rectangle',coordinating,'LineWidth',1);
imshow(RGB);
hold off
coord=forward(b1,c3);
%% Moving object segmentation

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
RGB = insertShape(b1,'rectangle',coordination,'LineWidth',1);
figure;
imshow(RGB);
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
RGB = insertShape(b1,'rectangle',coordinate2,'LineWidth',1);
figure;
imshow(RGB);
% object Boundary Refinement
backwa=zeros(256,320);
for i=1:size(coordinate2,1)
    c=coordinate2(i,1);
    r=coordinate2(i,2);
    backwa(r:r+7,c:c+7)=1;
end
figure,imshow(backwa);
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
RGB = insertShape(b1,'rectangle',coordinate,'LineWidth',1);
figure;
imshow(RGB);
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
RGB = insertShape(b1,'rectangle',s,'LineWidth',1);
figure;
imshow(RGB);
%% Moving object classification in HEVC compressed domain

feature=[];
for i=1:size(s,1)
    S=s(i,:);
    blockR=4;
    clockC=4;
    a=imcrop(rgb2gray(a1),S);
    b=imcrop(rgb2gray(b1),S);
    c=imcrop(rgb2gray(c3),S);
    for rowi = 1 : blockR : (round(S(4)/4))*4
        for coli = 1 : clockC : (round(S(3)/4))*4
            row1 = rowi;
            row2 = row1 + blockR - 1;
            col1 = coli;
            col2 = col1 + clockC - 1;
            % Extract out the block into a single subimage.
            currentback = a(row1:row2, col1:col2);
            current = b(row1:row2, col1:col2);
            currentforwa = c(row1:row2, col1:col2);
            len=current-currentforwa;
            MVlength=nnz(len);
            maxcurrentforwa=max(max(current-currentforwa));
            maxcurrentbackwa=max(max(currentback-current));
            feature=[feature;MVlength maxcurrentforwa maxcurrentbackwa];
        end
    end
    if size(feature,1)<25
        continue;
    end
    [k,C]=kmeans(double(feature),25);
    feat=zeros(1,25);
    for ii=1:size(feature,1)
        f=double(feature(ii,:));
        dist=zeros(1,25);
        for iii=1:length(C);
            C1=C(iii,:);
            distance=norm(f-C1);
            dist(iii)=distance;
        end
        minimum=find(dist==min(dist));
        feat(minimum)=feat(minimum)+1;
    end
end
load feature.mat
f1=feature2;
load feature1.mat
f2=feature2;
load feature2.mat
f3=feature2;
f=[f1;f2;f3];
save feat.mat f
feature2=[feature2;feat];
save feat.mat f
% % % 
% % % files = dir('Dataset\pedestrians\input\*.jpg');
% % % aviobj = VideoWriter('pedestrians.avi'); %creating a movie object
% % % open(aviobj);
% % % for i=1:numel(files) %number of images to be read
% % %     b = fullfile('Dataset','pedestrians','input',files(i).name);
% % %     a = imread(b);
% % %     a = uint8(a);%convert the images into unit8 type
% % %     M = im2frame(a);%convert the images into frames
% % % %     aviobj = getframe(M);%add the frames to the avi object created previously
% % %     writeVideo(aviobj,M);
% % %     fprintf('adding frame = %i\n', i);
% % % end
% % % disp('Closing movie file...')
% % % close(aviobj);
% % % disp('Playing movie file...')
% % % implay('pedestrians.avi');