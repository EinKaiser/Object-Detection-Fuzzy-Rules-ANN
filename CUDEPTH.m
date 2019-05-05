function [Avgdepth,CUdepth]=CUDEPTH(a,b,r,c)

CUdepth=0;
for iii=1:length(r)
    r1=r(iii);
    c1=c(iii);
    a1=a(r1:r1+7,c1:c1+7);
    b1=b(r1:r1+7,c1:c1+7);
    C1=a1-b1;
    MV1=C1(1:2,1:2);
    if sum(sum(MV1))>1
        CUdepth=CUdepth+2;
    end
    MV2=C1(3:4,1:2);
    if sum(sum(MV2))>1
        CUdepth=CUdepth+2;
    end
    MV3=C1(1:2,3:4);
    if sum(sum(MV3))>1
        CUdepth=CUdepth+2;
    end
    MV4=C1(3:4,3:4);
    if sum(sum(MV4))>1
        CUdepth=CUdepth+2;
    end
    MV5=C1(1:2,5:7);
    if sum(sum(MV5))>1
        CUdepth=CUdepth+2;
    end
    MV6=C1(1:2,7:8);
    if sum(sum(MV6))>1
        CUdepth=CUdepth+2;
    end
    MV7=C1(3,5);
    if sum(sum(MV7))>1
        CUdepth=CUdepth+3;
    end
    MV8=C1(3,6);
    if sum(sum(MV8))>1
        CUdepth=CUdepth+3;
    end
    MV9=C1(4,5);
    if sum(sum(MV9))>1
        CUdepth=CUdepth+3;
    end
    MV10=C1(4,6);
    if sum(sum(MV10))>1
        CUdepth=CUdepth+3;
    end
    MV11=C1(3:4,7:8);
    if sum(sum(MV11))>1
        CUdepth=CUdepth+2;
    end
    MV12=C1(5:8,1:4);
    if sum(sum(MV12))>1
        CUdepth=CUdepth+1;
    end
    MV13=C1(5:6,5:6);
    if sum(sum(MV13))>1
        CUdepth=CUdepth+2;
    end
    MV14=C1(5:6,7:8);
    if sum(sum(MV14))>1
        CUdepth=CUdepth+2;
    end
    MV15=C1(7:8,5:6);
    if sum(sum(MV15))>1
        CUdepth=CUdepth+2;
    end
    MV16=C1(7:8,7:8);
    if sum(sum(MV16))>1
        CUdepth=CUdepth+2;
    end
end
Avgdepth=floor(((CUdepth/length(r))+0.5)-1);
end
