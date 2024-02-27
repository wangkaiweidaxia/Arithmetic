%function [ProbotX1,ProbotY1,ProbotX2,ProbotY2,ProbotX01,ProbotY01,ProbotX02,ProbotY02,C01,C02,dg1,sg1,gB,gBx,gBy,p1B,p1Bx,p1By,p2B,p2Bx,p2By]=BSO(ProbotX1,ProbotY1,ProbotX2,ProbotY2,ProbotX01,ProbotY01,ProbotX02,ProbotY02,C01,C02,XI,YI,ZI,kz1,dg1,sg1,gB,gBx,gBy,p1B,p1Bx,p1By,p2B,p2Bx,p2By)
function[ProbotX1,ProbotY1,ProbotX2,ProbotY2,ProbotX3,ProbotY3,ProbotX4,ProbotY4,ProbotX5,ProbotY5,ProbotX01,ProbotY01,ProbotX02,ProbotY02,ProbotX03,ProbotY03,ProbotX04,ProbotY04,ProbotX05,ProbotY05,C01,C02,C03,C04,C05,dg1,sg1,gB,gBx,gBy,p1B,p1Bx,p1By,p2B,p2Bx,p2By,p3B,p3Bx,p3By,p4B,p4Bx,p4By,p5B,p5Bx,p5By]=BSO(ProbotX1,ProbotY1,ProbotX2,ProbotY2,ProbotX3,ProbotY3,ProbotX4,ProbotY4,ProbotX5,ProbotY5,ProbotX01,ProbotY01,ProbotX02,ProbotY02,ProbotX03,ProbotY03,ProbotX04,ProbotY04,ProbotX05,ProbotY05,C01,C02,C03,C04,C05,XI,YI,ZI,kz1,dg1,sg1,gB,gBx,gBy,p1B,p1Bx,p1By,p2B,p2Bx,p2By,p3B,p3Bx,p3By,p4B,p4Bx,p4By,p5B,p5Bx,p5By)
cc1=0.1;%����ϵ��
cc2=0.1;%����ϵ��
dg1=0.995*dg1+0.01;%���볤�ȣ��������������ľ���
sg1=0.995*sg1;%����

%��1ֻ��ţ�˶�
ProbotX01=ProbotX1;
ProbotY01=ProbotY1;
dir1=rands(2,1);%���һ������[x��y]�Ĵ�С
dirm1=dir1/norm(dir1);%����[x��y]ȡģ
clx1=-dirm1(1,1)*dg1*1+ProbotX1;
cly1=-dirm1(2,1)*dg1*0.5+ProbotY1;
crx1=dirm1(1,1)*dg1*1+ProbotX1;
cry1=dirm1(2,1)*dg1*0.5+ProbotY1;
cleft1=nongdu(clx1,cly1,XI,YI,ZI);
cright1=nongdu(crx1,cry1,XI,YI,ZI);
BA=sign(cleft1-cright1);
ProbotX1=ProbotX1-2*BA*sg1*dirm1(1,1)-5+cc1*rand*(gBx-ProbotX1)+cc2*rand*(p1Bx-ProbotX1);
ProbotY1=ProbotY1-BA*sg1*dirm1(2,1)-1+cc1*rand*(gBy-ProbotY1)+cc2*rand*(p1By-ProbotY1);
C001=C01;%��֮ǰ���Ũ�ȸ���C001��C01������ȡ��ǰŨ��

[C01,ProbotX1,ProbotY1]= bnongdu(ProbotX1,ProbotY1,ProbotX01,ProbotY01,XI,YI,ZI,C001,kz1); %����1��ǰŨ��
if(C01>p1B)
    p1B=C01;
    p1Bx=ProbotX1;
    p1By=ProbotY1;
end

%��2ֻ��ţ�˶�
ProbotX02=ProbotX2;
ProbotY02=ProbotY2;
dir2=rands(2,1);%���һ������[x��y]�Ĵ�С
dirm2=dir2/norm(dir2);%����[x��y]ȡģ
clx2=-dirm2(1,1)*dg1*1+ProbotX2;
cly2=-dirm2(2,1)*dg1*0.5+ProbotY2;
crx2=dirm2(1,1)*dg1*1+ProbotX2;
cry2=dirm2(2,1)*dg1*0.5+ProbotY2;
cleft2=nongdu(clx2,cly2,XI,YI,ZI);
cright2=nongdu(crx2,cry2,XI,YI,ZI);
BA=sign(cleft2-cright2);
ProbotX2=ProbotX2-2*BA*sg1*dirm2(1,1)-5+cc1*rand*(gBx-ProbotX2)+cc2*rand*(p2Bx-ProbotX2);
ProbotY2=ProbotY2-BA*sg1*dirm2(2,1)-1+cc1*rand*(gBy-ProbotY2)+cc2*rand*(p2By-ProbotY2);
C002=C02;%��֮ǰ���Ũ�ȸ���C002��C02������ȡ��ǰŨ��
[C02,ProbotX2,ProbotY2]= bnongdu(ProbotX2,ProbotY2,ProbotX02,ProbotY02,XI,YI,ZI,C002,kz1); %����1��ǰŨ��
if(C02>p2B)
    p2B=C02;
    p2Bx=ProbotX2;
    p2By=ProbotY2;
end

%��3ֻ��ţ�˶�
ProbotX03=ProbotX3;
ProbotY03=ProbotY3;
dir3=rands(2,1);%���һ������[x��y]�Ĵ�С
dirm3=dir3/norm(dir3);%����[x��y]ȡģ
clx3=-dirm3(1,1)*dg1*1+ProbotX3;
cly3=-dirm3(2,1)*dg1*0.5+ProbotY3;
crx3=dirm3(1,1)*dg1*1+ProbotX3;
cry3=dirm3(2,1)*dg1*0.5+ProbotY3;
cleft3=nongdu(clx3,cly3,XI,YI,ZI);
cright3=nongdu(crx3,cry3,XI,YI,ZI);
BA=sign(cleft3-cright3);
ProbotX3=ProbotX3-2*BA*sg1*dirm3(1,1)-5+cc1*rand*(gBx-ProbotX3)+cc2*rand*(p3Bx-ProbotX3);
ProbotY3=ProbotY3-BA*sg1*dirm3(2,1)-1+cc1*rand*(gBy-ProbotY3)+cc2*rand*(p3By-ProbotY3);
C003=C03;%��֮ǰ���Ũ�ȸ���C001��C01������ȡ��ǰŨ��
[C03,ProbotX3,ProbotY3]= bnongdu(ProbotX3,ProbotY3,ProbotX03,ProbotY03,XI,YI,ZI,C003,kz1); %����1��ǰŨ��
if(C03>p3B)
    p3B=C03;
    p3Bx=ProbotX3;
    p3By=ProbotY3;
end
%��4ֻ��ţ�˶�
ProbotX04=ProbotX4;
ProbotY04=ProbotY4;
dir4=rands(2,1);%���һ������[x��y]�Ĵ�С
dirm4=dir4/norm(dir4);%����[x��y]ȡģ
clx4=-dirm4(1,1)*dg1*1+ProbotX4;
cly4=-dirm4(2,1)*dg1*0.5+ProbotY4;
crx4=dirm4(1,1)*dg1*1+ProbotX4;
cry4=dirm4(2,1)*dg1*0.5+ProbotY4;
cleft4=nongdu(clx4,cly4,XI,YI,ZI);
cright4=nongdu(crx4,cry4,XI,YI,ZI);
BA=sign(cleft4-cright4);
ProbotX4=ProbotX4-2*BA*sg1*dirm4(1,1)-5+cc1*rand*(gBx-ProbotX4)+cc2*rand*(p4Bx-ProbotX4);
ProbotY4=ProbotY4-BA*sg1*dirm4(2,1)-1+cc1*rand*(gBy-ProbotY4)+cc2*rand*(p4By-ProbotY4);
C004=C04;%��֮ǰ���Ũ�ȸ���C001��C01������ȡ��ǰŨ��
[C04,ProbotX4,ProbotY4]= bnongdu(ProbotX4,ProbotY4,ProbotX04,ProbotY04,XI,YI,ZI,C004,kz1); %����1��ǰŨ��
if(C04>p4B)
    p4B=C04;
    p4Bx=ProbotX4;
    p4By=ProbotY4;
end
%��5ֻ��ţ�˶�
ProbotX05=ProbotX5;
ProbotY05=ProbotY5;
dir5=rands(2,1);%���һ������[x��y]�Ĵ�С
dirm5=dir5/norm(dir5);%����[x��y]ȡģ
clx5=-dirm5(1,1)*dg1*1+ProbotX5;
cly5=-dirm5(2,1)*dg1*0.5+ProbotY5;
crx5=dirm5(1,1)*dg1*1+ProbotX5;
cry5=dirm5(2,1)*dg1*0.5+ProbotY5;
cleft5=nongdu(clx5,cly5,XI,YI,ZI);
cright5=nongdu(crx5,cry5,XI,YI,ZI);
BA=sign(cleft5-cright5);
ProbotX5=ProbotX5-2*BA*sg1*dirm5(1,1)-5+cc1*rand*(gBx-ProbotX5)+cc2*rand*(p5Bx-ProbotX5);
ProbotY5=ProbotY5-BA*sg1*dirm5(2,1)-1+cc1*rand*(gBy-ProbotY5)+cc2*rand*(p5By-ProbotY5);
C005=C05;%��֮ǰ���Ũ�ȸ���C001��C01������ȡ��ǰŨ��
[C05,ProbotX5,ProbotY5]= bnongdu(ProbotX5,ProbotY5,ProbotX05,ProbotY05,XI,YI,ZI,C005,kz1); %����1��ǰŨ��
if(C05>p5B)
    p5B=C05;
    p5Bx=ProbotX5;
    p5By=ProbotY5;
end
%��ȡȫ������λ��
%gB=max(gB,C01,C02);%���ȫ�����ֵ
if gB<C01;
    gB=C01;%���ȫ�����ֵ��������
    gBx=ProbotX1;
    gBy=ProbotY1;
end
if gB<C03
    gB=C03;
    gBx=ProbotX3;
    gBy=ProbotY3;
end
if gB<C04
    gB=C04;
    gBx=ProbotX4;
    gBy=ProbotY4;
end
if gB<C05
    gB=C05;
    gBx=ProbotX5;
    gBy=ProbotY5;
end
if gB<C02
    gB=C02;
    gBx=ProbotX2;
    gBy=ProbotY2;
end
end


