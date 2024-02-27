%%**ˮ��Ⱦ��Դ����Ũ�ȳ�**%%
% ��������X_R1, Y_R1, C1��ˮ��Ⱦ��Դ����Ҫ�����ݣ�X_R1��Y_R1�ĵ�һ�к����һ��Ϊ�������꣬C1ΪX_R1��Y_R1��Ӧ���Ũ�ȡ�
clear
clc
close all;
data_o=load('qiantangjiang.txt');% ����ӵ����ĵ�����
data_p=[data_o(:,1)-min(data_o(:,1)),data_o(:,2)-min(data_o(:,2))];%x��y�������
L0=max(data_p(:,1))-min(data_p(:,1));%x�����ܳ���
xspace=100; % x������
yspace=25;  % y������
L=round(L0-round(mod(L0,xspace)));%ȡx�����xsapce��������
Q=1000; %��Ⱦ���ŷ����� ��λg/s;
%Q=10*1000*1000/7200; %��Ⱦ���ŷ����� ��λg/s;
dt=600;% ��λ�ͷ�ʱ�䣬�����ŷſ��Կ����ܶ��˲ʱԴǿ��� ��λs
M=Q*dt%  ��λʱ����Ⱦ���ͷ�����Դǿ ����λg 
b=10;       % �ŷ�λ�þ��밶�ߵľ��룬m
B=1000;      % �������,m
Dx=50;      % ������ɢ���� m2/s
Dy=0.5;     % ������ɢ���� m2/s
con_t=3600*3; % ����ʱ�䳤�� s
h=3;        % ˮ�m
u=0.5;      % ��������,m/s
v=0.015;        % ��������,m/s
kc=0.1/24/3600;        % ����ϵ��
x=0:xspace:L;  % x��������
y=-b:yspace:B-b; % y��������
Cs=6;      %  ˮ�ʳ����׼
n=2;           % ���߷������
ncol=L/xspace+1;% x������������
nrow=B/yspace+1;% y������������
ch=4;%����Ũ��

for i=1:ncol
    X1_temp=data_p(:,1);
    Y1_temp=data_p(:,2);
    for k=1:length(X1_temp)
      if x(i)>=X1_temp(k) && x(i)<=X1_temp(k+1) 
          yi(i)=(Y1_temp(k+1)-Y1_temp(k))/(X1_temp(k+1)-X1_temp(k))*(x(i)-X1_temp(k))+Y1_temp(k);
          break;%����ʵ�ʺӵ���x����ȥ��ֵ���������x(i)����Ӧ��y(i)
      end
    end
end

for i=1:ncol-1
  if (yi(i+1)-yi(i))~=0
        k1=-(x(i+1)-x(i))/(yi(i+1)-yi(i));
        jj=0;
        for j=B/2/yspace:-1:1
           jj=jj+1;
           X_R(jj,i)=x(i)+sqrt((j*yspace)*(j*yspace)/(k1*k1+1));
           Y_R(jj,i)=yi(i)+k1*sqrt((j*yspace)*(j*yspace)/(k1*k1+1));
        end
        for j=0:B/2/yspace
           X_R(j+B/2/yspace+1,i)=x(i)-sqrt((j*yspace)*(j*yspace)/(k1*k1+1));
           Y_R(j+B/2/yspace+1,i)=yi(i)-k1*sqrt((j*yspace)*(j*yspace)/(k1*k1+1));
        end 
  else 
        jj=0;
        for j=B/2/yspace:-1:1
          jj=jj+1;
          X_R(jj,i)=x(i);
          Y_R(jj,i)=yi(i)+j*yspace;
        end
        for j=0:B/2/yspace
          X_R(j+B/2/yspace+1,i)=x(i);
          Y_R(j+B/2/yspace+1,i)=yi(i)-j*yspace;
        end
  end
end
i=ncol  
if(yi(ncol-1)-yi(ncol))~=0
   k2=-(x(ncol-1)-x(ncol))/(yi(ncol-1)-yi(ncol));
   jj=0;
   for j=B/2/yspace:-1:1
     jj=jj+1;
     X_R(jj,i)=x(ncol)+sqrt((j*yspace)*(j*yspace)/(k2*k2+1));
     Y_R(jj,i)=yi(ncol)+k2*sqrt((j*yspace)*(j*yspace)/(k2*k2+1));   
   end

   for j=0:B/2/yspace
     X_R(jj,i)=x(ncol)-sqrt((j*yspace)*(j*yspace)/(k2*k2+1));
     Y_R(jj,i)=yi(ncol)-k2*sqrt((j*yspace)*(j*yspace)/(k2*k2+1));
   end
else 
    jj=0;
    for j=B/2/yspace:-1:1
      jj=jj+1;
      X_R(jj,i)=x(ncol);
      Y_R(jj,i)=yi(ncol)+j*yspace;   
    end
    for j=0:B/2/yspace
      X_R(jj,i)=x(ncol);
      Y_R(jj,i)=yi(ncol)-j*yspace;
    end
end        

[X_P_Temp,Y_P_Temp]=meshgrid(y);%ת�����Y�������
Y_P=Y_P_Temp;
for i=2:1:ncol
    dis(1)=0; 
    ay=(Y_R(floor(nrow/2+1)+1,i)-Y_R(floor(nrow/2+1)+1,i-1)).^2;
    ax=(X_R(1,i)-X_R(1,i-1)).^2;
    dis(i)=sqrt(ay+ax);
    X_P(1:nrow,1)=X_R(1:nrow,1);
    X_P(1:nrow,i)=X_P(1:nrow,i-1)+dis(i);
end

% ��ԴŨ�ȼ���
con_k=1;
for t=dt:dt:con_t %ÿ���ͷ�ʱ��Ϊdt���ͷ���ΪM=Q*dt,������ʱ�䣨������ʱ�䣩Ϊcon_t
  for i=1:ncol
    for j=1:nrow
        temp1=M/(4*pi*h*t*sqrt(Dx*Dy));
        temp2=exp(-(X_P(1,i)-u*t)^2/(4*Dx*t));
        temp3=exp(-(Y_P(j,1)-v*t)^2/(4*Dy*t));
        temp4=exp(-kc*t);
        C1(i,j)=temp1*temp2*temp3*temp4;
    end
 end
% ������Դ1Ũ�ȼ���
 for kk=1:n
  for i=1:ncol
    for j=1:nrow
        temp1=M/(4*pi*h*t*sqrt(Dx*Dy));
        temp2=exp(-(X_P(1,i)-u*t)^2/(4*Dx*t));
        temp3=exp(-(2*kk*b+Y_P(j,1)-v*t)^2/(4*Dy*t));
        temp4=exp(-kc*t);
        C2(i,j)=temp1*temp2*temp3*temp4;
    end
  end
 end
% Զ����Դ2Ũ�ȼ���
 for kk=1:n
  for i=1:ncol
    for j=1:nrow
        temp1=M/(4*pi*h*t*sqrt(Dx*Dy));
        temp2=exp(-(X_P(1,i)-u*t)^2/(4*Dx*t));
        temp3=exp(-(2*kk*B-2*kk*b-Y_P(j,1)-v*t)^2/(4*Dy*t));
        temp4=exp(-kc*t);
        C3(i,j)=temp1*temp2*temp3*temp4;
    end
  end
 end
 Ct(:,:,con_k)=(C1+C2+C3)';
 con_k=con_k+1;
end
C=Ct(:,:,1);

for k=2:1:con_k-1
  C=C+Ct(:,:,k);% ��ͬ�ŷ�ʱ���Ũ�ȳ��ۼ�
end
C=C+ch;
figure('color',[1 1 1]);
%X_R(1,40)=4200;
%X_R(1,41)=4210;
%Y_R(21,70)=6068;
%X_R(1,72)=7324;
%X_R(1,71)=7224;
X_R1=X_R(:,1:ncol-1);
Y_R1=Y_R(:,1:ncol-1);
C1=C(:,1:ncol-1);
%

%
[Ch,h]=contour(X_R1,Y_R1,C1);% X_R1, Y_R1, C1��ˮ��Ⱦ��Դ����Ҫ�����ݣ�X_R1��Y_R1�ĵ�һ�к����һ��Ϊ�������ꡣ
colormap cool;
colorbar;
hold on;
plot(X_R1(1,:),Y_R1(1,:),'c','linewidth',2);
plot(X_R1(nrow,:), Y_R1(nrow,:),'c','linewidth',2);
%grid on;
xlabel('x�������(m)','fontsize',12);
ylabel('y�������(m)','fontsize',12);
%imagesc(X_R,Y_R,C);
%plot(data_p(:,1),data_p(:,2),'*r')

%plot(x,yi,'og');
%plot(X_R,Y_R,'ob')
%plot(X_R1,Y_R1,'.y')

max(max(C))
%xlswrite('qianx.xlsx',X_R)
%xlswrite('qiany.xlsx',Y_R)
%xlswrite('qianc.xlsx',C)
[ x, y ,c]=textread('nongdu10003.csv','%f,%f,%f','headerlines',1);%aΪ��� xΪx�� yΪyƽ�� cΪŨ��
%[ x ,y, c]=textread('qiantangjiangnongdu.txt','%f,%f,%f','headerlines',1);%aΪ��� bΪy�� cΪzƽ�� dΪx�� eΪŨ��
 x=x';
 y=y';
 c=c';
m=length(x)/2;%��������
n=length(y)/2;%��������
XI=linspace(min(x),max(x),m); %������Ҫ��X����Ϊm��
YI=linspace(min(y),max(y),n); %������Ҫ��Y����Ϊn��
ZI=griddata(x,y,c,XI,YI.'); %���ZI�Ǹ�nxm�ľ���
imagesc(XI,YI,ZI);colorbar; %ʹ��imagesc�����󻭳�ͼ��
plot(X_R1(1,:),Y_R1(1,:),'c','linewidth',1);
plot(X_R1(nrow,:), Y_R1(nrow,:),'c','linewidth',1);
%axis([400 1200,-500 0]);
global ProbotX
global ProbotY

ProbotX=700;ProbotY=-300; %�������
ProbotXD=800;ProbotYD=-350; %�������������
%80 30
Nstep=200; %���Ʋ�����������ѭ��
Rpx=zeros(1,Nstep); 
Rpy=zeros(1,Nstep);
RpxD=zeros(1,Nstep); 
RpyD=zeros(1,Nstep);
plot(ProbotX,ProbotY,'g*');
hold on;
C0= nongdu(ProbotX,ProbotY ,XI,YI,ZI);
C0D= nongdu(ProbotXD,ProbotYD ,XI,YI,ZI);%������Ũ��
plot(ProbotX,ProbotY,'r.');
ii=0;%V�����жϵ�ǰ�ж�ģʽ
ddx=10;%NDG��ÿһ����X���ϵľ���
ddy=10;%NDG��ÿһ����Y���ϵľ���
s1=1;%Six��1��ʱ��2˳ʱ��
s2=6000;%Six�������ж�
C4=6;%Six��n-1��Ũ��
C5=5;%Six��n-2��Ũ��
aa=0;%NSix���ڵ�һ�����ߣ�����Ķ�
dd=10;%NSix�㷨��������
dg=10;%BAS���������ʼ����
sg=10;%BAS��ʼ����
dgD=5;%���������������ʼ����
sgD=5;%�������ʼ����
for kz=1:Nstep
    ProbotX0=ProbotX;%��һ���x����
    ProbotY0=ProbotY;%��һ���y����
if C0<4%ʹ�������㷨��Ũ����ֵ
       %[ProbotX,ProbotY]=NEWZ(ProbotX,ProbotY,kz);
        %[ProbotX,ProbotY]=Z(ProbotX,ProbotY,kz);
      [ProbotX,ProbotY,ii]=V(ProbotX,ProbotY,ii,X_R1,Y_R1);
else
       
       %[ProbotX,ProbotY]=DG(ProbotX,ProbotY,XI,YI,ZI);%����Ũ���ݶ��㷨
       %[ProbotX,ProbotY]=NDG(ProbotX,ProbotY,XI,YI,ZI,ddx,ddy);%�Ľ��䲽��Ũ���ݶ��㷨
       %[ProbotX,ProbotY,s1,s2,C0,C4,C5]=Six(ProbotX,ProbotY,XI,YI,ZI,s1,s2,C0,C4,C5);%�������㷨
       %[ProbotX,ProbotY,s1,s2,C0,C4,C5,aa,dd]=NSix(ProbotX,ProbotY,XI,YI,ZI,s1,s2,C0,C4,C5,aa,dd);%�Ľ��������㷨
       %[ProbotX,ProbotY,dg,sg]=BAS(ProbotX,ProbotY,XI,YI,ZI,dg,sg);%��ţ���㷨
       [ProbotX,ProbotY,dg,sg]=NBAS(ProbotX,ProbotY,XI,YI,ZI,dg,sg)%�Ľ���ţ���㷨
end
     
    
    Rpx(kz)=ProbotX; 
    Rpy(kz)=ProbotY; 
    C5=C4;
    C4=C0;%��һ��Ũ��
   %C0= nongdu(ProbotX,ProbotY ,XI,YI,ZI);
     [C0,ProbotX,ProbotY]= bnongdu(ProbotX,ProbotY,ProbotX0,ProbotY0,XI,YI,ZI,C4,kz);%�Ľ�Ũ�ȣ��������߽�
    CM(1,kz)=C0;
end

for kzD=1:Nstep
    ProbotX0D=ProbotXD;%��һ���x����
    ProbotY0D=ProbotYD;%��һ���y����
if C0D<4%ʹ�������㷨��Ũ����ֵ
       %[ProbotX,ProbotY]=NEWZ(ProbotX,ProbotY,kz);
        %[ProbotX,ProbotY]=Z(ProbotX,ProbotY,kz);
      [ProbotXD,ProbotYD,ii]=V(ProbotX,ProbotY,ii,X_R1D,Y_R1D);
else
       
       %[ProbotX,ProbotY]=DG(ProbotX,ProbotY,XI,YI,ZI);%����Ũ���ݶ��㷨
       %[ProbotX,ProbotY]=NDG(ProbotX,ProbotY,XI,YI,ZI,ddx,ddy);%�Ľ��䲽��Ũ���ݶ��㷨
       %[ProbotX,ProbotY,s1,s2,C0,C4,C5]=Six(ProbotX,ProbotY,XI,YI,ZI,s1,s2,C0,C4,C5);%�������㷨
       %[ProbotX,ProbotY,s1,s2,C0,C4,C5,aa,dd]=NSix(ProbotX,ProbotY,XI,YI,ZI,s1,s2,C0,C4,C5,aa,dd);%�Ľ��������㷨
       [ProbotXD,ProbotYD,dgD,sgD]=BAS(ProbotXD,ProbotYD,XI,YI,ZI,dgD,sgD);%��ţ���㷨
       %[ProbotXD,ProbotYD,dgD,sgD]=NBAS(ProbotXD,ProbotYD,XI,YI,ZI,dgD,sgD)%�Ľ���ţ���㷨
end
     
    C4D=4;
    RpxD(kz)=ProbotXD; 
    RpyD(kz)=ProbotYD; 
    C5D=C4D;
    C4D=C0D;%��һ��Ũ��
   %C0= nongdu(ProbotX,ProbotY ,XI,YI,ZI);
     [C0D,ProbotXD,ProbotYD]= bnongdu(ProbotXD,ProbotYD,ProbotX0D,ProbotY0D,XI,YI,ZI,C4D,kzD);%�Ľ�Ũ�ȣ��������߽�
    CMD(1,kzD)=C0D;
end
%��ţȺ�㷨��ʼ
% global ProbotX1
% global ProbotY1
% global ProbotX01
% global ProbotY01
% global ProbotX2
% global ProbotY2
% global ProbotX02
% global ProbotY02
% global ProbotX3
% global ProbotY3
% global ProbotX03
% global ProbotY03
% global ProbotX4
% global ProbotY4
% global ProbotX04
% global ProbotY04
% global ProbotX5
% global ProbotY5
% global ProbotX05
% global ProbotY05
% global gB%ȫ�����ֵ
% global gBx%ȫ�����ֵx����
% global gBy%ȫ�����ֵy����
% global p1B%����1���ֵ
% global p1Bx%����1���ֵx����
% global p1By%����1���ֵy����
% global p2B%����2���ֵ
% global p2Bx%����2���ֵx����
% global p2By%����2���ֵy����
% global p3B%����3���ֵ
% global p3Bx%����3���ֵx����
% global p3By%����3���ֵy����
% global p4B%����4���ֵ
% global p4Bx%����4���ֵx����
% global p4By%����4���ֵy����
% global p5B%����5���ֵ
% global p5Bx%����5���ֵx����
% global p5By%����5���ֵy����
% global dg1;
% global sg1;
% dg1=6;%BSO���������ʼ����
% sg1=12;%BSO��ʼ����
% ProbotX1=1400;ProbotY1=-200;ProbotX01=ProbotX1;ProbotY01=ProbotY1;
% ProbotX2=1200;ProbotY2=-300;ProbotX02=ProbotX2;ProbotY02=ProbotY2;
% ProbotX3=2000;ProbotY3=0;ProbotX03=ProbotX3;ProbotY03=ProbotY3;
% ProbotX4=500;ProbotY4=-300;ProbotX04=ProbotX4;ProbotY04=ProbotY4;
% ProbotX5=400;ProbotY5=-100;ProbotX05=ProbotX5;ProbotY05=ProbotY5;
% plot(ProbotX1,ProbotY1,'b*');
% plot(ProbotX2,ProbotY2,'y*');
% plot(ProbotX3,ProbotY3,'r*');
% plot(ProbotX4,ProbotY4,'g*');
% plot(ProbotX5,ProbotY5,'w*');
% p1Bx=ProbotX1;p1By=ProbotY1;p2Bx=ProbotX2;p2By=ProbotY2;
% p3Bx=ProbotX3;p3By=ProbotY3;p4Bx=ProbotX4;p4By=ProbotY4;p5Bx=ProbotX5;p5By=ProbotY5;
% 
% C01= nongdu(ProbotX1,ProbotY1 ,XI,YI,ZI);%����1��ǰŨ��
% p1B=C01;
% C001=C01;%����1�ϸ�λ��Ũ��
% C02= nongdu(ProbotX2,ProbotY2 ,XI,YI,ZI);%����2��ǰŨ��
% p2B=C02;
% C002=C02;%����2�ϸ�λ��Ũ��
% C03= nongdu(ProbotX3,ProbotY3 ,XI,YI,ZI);%����2��ǰŨ��
% p3B=C03;
% C003=C03;%����3�ϸ�λ��Ũ��
% C04= nongdu(ProbotX4,ProbotY4 ,XI,YI,ZI);%����4��ǰŨ��
% p4B=C04;
% C004=C04;%����4�ϸ�λ��Ũ��
% C05= nongdu(ProbotX5,ProbotY5 ,XI,YI,ZI);%����5��ǰŨ��
% p5B=C05;
% C005=C05;%����5�ϸ�λ��Ũ��
% gB=max([C01,C02,C03,C04,C05]);%���ȫ�����ֵ
% if gB==C01%���ȫ�����ֵ��������
%     gBx=ProbotX1;
%     gBy=ProbotY1;
% end
% if gB==C02
%     gBx=ProbotX2;
%     gBy=ProbotY2;
% end
% if gB==C03
%     gBx=ProbotX3;
%     gBy=ProbotY3;
% end
% if gB==C04
%     gBx=ProbotX4;
%     gBy=ProbotY4;
% end
% if gB==C05
%     gBx=ProbotX5;
%     gBy=ProbotY5;
% end
% for kz1=1:Nstep
%     %[ProbotX1,ProbotY1,ProbotX2,ProbotY2,ProbotX01,ProbotY01,ProbotX02,ProbotY02,C01,C02,dg1,sg1,gB,gBx,gBy,p1B,p1Bx,p1By,p2B,p2Bx,p2By]=BSO(ProbotX1,ProbotY1,ProbotX2,ProbotY2,ProbotX01,ProbotY01,ProbotX02,ProbotY02,C01,C02,XI,YI,ZI,kz1,dg1,sg1,gB,gBx,gBy,p1B,p1Bx,p1By,p2B,p2Bx,p2By);%�Ľ���ţȺ�㷨
%     [ProbotX1,ProbotY1,ProbotX2,ProbotY2,ProbotX3,ProbotY3,ProbotX4,ProbotY4,ProbotX5,ProbotY5,ProbotX01,ProbotY01,ProbotX02,ProbotY02,ProbotX03,ProbotY03,ProbotX04,ProbotY04,ProbotX05,ProbotY05,C01,C02,C03,C04,C05,dg1,sg1,gB,gBx,gBy,p1B,p1Bx,p1By,p2B,p2Bx,p2By,p3B,p3Bx,p3By,p4B,p4Bx,p4By,p5B,p5Bx,p5By]=BSO(ProbotX1,ProbotY1,ProbotX2,ProbotY2,ProbotX3,ProbotY3,ProbotX4,ProbotY4,ProbotX5,ProbotY5,ProbotX01,ProbotY01,ProbotX02,ProbotY02,ProbotX03,ProbotY03,ProbotX04,ProbotY04,ProbotX05,ProbotY05,C01,C02,C03,C04,C05,XI,YI,ZI,kz1,dg1,sg1,gB,gBx,gBy,p1B,p1Bx,p1By,p2B,p2Bx,p2By,p3B,p3Bx,p3By,p4B,p4Bx,p4By,p5B,p5Bx,p5By);%�Ľ���ţȺ�㷨
%     CM0(1,kz1)=gB;
%     CM1(1,kz1)=C01;
%     CM2(1,kz1)=C02;
%     CM3(1,kz1)=C03;
%     CM4(1,kz1)=C04;
%     CM5(1,kz1)=C05;
%     
%     Rpx1(kz1)=ProbotX1; 
%     Rpy1(kz1)=ProbotY1; 
%     Rpx2(kz1)=ProbotX2; 
%     Rpy2(kz1)=ProbotY2; 
%     Rpx3(kz1)=ProbotX3; 
%     Rpy3(kz1)=ProbotY3; 
%     Rpx4(kz1)=ProbotX4; 
%     Rpy4(kz1)=ProbotY4; 
%     Rpx5(kz1)=ProbotX5; 
%     Rpy5(kz1)=ProbotY5; 
% end
% plot(Rpx1,Rpy1,'b-');
% plot(ProbotX1,ProbotY1,'b*');
% plot(Rpx2,Rpy2,'y-');
% plot(ProbotX2,ProbotY2,'y*');
% plot(Rpx3,Rpy3,'r-');
% plot(ProbotX3,ProbotY3,'r*');
% plot(Rpx4,Rpy4,'g-');
% plot(ProbotX4,ProbotY4,'g*');
% plot(Rpx5,Rpy5,'w-');
% plot(ProbotX5,ProbotY5,'w*');
% % 

%��ţȺ�㷨����
plot(Rpx,Rpy,'r-');
plot(ProbotX,ProbotY,'g*');
plot(X_R1(1,:),Y_R1(1,:),'c','linewidth',1);
plot(X_R1(nrow,:), Y_R1(nrow,:),'c','linewidth',1);
%xlabel('Ũ�� ��mg/L��','fontsize',12);
%plot(X_R1(1,:),Y_R1(1,:));
%plot(X_R1(nrow,:), Y_R1(nrow,:));
%plot(X_R1(1,:),Y_R1(1,:),'b*')
%plot(X_R1(nrow,:), Y_R1(nrow,:),'b*');
%axis([0 4000,-500 2000]);
%[maxVal,maxIndex] = max(C0);
figure(2);%��ʼ���Ƶ���������Ũ�ȹ�ϵ

for cm=1:kz
    if CM(1,cm)==0
        CM(1,cm)=CM(1,cm-1);
    end;
end;
plot(CM(1,:),'k');
hold on;
%�������ͼ��ʼ
for cmD=1:kzD
    if CMD(1,cmD)==0
        CMD(1,cmD)=CMD(1,cmD-1);
    end;
end;
plot(CMD(1,:),'b');

%�������ͼ����
%grid on;
xlabel('��������','fontsize',12);
ylabel('Ũ�� ��mg/L��','fontsize',12);

%title('�Ľ���ţ�������㷨������');
legend('�Ľ���ţ���㷨','��׼��ţ���㷨')
axis([0 500,0 35]);


%��ţȺ��ͼ��ʼ
% figure(3);
% subplot(2,3,1);%��ʼ���Ƶ���������ȫ�����Ũ��Ũ�ȹ�ϵ
% for cm=1:kz1
%     if CM0(1,cm)==0
%         CM0(1,cm)=CM0(1,cm-1);
%     end;
% end;
% plot(CM0(1,:));
% grid on;
% xlabel('��������','fontsize',12);
% ylabel('Ũ�� mg/L','fontsize',12);
% title('��ţȺ�����㷨���Ũ��ֵ������');
% 
% 
% subplot(2,3,2);%��ʼ���Ƹ���1����������Ũ�ȹ�ϵ
% for cm=1:kz1
%     if CM1(1,cm)==0
%         CM1(1,cm)=CM1(1,cm-1);
%     end;
% end;
% plot(CM1(1,:));
% grid on;
% xlabel('��������','fontsize',12);
% ylabel('Ũ�� mg/L','fontsize',12);
% title('��ţȺ�����㷨����1������');
% 
% 
% subplot(2,3,3);%��ʼ���Ƹ���2����������Ũ�ȹ�ϵ
% for cm=1:kz1
%     if CM2(1,cm)==0
%         CM2(1,cm)=CM2(1,cm-1);
%     end;
% end;
% plot(CM2(1,:));
% grid on;
% xlabel('��������','fontsize',12);
% ylabel('Ũ�� mg/L','fontsize',12);
% title('��ţȺ�����㷨����2������');
% 
% subplot(2,3,4);%��ʼ���Ƹ���3����������Ũ�ȹ�ϵ
% for cm=1:kz1
%     if CM3(1,cm)==0
%         CM3(1,cm)=CM3(1,cm-1);
%     end;
% end;
% plot(CM3(1,:));
% grid on;
% xlabel('��������','fontsize',12);
% ylabel('Ũ�� mg/L','fontsize',12);
% title('��ţȺ�����㷨����3������');
% 
% subplot(2,3,5);%��ʼ���Ƹ���4����������Ũ�ȹ�ϵ
% for cm=1:kz1
%     if CM4(1,cm)==0
%         CM4(1,cm)=CM4(1,cm-1);
%     end;
% end;
% plot(CM4(1,:));
% grid on;
% xlabel('��������','fontsize',12);
% ylabel('Ũ�� mg/L','fontsize',12);
% title('��ţȺ�����㷨����4������');
% 
% subplot(2,3,6);%��ʼ���Ƹ���5����������Ũ�ȹ�ϵ
% for cm=1:kz1
%     if CM5(1,cm)==0
%         CM5(1,cm)=CM5(1,cm-1);
%     end;
% end;
% plot(CM5(1,:));
% grid on;
% xlabel('��������','fontsize',12);
% ylabel('Ũ�� mg/L','fontsize',12);
% title('��ţȺ�����㷨����5������');
%��ţȺ��ͼ����


%��������Ϊ�����㷨�����Ũ�ȵ㣬��һ�����ʳ���������Ũ�ȳ�������󵥶�����WOA�ļ����е�main.m��ִ��
% SearchAgents_no=30; % Number of search agents
% Function_name='F24'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
% Max_iteration=500; % Maximum numbef of iterations% Load details of the selected benchmark function
% [lb,ub,dim,fobj]=Get_Functions_details(Function_name);
% [Best_score,Best_pos,WOA_cg_curve]=WOA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj,XI,YI,ZI);