%%**水污染溯源背景浓度场**%%
% 程序最后的X_R1, Y_R1, C1是水污染溯源所需要的数据，X_R1和Y_R1的第一行和最后一行为岸线坐标，C1为X_R1和Y_R1对应点的浓度。
clear
clc
close all;
data_o=load('qiantangjiang.txt');% 导入河道中心点坐标
data_p=[data_o(:,1)-min(data_o(:,1)),data_o(:,2)-min(data_o(:,2))];%x和y坐标归零
L0=max(data_p(:,1))-min(data_p(:,1));%x方向总长度
xspace=100; % x方向间距
yspace=25;  % y方向间距
L=round(L0-round(mod(L0,xspace)));%取x方向的xsapce的整数倍
Q=1000; %污染物排放速率 单位g/s;
%Q=10*1000*1000/7200; %污染物排放速率 单位g/s;
dt=600;% 单位释放时间，连续排放可以看作很多个瞬时源强组成 单位s
M=Q*dt%  单位时间污染物释放量，源强 ，单位g 
b=10;       % 排放位置距离岸边的距离，m
B=1000;      % 河流宽度,m
Dx=50;      % 横向扩散参数 m2/s
Dy=0.5;     % 纵向扩散参数 m2/s
con_t=3600*3; % 计算时间长度 s
h=3;        % 水深，m
u=0.5;      % 纵向流速,m/s
v=0.015;        % 横向流速,m/s
kc=0.1/24/3600;        % 降解系数
x=0:xspace:L;  % x方向网格
y=-b:yspace:B-b; % y方向网格
Cs=6;      %  水质超标标准
n=2;           % 岸边反射次数
ncol=L/xspace+1;% x方向网格数量
nrow=B/yspace+1;% y方向网格数量
ch=4;%本底浓度

for i=1:ncol
    X1_temp=data_p(:,1);
    Y1_temp=data_p(:,2);
    for k=1:length(X1_temp)
      if x(i)>=X1_temp(k) && x(i)<=X1_temp(k+1) 
          yi(i)=(Y1_temp(k+1)-Y1_temp(k))/(X1_temp(k+1)-X1_temp(k))*(x(i)-X1_temp(k))+Y1_temp(k);
          break;%根据实际河道的x坐标去插值网格点坐标x(i)所对应的y(i)
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

[X_P_Temp,Y_P_Temp]=meshgrid(y);%转换后的Y坐标分量
Y_P=Y_P_Temp;
for i=2:1:ncol
    dis(1)=0; 
    ay=(Y_R(floor(nrow/2+1)+1,i)-Y_R(floor(nrow/2+1)+1,i-1)).^2;
    ax=(X_R(1,i)-X_R(1,i-1)).^2;
    dis(i)=sqrt(ay+ax);
    X_P(1:nrow,1)=X_R(1:nrow,1);
    X_P(1:nrow,i)=X_P(1:nrow,i-1)+dis(i);
end

% 真源浓度计算
con_k=1;
for t=dt:dt:con_t %每次释放时间为dt，释放量为M=Q*dt,计算总时间（排污总时间）为con_t
  for i=1:ncol
    for j=1:nrow
        temp1=M/(4*pi*h*t*sqrt(Dx*Dy));
        temp2=exp(-(X_P(1,i)-u*t)^2/(4*Dx*t));
        temp3=exp(-(Y_P(j,1)-v*t)^2/(4*Dy*t));
        temp4=exp(-kc*t);
        C1(i,j)=temp1*temp2*temp3*temp4;
    end
 end
% 近岸虚源1浓度计算
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
% 远岸虚源2浓度计算
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
  C=C+Ct(:,:,k);% 不同排放时间的浓度场累计
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
[Ch,h]=contour(X_R1,Y_R1,C1);% X_R1, Y_R1, C1是水污染溯源所需要的数据，X_R1和Y_R1的第一行和最后一行为岸线坐标。
colormap cool;
colorbar;
hold on;
plot(X_R1(1,:),Y_R1(1,:),'c','linewidth',2);
plot(X_R1(nrow,:), Y_R1(nrow,:),'c','linewidth',2);
%grid on;
xlabel('x方向距离(m)','fontsize',12);
ylabel('y方向距离(m)','fontsize',12);
%imagesc(X_R,Y_R,C);
%plot(data_p(:,1),data_p(:,2),'*r')

%plot(x,yi,'og');
%plot(X_R,Y_R,'ob')
%plot(X_R1,Y_R1,'.y')

max(max(C))
%xlswrite('qianx.xlsx',X_R)
%xlswrite('qiany.xlsx',Y_R)
%xlswrite('qianc.xlsx',C)
[ x, y ,c]=textread('nongdu10003.csv','%f,%f,%f','headerlines',1);%a为序号 x为x轴 y为y平面 c为浓度
%[ x ,y, c]=textread('qiantangjiangnongdu.txt','%f,%f,%f','headerlines',1);%a为序号 b为y轴 c为z平面 d为x轴 e为浓度
 x=x';
 y=y';
 c=c';
m=length(x)/2;%划分区域
n=length(y)/2;%划分区域
XI=linspace(min(x),max(x),m); %根据需要将X划分为m分
YI=linspace(min(y),max(y),n); %根据需要将Y划分为n分
ZI=griddata(x,y,c,XI,YI.'); %最后ZI是个nxm的矩阵
imagesc(XI,YI,ZI);colorbar; %使用imagesc将矩阵画成图像
plot(X_R1(1,:),Y_R1(1,:),'c','linewidth',1);
plot(X_R1(nrow,:), Y_R1(nrow,:),'c','linewidth',1);
%axis([400 1200,-500 0]);
global ProbotX
global ProbotY

ProbotX=700;ProbotY=-300; %起点坐标
ProbotXD=800;ProbotYD=-350; %对照组起点坐标
%80 30
Nstep=200; %控制步数，以免死循环
Rpx=zeros(1,Nstep); 
Rpy=zeros(1,Nstep);
RpxD=zeros(1,Nstep); 
RpyD=zeros(1,Nstep);
plot(ProbotX,ProbotY,'g*');
hold on;
C0= nongdu(ProbotX,ProbotY ,XI,YI,ZI);
C0D= nongdu(ProbotXD,ProbotYD ,XI,YI,ZI);%对照组浓度
plot(ProbotX,ProbotY,'r.');
ii=0;%V用于判断当前行动模式
ddx=10;%NDG中每一步在X轴上的距离
ddy=10;%NDG中每一步在Y轴上的距离
s1=1;%Six中1逆时针2顺时针
s2=6000;%Six中用于判断
C4=6;%Six中n-1点浓度
C5=5;%Six中n-2点浓度
aa=0;%NSix用于第一步行走，无需改动
dd=10;%NSix算法步长倍率
dg=10;%BAS左右两须初始距离
sg=10;%BAS初始步长
dgD=5;%对照组左右两须初始距离
sgD=5;%对照组初始步长
for kz=1:Nstep
    ProbotX0=ProbotX;%上一点的x坐标
    ProbotY0=ProbotY;%上一点的y坐标
if C0<4%使用搜索算法的浓度阈值
       %[ProbotX,ProbotY]=NEWZ(ProbotX,ProbotY,kz);
        %[ProbotX,ProbotY]=Z(ProbotX,ProbotY,kz);
      [ProbotX,ProbotY,ii]=V(ProbotX,ProbotY,ii,X_R1,Y_R1);
else
       
       %[ProbotX,ProbotY]=DG(ProbotX,ProbotY,XI,YI,ZI);%基础浓度梯度算法
       %[ProbotX,ProbotY]=NDG(ProbotX,ProbotY,XI,YI,ZI,ddx,ddy);%改进变步长浓度梯度算法
       %[ProbotX,ProbotY,s1,s2,C0,C4,C5]=Six(ProbotX,ProbotY,XI,YI,ZI,s1,s2,C0,C4,C5);%六边形算法
       %[ProbotX,ProbotY,s1,s2,C0,C4,C5,aa,dd]=NSix(ProbotX,ProbotY,XI,YI,ZI,s1,s2,C0,C4,C5,aa,dd);%改进六边形算法
       %[ProbotX,ProbotY,dg,sg]=BAS(ProbotX,ProbotY,XI,YI,ZI,dg,sg);%天牛须算法
       [ProbotX,ProbotY,dg,sg]=NBAS(ProbotX,ProbotY,XI,YI,ZI,dg,sg)%改进天牛须算法
end
     
    
    Rpx(kz)=ProbotX; 
    Rpy(kz)=ProbotY; 
    C5=C4;
    C4=C0;%上一点浓度
   %C0= nongdu(ProbotX,ProbotY ,XI,YI,ZI);
     [C0,ProbotX,ProbotY]= bnongdu(ProbotX,ProbotY,ProbotX0,ProbotY0,XI,YI,ZI,C4,kz);%改进浓度，防跳出边界
    CM(1,kz)=C0;
end

for kzD=1:Nstep
    ProbotX0D=ProbotXD;%上一点的x坐标
    ProbotY0D=ProbotYD;%上一点的y坐标
if C0D<4%使用搜索算法的浓度阈值
       %[ProbotX,ProbotY]=NEWZ(ProbotX,ProbotY,kz);
        %[ProbotX,ProbotY]=Z(ProbotX,ProbotY,kz);
      [ProbotXD,ProbotYD,ii]=V(ProbotX,ProbotY,ii,X_R1D,Y_R1D);
else
       
       %[ProbotX,ProbotY]=DG(ProbotX,ProbotY,XI,YI,ZI);%基础浓度梯度算法
       %[ProbotX,ProbotY]=NDG(ProbotX,ProbotY,XI,YI,ZI,ddx,ddy);%改进变步长浓度梯度算法
       %[ProbotX,ProbotY,s1,s2,C0,C4,C5]=Six(ProbotX,ProbotY,XI,YI,ZI,s1,s2,C0,C4,C5);%六边形算法
       %[ProbotX,ProbotY,s1,s2,C0,C4,C5,aa,dd]=NSix(ProbotX,ProbotY,XI,YI,ZI,s1,s2,C0,C4,C5,aa,dd);%改进六边形算法
       [ProbotXD,ProbotYD,dgD,sgD]=BAS(ProbotXD,ProbotYD,XI,YI,ZI,dgD,sgD);%天牛须算法
       %[ProbotXD,ProbotYD,dgD,sgD]=NBAS(ProbotXD,ProbotYD,XI,YI,ZI,dgD,sgD)%改进天牛须算法
end
     
    C4D=4;
    RpxD(kz)=ProbotXD; 
    RpyD(kz)=ProbotYD; 
    C5D=C4D;
    C4D=C0D;%上一点浓度
   %C0= nongdu(ProbotX,ProbotY ,XI,YI,ZI);
     [C0D,ProbotXD,ProbotYD]= bnongdu(ProbotXD,ProbotYD,ProbotX0D,ProbotY0D,XI,YI,ZI,C4D,kzD);%改进浓度，防跳出边界
    CMD(1,kzD)=C0D;
end
%天牛群算法开始
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
% global gB%全局最大值
% global gBx%全局最大值x坐标
% global gBy%全局最大值y坐标
% global p1B%个体1最大值
% global p1Bx%个体1最大值x坐标
% global p1By%个体1最大值y坐标
% global p2B%个体2最大值
% global p2Bx%个体2最大值x坐标
% global p2By%个体2最大值y坐标
% global p3B%个体3最大值
% global p3Bx%个体3最大值x坐标
% global p3By%个体3最大值y坐标
% global p4B%个体4最大值
% global p4Bx%个体4最大值x坐标
% global p4By%个体4最大值y坐标
% global p5B%个体5最大值
% global p5Bx%个体5最大值x坐标
% global p5By%个体5最大值y坐标
% global dg1;
% global sg1;
% dg1=6;%BSO左右两须初始距离
% sg1=12;%BSO初始步长
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
% C01= nongdu(ProbotX1,ProbotY1 ,XI,YI,ZI);%个体1当前浓度
% p1B=C01;
% C001=C01;%个体1上个位置浓度
% C02= nongdu(ProbotX2,ProbotY2 ,XI,YI,ZI);%个体2当前浓度
% p2B=C02;
% C002=C02;%个体2上个位置浓度
% C03= nongdu(ProbotX3,ProbotY3 ,XI,YI,ZI);%个体2当前浓度
% p3B=C03;
% C003=C03;%个体3上个位置浓度
% C04= nongdu(ProbotX4,ProbotY4 ,XI,YI,ZI);%个体4当前浓度
% p4B=C04;
% C004=C04;%个体4上个位置浓度
% C05= nongdu(ProbotX5,ProbotY5 ,XI,YI,ZI);%个体5当前浓度
% p5B=C05;
% C005=C05;%个体5上个位置浓度
% gB=max([C01,C02,C03,C04,C05]);%获得全局最大值
% if gB==C01%获得全局最大值所在坐标
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
%     %[ProbotX1,ProbotY1,ProbotX2,ProbotY2,ProbotX01,ProbotY01,ProbotX02,ProbotY02,C01,C02,dg1,sg1,gB,gBx,gBy,p1B,p1Bx,p1By,p2B,p2Bx,p2By]=BSO(ProbotX1,ProbotY1,ProbotX2,ProbotY2,ProbotX01,ProbotY01,ProbotX02,ProbotY02,C01,C02,XI,YI,ZI,kz1,dg1,sg1,gB,gBx,gBy,p1B,p1Bx,p1By,p2B,p2Bx,p2By);%改进天牛群算法
%     [ProbotX1,ProbotY1,ProbotX2,ProbotY2,ProbotX3,ProbotY3,ProbotX4,ProbotY4,ProbotX5,ProbotY5,ProbotX01,ProbotY01,ProbotX02,ProbotY02,ProbotX03,ProbotY03,ProbotX04,ProbotY04,ProbotX05,ProbotY05,C01,C02,C03,C04,C05,dg1,sg1,gB,gBx,gBy,p1B,p1Bx,p1By,p2B,p2Bx,p2By,p3B,p3Bx,p3By,p4B,p4Bx,p4By,p5B,p5Bx,p5By]=BSO(ProbotX1,ProbotY1,ProbotX2,ProbotY2,ProbotX3,ProbotY3,ProbotX4,ProbotY4,ProbotX5,ProbotY5,ProbotX01,ProbotY01,ProbotX02,ProbotY02,ProbotX03,ProbotY03,ProbotX04,ProbotY04,ProbotX05,ProbotY05,C01,C02,C03,C04,C05,XI,YI,ZI,kz1,dg1,sg1,gB,gBx,gBy,p1B,p1Bx,p1By,p2B,p2Bx,p2By,p3B,p3Bx,p3By,p4B,p4Bx,p4By,p5B,p5Bx,p5By);%改进天牛群算法
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

%天牛群算法结束
plot(Rpx,Rpy,'r-');
plot(ProbotX,ProbotY,'g*');
plot(X_R1(1,:),Y_R1(1,:),'c','linewidth',1);
plot(X_R1(nrow,:), Y_R1(nrow,:),'c','linewidth',1);
%xlabel('浓度 （mg/L）','fontsize',12);
%plot(X_R1(1,:),Y_R1(1,:));
%plot(X_R1(nrow,:), Y_R1(nrow,:));
%plot(X_R1(1,:),Y_R1(1,:),'b*')
%plot(X_R1(nrow,:), Y_R1(nrow,:),'b*');
%axis([0 4000,-500 2000]);
%[maxVal,maxIndex] = max(C0);
figure(2);%开始绘制迭代次数和浓度关系

for cm=1:kz
    if CM(1,cm)==0
        CM(1,cm)=CM(1,cm-1);
    end;
end;
plot(CM(1,:),'k');
hold on;
%对照组绘图开始
for cmD=1:kzD
    if CMD(1,cmD)==0
        CMD(1,cmD)=CMD(1,cmD-1);
    end;
end;
plot(CMD(1,:),'b');

%对照组绘图结束
%grid on;
xlabel('迭代次数','fontsize',12);
ylabel('浓度 （mg/L）','fontsize',12);

%title('改进天牛须搜索算法仿真结果');
legend('改进天牛须算法','标准天牛须算法')
axis([0 500,0 35]);


%天牛群绘图开始
% figure(3);
% subplot(2,3,1);%开始绘制迭代次数和全局最高浓度浓度关系
% for cm=1:kz1
%     if CM0(1,cm)==0
%         CM0(1,cm)=CM0(1,cm-1);
%     end;
% end;
% plot(CM0(1,:));
% grid on;
% xlabel('迭代次数','fontsize',12);
% ylabel('浓度 mg/L','fontsize',12);
% title('天牛群搜索算法最高浓度值仿真结果');
% 
% 
% subplot(2,3,2);%开始绘制个体1迭代次数和浓度关系
% for cm=1:kz1
%     if CM1(1,cm)==0
%         CM1(1,cm)=CM1(1,cm-1);
%     end;
% end;
% plot(CM1(1,:));
% grid on;
% xlabel('迭代次数','fontsize',12);
% ylabel('浓度 mg/L','fontsize',12);
% title('天牛群搜索算法个体1仿真结果');
% 
% 
% subplot(2,3,3);%开始绘制个体2迭代次数和浓度关系
% for cm=1:kz1
%     if CM2(1,cm)==0
%         CM2(1,cm)=CM2(1,cm-1);
%     end;
% end;
% plot(CM2(1,:));
% grid on;
% xlabel('迭代次数','fontsize',12);
% ylabel('浓度 mg/L','fontsize',12);
% title('天牛群搜索算法个体2仿真结果');
% 
% subplot(2,3,4);%开始绘制个体3迭代次数和浓度关系
% for cm=1:kz1
%     if CM3(1,cm)==0
%         CM3(1,cm)=CM3(1,cm-1);
%     end;
% end;
% plot(CM3(1,:));
% grid on;
% xlabel('迭代次数','fontsize',12);
% ylabel('浓度 mg/L','fontsize',12);
% title('天牛群搜索算法个体3仿真结果');
% 
% subplot(2,3,5);%开始绘制个体4迭代次数和浓度关系
% for cm=1:kz1
%     if CM4(1,cm)==0
%         CM4(1,cm)=CM4(1,cm-1);
%     end;
% end;
% plot(CM4(1,:));
% grid on;
% xlabel('迭代次数','fontsize',12);
% ylabel('浓度 mg/L','fontsize',12);
% title('天牛群搜索算法个体4仿真结果');
% 
% subplot(2,3,6);%开始绘制个体5迭代次数和浓度关系
% for cm=1:kz1
%     if CM5(1,cm)==0
%         CM5(1,cm)=CM5(1,cm-1);
%     end;
% end;
% plot(CM5(1,:));
% grid on;
% xlabel('迭代次数','fontsize',12);
% ylabel('浓度 mg/L','fontsize',12);
% title('天牛群搜索算法个体5仿真结果');
%天牛群绘图结束


%以下五行为鲸鱼算法求最高浓度点，有一定概率出错，建议在浓度场构建完后单独运行WOA文件夹中的main.m来执行
% SearchAgents_no=30; % Number of search agents
% Function_name='F24'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
% Max_iteration=500; % Maximum numbef of iterations% Load details of the selected benchmark function
% [lb,ub,dim,fobj]=Get_Functions_details(Function_name);
% [Best_score,Best_pos,WOA_cg_curve]=WOA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj,XI,YI,ZI);