function [ProbotX,ProbotY,dg,sg]=NBAS(ProbotX,ProbotY,XI,YI,ZI,dg,sg)
dir=rands(2,1);%随机一个向量[x，y]的大小
dirm=dir/norm(dir);%向量[x，y]取模
clx1=-dirm(1,1)*dg*1+ProbotX;
cly1=-dirm(2,1)*dg*0.5+ProbotY;%dg*0.5
crx1=dirm(1,1)*dg*1+ProbotX;
cry1=dirm(2,1)*dg*0.5+ProbotY;%dg*0.5

cleft=nongdu(clx1,cly1,XI,YI,ZI);
cright=nongdu(crx1,cry1,XI,YI,ZI);
%if C0<cleft || C0<cright 
BA=sign(cleft-cright);
ProbotX=ProbotX-2*BA*sg*dirm(1,1)-0.5;%-3
ProbotY=ProbotY-BA*sg*dirm(2,1);
dg=0.995*dg+0.01;%触须长度，即传感器离中心距离
sg=0.998*sg;%步长
%end