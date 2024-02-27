function [ProbotX,ProbotY,dg,sg]=BAS(ProbotX,ProbotY,XI,YI,ZI,dg,sg)
dir=rands(2,1);
dirm=dir/norm(dir);
clx1=-dirm(1,1)*dg*1+ProbotX;
cly1=-dirm(2,1)*dg*1+ProbotY;
crx1=dirm(1,1)*dg*1+ProbotX;
cry1=dirm(2,1)*dg*1+ProbotY;

cleft=nongdu(clx1,cly1,XI,YI,ZI);
cright=nongdu(crx1,cry1,XI,YI,ZI);
BA=sign(cleft-cright);
ProbotX=ProbotX-BA*sg*dirm(1,1);
ProbotY=ProbotY-BA*sg*dirm(2,1);
dg=0.995*dg+0.01;
sg=0.995*sg;