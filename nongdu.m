%返回当前点位的浓度
function  C0= nongdu(ProbotX,ProbotY ,XI,YI,ZI)
max0=length(XI);
	for i=1:max0
        if(ProbotX>=XI(i)&&ProbotX<XI(i+1))
            for j=1:max0
                if(ProbotY>=YI(j)&&ProbotY<YI(j+1))
                    C0=ZI(j,i);
                end
            end
        end
  end
end
