function [Out, Costs] = cc_outliers(A,k,l,Nx,Ny,Qx,Qy,Dnz,isSelfGraph)
% How much decrease in cost do I get by removing an edge. Return the top
% few such edges.

HowManyOutliers = 50;
Out = zeros(2,HowManyOutliers);
Costs = zeros(1,HowManyOutliers);

for i=1:HowManyOutliers
  BestCostFall = -Inf;
  Bestm1 = -1;
  Bestm2 = -1;
  for m1=1:k
    for m2=1:l
      sizerect = Nx(m1)*Ny(m2);
      num1 = Dnz(m1,m2);
      num2 = sizerect - num1;
      if(Dnz(m1,m2) < 0.5*sizerect)
        num3 = num1-1;
      else
        num3 = num1+1;
      end
      num4 = sizerect-num3;
      if(num1>0 && num1<sizerect)  % Else no possibility of improvement
        CostFall = num1*log(sizerect/num1)+num2*log(sizerect/num2);
	if(num3>0)
          CostFall = CostFall - num3*log(sizerect/num3);
	end
	if(num4>0)
          CostFall = CostFall - num4*log(sizerect/num4);
        end
        if(CostFall > BestCostFall)
          BestCostFall = CostFall;
	  Bestm1 = m1;
	  Bestm2 = m2;
	  Bestnum3 = num3;
        end
      end
    end
  end

  if(Bestm1>0 && Bestm2>0)
    Out(1,i) = Bestm1;
    Out(2,i) = Bestm2;
    Costs(i) = BestCostFall;
    Dnz(Bestm1,Bestm2) = Bestnum3;
  end
end 

