clear, clc
w1= (0:0.1:1);
w2= (0:0.1:1);
w3= (0:0.1:1);
[W1, W2, W3] = ndgrid(w1, w2, w3);
combs = [W1(:), W2(:), W3(:)];
combs(:,4)=combs(:,1)+combs(:,2)+combs(:,3);

combs=combs(combs(:,4)==1,1:3);
