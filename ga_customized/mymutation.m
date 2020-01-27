function [parent] = mymutation(parent1)

dimCol=length(parent1);
prob_mut=.45;

num=round(prob_mut*dimCol);
k=randperm(dimCol);
k1=k(1:num);

for i=1:num
    if parent1(k1(i)) == 1
        parent1(k1) = 0;
    elseif parent1(k1(i)) == 0
        parent1(k1) = 1;
    end
end
parent=parent1;
