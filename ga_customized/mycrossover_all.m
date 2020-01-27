function [parent_1, parent_2] = mycrossover_all(parent1, parent2, nCombo)

dimcol=length(parent1)-1;
mat1=vec2mat(parent1,nCombo);

[m,n]=size(mat1);
mut_num = round(m/4);
%mut_num=2;
pos = sort(randperm(dimcol,2*mut_num));
z=1;
while z<=mut_num
    temp=parent1(1,pos(z):pos(z+1));
    parent1(pos(z):pos(z+1))=parent2(pos(z):pos(z+1));
    parent2(pos(z):pos(z+1))=temp;
    z=z+2;
end
parent_1=parent1;
parent_2=parent2;
