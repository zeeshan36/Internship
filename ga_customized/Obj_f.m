% function [TotalV,Tfvio] = Obj_f(x)
function [fvalue,sum_sj,sum_dui,sum_dli] = Obj_f(x,nShift,nCombo,dli,dui,sj,t)
%Time halo violation cost of using matrix Xij

% nvars=nShift*nCombo;
% b=[-dli;dui;-sj;sj];

% Ar=kron(eye(nShift),ones(1,nCombo));%Summation along rows of Xij can be expressed as a matrix-vector multiplication A1*XijT(:) 
% Ac=kron(ones(1,nShift),eye(nCombo));% Summation along columns can be expressed the same way with

  %Note:  for integer programming, Aeq and beq must be empty ([]), so
  %equality contraints are rewritten to si<=Ac<=sj
% A=[-Ar;Ar;-Ac;Ac];

%Convert the variable vector x to matrices with nCombo columns.
Xij=vec2mat(x,nCombo);
  %Compute cij for the initial xij_Feasible matrix
Cij=[Xij;Xij(1:t-1,:)];%add the remaining t-1 element's cost into matrix Xij
[m,n]=size(Xij);
cij=zeros(m,n);

for j=1:nCombo
    for i=1:nShift
        cij(i,j)=sum(Cij((i:i+t-1),j));
    end
end

% sum_f=0;
value=cij.*Xij;
fvalue=sum(value(:));
% for i=1:nShift
%     for j=1:nCombo
%         sum_f=sum_f+cij(i,j)*Xij(i,j);
%     end
% end
% sum_f
sum_sj=0;
sum_dli=0;
sum_dui=0;

for j=1:nCombo
    if sum(Xij(:,j))~=sj(j)
        sum_sj=sum_sj + abs(sum(Xij(:,j))-sj(j));
    end
end
for i=1:nShift
    if sum(Xij(i,:))>dui(i)
        sum_dui = sum_dui + (sum(Xij(i,:))-dui(i));
    end
    if sum(Xij(i,:))<dli(i)
        sum_dli = sum_dli + (dli(i)-sum(Xij(i,:)));
    end
end

%ConstraintVio=sum(b-A*x');
cons = sum_sj + sum_dui + sum_dli;
%Tfval=fvalue + cons;

% The double summation for Cij and Xij can be epxressed can be expressed C(:).'(X(:)).
end

