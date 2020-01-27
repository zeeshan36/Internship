% function [TotalV,Tfvio] = Obj_f(x)
function [TotalV,Tfvio,cons] = Obj_f(x)
%Time halo violation cost of using matrix Xij

%   Give parameters same as GA_Test.m 
% nShift=60;
% nCombo=150;
% t=5;
% 
global nShift
global nCombo
global t
global dli
global dui
global sj

% nShift=10;
% nCombo=20;
% t=7;
% Intlpopsize=400;
% nvars=nShift*nCombo;
% b=[-dli;dui;-sj;sj];

% filename = '200.xlsx';
% sheet = 'A';
% 
% dli=xlsread(filename,sheet,'A2:A11');
% dui=xlsread(filename,sheet,'B2:B11');
% sj=xlsread(filename,sheet,'C2:C21');

% b=[-dli;dui;-sj;sj];
% Ar=kron(eye(nShift),ones(1,nCombo));%Summation along rows of Xij can be expressed as a matrix-vector multiplication A1*XijT(:) 
% Ac=kron(ones(1,nShift),eye(nCombo));% Summation along columns can be expressed the same way with

  %Note:  for integer programming, Aeq and beq must be empty ([]), so
  %equality contraints are rewritten to si<=Ac<=sj
%A=[-Ar;Ar;-Ac;Ac];

% filename = 'Book1.xlsx';
% sheet = 1;
% xlRange = 'A1:A420';
% b = xlsread(filename,sheet,xlRange);

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

fViolation=cij.*Xij;
Tfvio=sum(fViolation(:));

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

TotalV=Tfvio+cons;
%ConstraintVio=sum(b-A*x');
% TotalV=Tfvio+(ConstraintVio);

% The double summation for Cij and Xij can be epxressed can be expressed C(:).'(X(:)).
end
