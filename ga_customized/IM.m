function [InitialMatrix] = IM(Intlpopsize,nShift,nCombo,dli,dui,sj)
%UNTITLED3 Summary of this function goes here
%   Give parameters same as GA_Test.m 

nvars=nCombo*nShift;
%disp("in IM");
InitialMatrix=zeros(Intlpopsize,nvars);
for g=1:Intlpopsize
    
%----------------------------------------------------------------------------------------
%Create an initial xij matrix, each column of which consists of a random binary xij
%sequence, but sum of ones in each column equals to sj of that combo
    Init_xij = zeros(nShift,nCombo);
    for j=1:nCombo
    randRows=randperm(nShift); %a row vector containing a random permutation of the integers from 1 to nShift inclusive.
    rowsWithOne=randRows(1:sj(j)); %row indices having 1, and sum up to sj
    Init_xij(rowsWithOne,j)=1;
    end
%----------------------------------------------------------------------------------------
%Adjust the initial xij matrix to make it feasible, satisfying horizontal
%LHS (dli) and RHS (dui) constraints
    Init_xij_Feasible=Init_xij;
    w=1;
    while w
        RowSum=sum(Init_xij_Feasible,2); %create a column vector containing the sum of each row
        CheckLB=lt(RowSum,dli); %if RowSum <dli, true
        CheckUB=gt(RowSum,dui); %if RowSum >dui, true
          if ~any(CheckLB)&&~any(CheckUB) %if any element in CheckLB and CheckUB is zero
            break,
          else
            [~,RowIdxMin]=min(RowSum);
            [~,RowIdxMax]=max(RowSum);
            ColIdx=find(Init_xij_Feasible(RowIdxMax,:) & ~Init_xij_Feasible(RowIdxMin,:),1); % returns the first 1 column index corresponding to the nonzero elements in row RowIdxMax and zero elements in row RowIdxMin.
            %swap the min and max elements
            [Init_xij_Feasible(RowIdxMin,ColIdx),Init_xij_Feasible(RowIdxMax,ColIdx)]=deal(Init_xij_Feasible(RowIdxMax,ColIdx),Init_xij_Feasible(RowIdxMin,ColIdx));
          end 
          
    w=w+1;
    end
%disp("Im out\n");
InitialMatrix(g,:)=reshape(Init_xij_Feasible,1,[]); % convert matrix to row vector
end