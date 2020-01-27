%MPE resource scheduling model (MPE-RSM)solving by Dantzig-Wolfe
%Decompostiion and Column Generation (also referred to as Branch and Price)

tic
%----------------------------------------------------------------------------------------
%Notations
%cij: cost of shift i on site combination (combo) j
%xij:association between shift i and combo j
%dli and dui:lower and upper bounds on available visits during shift i
%sj:total visits on combo j
%t:number of consecutive shifts observing time halo effect


%----------------------------------------------------------------------------------------
%Step I - Restricted Master Problem (RMP) considering 1 initial known shift assignment (column) for
%each combo j

%Parameter Setting
% nShift=4;
% nCombo=5;
nShift=60;
nCombo=150;
t=7;
Intlpopsize=100;
nvars=nShift*nCombo;
% sj=[1;2;3;3;2];
% dli=ones(4,1);
% dui=5*ones(4,1);
filename = 'Book1.xlsx';
sheet = 1;
xlRange = 'A1:A420';
% b = xlsread(filename,sheet,xlRange);

dli=xlsread(filename,sheet,'F2:F61');
dui=xlsread(filename,sheet,'G2:G61');
sj=xlsread(filename,sheet,'H2:H151');
b1=[-dli;dui;-sj;sj];

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
%{
%----------------------------------------------------------------------------------------
%Compute cij for the initial xij_Feasible matrix
Init_cij = zeros(nShift,nCombo);
Init_xij_Feasible_forCost=[Init_xij_Feasible;Init_xij_Feasible(1:t-1,:)];%compute the remaining t-1 element's cost
for j=1:nCombo
    for i=1:nShift
        Init_cij(i,j)=sum(Init_xij_Feasible_forCost((i:i+t-1),j));
    end
end

%----------------------------------------------------------------------------------------
% Build cost coeffiicent in the objective function of the RMP
Init_cij_T=Init_cij'; %transpose cij
for j=1:nCombo
    ObjC_InitRMP(j,1)=Init_cij_T(j,:)*Init_xij_Feasible(:,j);
end

%----------------------------------------------------------------------------------------
% Build constraint matrix A, LHS, and RHS
D=diag(ones(1,nCombo)); % create a 145*145 square diagonal matrix D with ones on the main diagonal.
A=[Init_xij_Feasible;-Init_xij_Feasible;D];
LHS=[dli;-dui;ones(nCombo,1)];
RHS=inf([2*nShift+nCombo,1]);

%----------------------------------------------------------------------------------------
% Implement Cplex Engine to solve RMP
%1)Initialize the CPLEX object
RMPOpt=Cplex();
RMPOpt.Model.sense='minimize';
% 2)Fill in the data for the problem by populate by column
RMPOpt.addRows(LHS, [], RHS);
RMPOpt.addCols(ObjC_InitRMP, A);
% 3)Optimize the problem
RMPOpt.solve();
%4)report(RMPOpt);
Init_RMP_Solution_objval=RMPOpt.Solution.objval;
% disp ('Initial RMP Objective Value=');  
% disp (RMPOpt.Solution.objval);
% disp ('optimal value of variables in RMP=');
% disp (RMPOpt.Solution.x'); 
% disp ('Dual=');
% disp (RMPOpt.Solution.dual);
RMP_dual=RMPOpt.Solution.dual;


%----------------------------------------------------------------------------------------
%Step II - Restricted Master Problem (RMP) considering 1 initial known shift assignment (column) for
%each combo j

%----------------------------------------------------------------------------------------
%Compute Hessian Matrix Q of the objective function of the pricing subproblem:(1/2)X'QX for f=sum (i in nShift-t+1)(x[i]+x[i+1]+x[i+2])*x[i]
x = sym('x',[1 nShift+t-1]); % vector of variables
f = 0;
for i = 1:nShift
    f = f+x(i)*sum(x(i:i+t-1));
end

for i=1:t-1
    f=subs(f,x(nShift+i),x(i));
end
symhessian=jacobian(gradient(f,x),x);
numhessian=double(symhessian);% convert symbolic matrix into numeric one 
numhessian=numhessian(1:nShift,1:nShift);
%----------------------------------------------------------------------------------------
%Solve Pricing subproblems
k=1;
while k
    %get optimal dual multipliers of the RMP
    RMPOpt.solve();
    RMP_x=RMPOpt.Solution.x'; 
    RMP_dual=RMPOpt.Solution.dual;
    
    %q'X component in the objective function of the pricing subproblems
    q=RMP_dual(nShift+1:2*nShift)-RMP_dual(1:nShift);
    
    
    %Solve each Pricing subproblem iteratively
    for j=1:nCombo
        % Initialize the CPLEX object:(1/2)X'QX+q'X
        SubOpt=Cplex();
        SubOpt.Model.sense='minimize';
        % Fill in the data for the problem with populatebyrow
        SubOpt.addCols(q, [], zeros(nShift,1), ones(nShift,1)); 
        SubOpt.Model.Q = numhessian;
        SubOpt.addRows(sj(j), ones(1,nShift), sj(j));
        SubOpt.Model.ctype = char (ones ([1, length(SubOpt.Model.obj)]) * ('I'));
        
        % Set cplex parameter to solve non-convex QP problem
        SubOpt.Param.optimalitytarget.Cur=3; % When this parameter is set to 3, if the problem type is QP, CPLEX first changes the problem type to MIQP. CPLEX then solves the problem to global optimality. In problems of type QP or MIQP, this setting interacts with linearization switch for QP, MIQP.
        
        % Set cplex parameter to improve MIP solving speed (parameters are
        % decided by implementing Cplex's automatic tunning tool)
        %SubOpt.Param.mip.limits.eachcutlimit.Cur=0; %Sets a limit for each type of cut; A setting of 0 means no cuts.
        %SubOpt.Param.mip.limits.cutsfactor.Cur=1; %A CutsFactor of 1.0 or less means that no cuts will be generated.
        %SubOpt.Param.mip.strategy.heuristicfreq.Cur=-1; %heuristic off
        SubOpt.Param.mip.cuts.flowcovers.Cur=1;
        SubOpt.Param.mip.cuts.mircut.Cur=1;
        SubOpt.Param.mip.strategy.backtrack.Cur=0.1;
        SubOpt.Param.mip.strategy.heuristicfreq.Cur=100;
        SubOpt.Param.mip.strategy.presolvenode.Cur=2;
        SubOpt.Param.mip.strategy.probe.Cur=3;
        SubOpt.Param.mip.strategy.variableselect.Cur=4;
        SubOpt.Param.mip.tolerances.mipgap.Cur=0.05;
        % Optimize the problem
        SubOpt.solve();
        
        % Write the solution
       fprintf ('\nSolution status = %s\n', SubOpt.Solution.statusstring);
       fprintf ('Solution value = %f\n', SubOpt.Solution.objval);
       disp ('Values = ');
       disp (SubOpt.Solution.x');
%      disp ('Slacks = ');
%      disp (SubOpt.Model.rhs - SubOpt.Solution.ax);
       
       %save the new solution to subproblem j
       Sub_x(:,j)=SubOpt.Solution.x;
       Sub_objval(j)=SubOpt.Solution.objval;
    end
    
    %----------------------------------------------------------------------------------------
    %Check the reduced cost for all subproblems' new identified cols, and
    %find the one with the minimum reduced cost
    ReducedCost=Sub_objval'-RMP_dual(2*nShift+1:end);
    %report
    disp ('Reduced costs of Subproblems');
    disp (ReducedCost);
    
    % if all columns' reduced costs are nonnegative, stop
    if all(ReducedCost>=-1.0e-4)==1 %1 is ture
        break,
    end
    
    %----------------------------------------------------------------------------------------
    %record the newly identified column with the most negative reduced cost
    [M,I]=min(ReducedCost);
    B=zeros(nCombo,1);
    B(I)=1;
    newPatt=Sub_x(:,I);
    newPatt_forCost=[newPatt;newPatt(1:t-1,:)];
    
    Sub_ci=zeros(1,nShift);
    for i=1:nShift
        Sub_ci(i)=sum(newPatt_forCost(i:i+t-1));
    end
    
    %report
    disp ('UseSubProblemNo=');
    disp (I);
    disp ('Reduced cost=');
    disp (M);
    disp ('Use=');
    disp (newPatt);
    
    %----------------------------------------------------------------------------------------
    %add columns to the RMP
    ObjC_AddedCol=Sub_ci*newPatt;
    RMPOpt.addCols(ObjC_AddedCol, [newPatt;-newPatt;B]);
    
    % Write the solution
    Rprt_AddedCols(:,k)=[I;M;ObjC_AddedCol;newPatt;-newPatt;B];
    
    k = k+1;
end

toc
%-----------------------------------------------------------------------------
% %Reporting
Rprt_RMP_LP_Time=toc;
% 
% Rprt_RMP_LP_Solution_objval=RMPOpt.Solution.objval;
% Rprt_RMP_LP_Solution=RMPOpt.Solution.x';
% Rprt_SchedulePool_SN=[[(1:nCombo);Init_xij_Feasible]  [Rprt_AddedCols(1,:);Rprt_AddedCols(4:3+nShift,:)]];
% %Rprt_SchedulePool([1],:)=[];
% 
% %Compute cij for the final xij pool matrix
% [sm,sn]=size(Rprt_SchedulePool);
% Pool_cij = zeros(sm,sn);
% Rprt_SchedulePool_forCost=[Rprt_SchedulePool;Rprt_SchedulePool(1:t-1,:)];%comp
% for j=1:sn
%     for i=1:sm
%         Pool_cij(i,j)=sum(Rprt_SchedulePool_forCost(i:i+t-1,j));
%     end
% end
% 
% % Build cost coeffiicent in the objective function of the RMP
% Pool_cij_T=Pool_cij'; %transpose cij
% for j=1:sn
%     ObjC_PoolRMP(j,1)=Pool_cij_T(j,:)*Rprt_SchedulePool(:,j);
% end
% CheckObjval_LP_RMP=Rprt_RMP_LP_Solution*ObjC_PoolRMP;
% 
% Rprt_SchedulePool_cost_SN=[ObjC_PoolRMP';Rprt_SchedulePool_SN];


%----------------------------------------------------------------------------------------
%Solve RMP again with all the identified newly entering column pool and
%find integer solutions
tic
RMPOpt.Model.ctype = char (ones ([1, length(RMPOpt.Model.obj)]) * ('I'));

%Set cplex parameter to improve MIP solving speed (parameters are decided by implementing Cplex's automatic tunning tool)
% RMPOpt.Param.mip.limits.cutsfactor.Cur=10; %less cut
RMPOpt.Param.mip.strategy.heuristicfreq.Cur=-1; %heuristic off

% Optimize the problem
RMPOpt.solve();

%report
toc
Rprt_RMP_IP_Time=toc;
Time=Rprt_RMP_IP_Time+Rprt_RMP_LP_Time;
Rprt_RMP_IP_Solution_objval=RMPOpt.Solution.objval;
Rprt_RMP_IP_Solution=RMPOpt.Solution.x';

fprintf ('Best integer Solution in MP %f\n', RMPOpt.Solution.objval);
% disp ('Best Schedule=');
% disp (RMPOpt.Solution.x');

% %write the optimal solution
% % Rprt_BestSchedule=RMPOpt.Solution.x';
% % Rprt_SchedulePool=[[(1:nCombo);Init_xij_Feasible]  [Rprt_AddedCols(1,:);Rprt_AddedCols(4:3+nShift,:)]];
% Idx=find(Rprt_RMP_IP_Solution);
% Rprt_BestScheduleTable=Rprt_SchedulePool_cost_SN(:,Idx);
% Rprt_BestScheduleTable = sortrows(Rprt_BestScheduleTable',2)'; % sort according to second row which is comb # (first row is cost of each column)
% CheckObjval_IP_RMP=Rprt_RMP_IP_Solution*ObjC_PoolRMP;
% %check feasibility
% fprintf ('Feasible horizontally if returns 1 (true) %f\n', all(sum(Rprt_BestScheduleTable(2:end,:),2)>=dli & sum(Rprt_BestScheduleTable(2:end,:),2)<=dui));
% fprintf ('Feasible vertically if returns 1 (true) %f\n', isequal(sum(Rprt_BestScheduleTable(2:end,:),1)',sj));

%----------------------------------------------------------------------------------------
%}
