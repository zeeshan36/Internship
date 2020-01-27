%Parameter Setting
%Create GAP Set
%disp("started");
clear all
% global iteration
% global halo_t
% global fval
% global timeval
global nShift
global nCombo
global t
global dli
global dui
global sj
% prob_set=[75 200];
% shift=[5 10];
% combo=[15 20];
prob_set=[75 200 1000 4000 9000];
shift=[5 10 20 40 60];
combo=[15 20 50 100 150];
iteration=10;
halo_t=[3 5 7 10];
sol=cell(length(prob_set),length(halo_t));
Intlpopsize=500;
run_time=[500 500 1000 10000 10000];
% fval=zeros(length(halo_t),iteration);
% timeval=zeros(length(halo_t),iteration);

for prob=1:length(prob_set)
    A=[];
    MM=[];
    if prob_set(prob)==75
        halo_t=[3 5];
    else
        halo_t=[3 5 7 10];
    end
    for gap=1:length(halo_t)
        for itter=1:iteration
            tic
            nShift=shift(prob);
            nCombo=combo(prob);
            t=halo_t(gap);
            nvars=nShift*nCombo;

            filename = 'DATA.xlsx';
            sheet = int2str(prob_set(prob));
            range_dl=strcat('A2:A',int2str(nShift+1));
            range_du=strcat('B2:B',int2str(nShift+1));
            range_sj=strcat('C2:C',int2str(nCombo+1));
            dli=xlsread(filename,sheet,range_dl);
            dui=xlsread(filename,sheet,range_du);
            sj=xlsread(filename,sheet,range_sj);
            
            b=[-dli;dui;-sj;sj];
            lb=zeros(1,nvars);
            ub=ones(1,nvars);
            IntCon=[1:nvars];
            %build constraint matrix
              %ga evaluates the matrix product A*x as if x is transposed (A*x').
            Ar=kron(eye(nShift),ones(1,nCombo));%Summation along rows of Xij can be expressed as a matrix-vector multiplication A1*XijT(:) 
            Ac=kron(ones(1,nShift),eye(nCombo));% Summation along columns can be expressed the same way with

          %Note:  for integer programming, Aeq and beq must be empty ([]), so
          %equality contraints are rewritten to si<=Ac<=sj
        A=[-Ar;Ar;-Ac;Ac];
        nonlcon = [];
        MM=IM(Intlpopsize,nShift,nCombo,dli,dui,sj);
        %PopulationSize: Positive integer | {50} whennumberOfVariables <= 5, {200} otherwise |{min(max(10*nvars,40),100)} for mixed-integer problems

        opts = optimoptions(@ga, ...
                            'PopulationSize', Intlpopsize, ...
                            'Generations', 2000, ...
                            'MaxStallGenerations',200,...%1e10%
                            'FunctionTolerance',1e-8, ...
                            'TolCon', 1e-8,...%
                            'ConstraintTolerance',1e-8, ...
                            'InitialPopulationMatrix',MM,...%
                            'MaxTime',run_time(prob),...
                            'PlotFcns', {@gaplotbestf,@gaplotstopping,@gaplotbestindiv});

        [x,faval,exitflag]=ga(@Obj_f,nvars,A,b,[],[],lb,ub,nonlcon,IntCon,opts);
        toc
        Time=toc;
        [TotalV,Tfvio,cons] = Obj_f(x);
        sol{prob,gap}(1,itter)=Tfvio;
        sol{prob,gap}(2,itter)=cons;
        sol{prob,gap}(3,itter)=Time;
       %sol{4}(gap,itter)=
        %average_time=average(sol{3})
        end
    end
end