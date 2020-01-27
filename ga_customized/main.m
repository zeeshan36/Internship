tic
nShift=20;
nCombo=50;
t=5;
popsize=500;
nvars=nShift*nCombo;
prob_cross=0.9;
iteration=350;

filename = '2.xlsx';
sheet = 'A';
dli=xlsread(filename,sheet,'A2:A21');
dui=xlsread(filename,sheet,'B2:B21');
sj=xlsread(filename,sheet,'C2:C51');
b=[-dli;dui;-sj;sj];

x=1:iteration;
y=zeros(1,iteration);

p = 6.5;
q = 3.5;
r = 3.5;
parentChromosomeN=zeros(popsize,nvars+1);
parentChromosomeN(:,1:nvars) = IM(popsize,nShift,nCombo,dli,dui,sj);


%total_f = zeros(1,popsize);
for i=1:popsize
    [f,sum_sj,sum_dui,sum_dli] = Obj_f(parentChromosomeN(i,1:nvars),nShift,nCombo,dli,dui,sj,t);
    parentChromosomeN(i,nvars+1) = 1*f + p*(sum_sj)^1 + q*(sum_dui)^1 + r*(sum_dli)^1;
end
%parentChromosomeN(:,nvars+1)
itter=1;
while itter<=iteration
	matingPool= zeros(1,popsize);
	k=1;
    while k<=popsize
        randNum1=randi([1 popsize],1,1);
        randNum2=randi([1 popsize],1,1);
        while (randNum1 == randNum2)
            randNum2=randi([1 popsize],1,1);
        end
        if (parentChromosomeN(randNum1,nvars+1)<parentChromosomeN(randNum2,nvars+1))
            matingPool(k)=randNum1;
        else
            matingPool(k)=randNum2;
        end
        k=k+1;
    end

	childrenChromosomeN=zeros(popsize,nvars+1);
    m=1;
	while (m<round(prob_cross*popsize))
		chromosomeId1=matingPool(m);
		chromosomeId2=matingPool(m+1);
		parent1=parentChromosomeN(chromosomeId1,1:nvars);
		parent2=parentChromosomeN(chromosomeId2,1:nvars);
        [parent1, parent2] = mycrossover_all(parent1,parent2,nCombo);
		childrenChromosomeN(m,1:nvars)=parent1;
		childrenChromosomeN(m+1,1:nvars)=parent2;
		m=m+2;
	end
    m=1;
	while (m<=50)
		chromosomeId1=matingPool(m);
		parent1=parentChromosomeN(chromosomeId1,1:nvars);
        [parent1] = mymutation(parent1);
		childrenChromosomeN(m,1:nvars)=parent1;
		m=m+1;
	end


% 	x=zeros(1,2);                %* initializing 1x2 matrix with all zero;
% 	for k=1:50
% 		for i=1:8                %* binary to decimal conversion;
% 			x(1)= x(1)+childrenChromosomeN(k,9-i)*(2^(i-1));
% 			x(2)= x(2)+childrenChromosomeN(k,17-i)*(2^(i-1));
% 		end
% 		childrenChromosomeN(k,17)=shubert(x);
% 	end
    for i=1:popsize
        [f,sum_sj,sum_dui,sum_dli] = Obj_f(childrenChromosomeN(i,:),nShift,nCombo,dli,dui,sj,t);
        childrenChromosomeN(i,nvars+1) = 1*f + p*(sum_sj)^1 + q*(sum_dui)^1 + r*(sum_dli)^1;
    end
	%%---------/*Selecting elite chromosomes*/----------%%

	combinedMatrix=[parentChromosomeN;childrenChromosomeN];
    sortedMatrix = sortrows(combinedMatrix,nvars+1);
	parentChromosomeN=sortedMatrix(1:popsize,1:nvars+1);
    y(itter)=parentChromosomeN(1,nvars+1);
	itter=itter+1
end
toc
%time=toc

% x=zeros(1,2);             %* initializing 1x2 matrix with all zero;
% for i=1:8                 %;* binary to decimal conversion;
% 	x(1)= x(1)+parentChromosomeN(1,9-i)*(2^(i-1));
% 	x(2)= x(2)+parentChromosomeN(1,17-i)*(2^(i-1));
% end
 scatter(x,y);
 disp ("Final Solution");
 disp(parentChromosomeN(1,nvars+1));
 [f,sum_sj,sum_dui,sum_dli] = Obj_f(parentChromosomeN(1,1:nvars),nShift,nCombo,dli,dui,sj,t);
 %disp (x(1));
 %disp (x(2));
 