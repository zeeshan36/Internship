/*********************************************
 * OPL 12.8.0.0 Model
 * Author: li18
 * Creation Date: 2019年4月25日 at 下午3:37:34
 *********************************************/
int nc=...; // # of combinations
int ns=...; // # of shifts 

float temp;
execute{
var before = new Date();
temp = before.getTime();
}

range comb=1..nc; 
range shift=1..ns; 
//define decision variables
dvar boolean x[shift][comb];

//define parameters 
int CombSF[comb]=...;
//dexpr int k[i in shift]=1+(i-1)%ns; k=t-1
dexpr int c[i in shift, j in comb]= sum (k in 0..6)x[1+(i+k-1)%ns][j];
int dl[shift]=...;
int du[shift]=...;

//define objective function
minimize sum(i in shift, j in comb)c[i][j]*x[i][j];

//define constraintss
subject to 
{ 
Con01:
forall(i in shift)
  dl[i]<=sum(j in comb)x[i][j]<=du[i];

Con02:
forall(j in comb)
  sum(i in shift)x[i][j]==CombSF[j];
}

// solve the model
execute{
var after = new Date();
writeln("solving time ~= ",after.getTime()-temp);
}