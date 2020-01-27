clc

%%
%%THIS CODE CONTAINS ALL THE TOPOLOGICAL MEASURES

% NodeStruct=shaperead('Nodes123')
% EdgeStruct=shaperead('Links12_f2l')

%X and Y positions of nodes
nX=([NodeStruct.X1]');
nY=([NodeStruct.Y1]');
%%
%Assign start and end node
for i=1:length(EdgeStruct)
    %Find start and end node based on integer/rounded position of ndoe and edge
    %locations; using find means that nEnd/nStart would be ID of the connected
    %node
    EdgeStruct(i).Id=i;
    % if multiple type nodes overlap eachother keep the community
    e= find((round((EdgeStruct(i).Xout),4))==round(nX,4) & (round((EdgeStruct(i).Yout),4))==round(nY,4));
    if length(e) >1
        for ei= 1:length(e)
            if NodeStruct(e(ei)).Type==3
                e = e(ei);
                break
            elseif NodeStruct(e(ei)).Type==1
                e = e(ei);
                break
            end
        end
    end
    EdgeStruct(i).EndNode= e;
    
    s=find((round((EdgeStruct(i).Xin),4))==round(nX,4) & (round((EdgeStruct(i).Yin),4))==round(nY,4));
    if length(s) >1
        for si= 1:length(s)
            if NodeStruct(s(si)).Type==3
                s = s(si);
                break
            elseif NodeStruct(s(si)).Type==1
                s = s(si);
                break
            end
        end
    end
    
    EdgeStruct(i).StartNode=s;
    
end

%% find demand and supply set
SN=0;

for i=1:length(NodeStruct)
    if (NodeStruct(i).Type== 3 & NodeStruct(i).CSDpop2016>=1000)
        SN=SN+1;
        ShelterNodes(1,SN)=i;
    end
end
CN=0;
for i=1:length(NodeStruct)
    if (NodeStruct(i).Type == 3 )
        CN=CN+1;
        CommunityNodes(1,CN)=i;
    end
end
%%
EdgeTable=table([EdgeStruct.StartNode;EdgeStruct.EndNode]',[EdgeStruct.length]',[EdgeStruct.Type]',[EdgeStruct.Id]','VariableNames',{'EndNodes','Weight','Type','ID'});
NodeTable=table([NodeStruct.X]',[NodeStruct.Y]',[NodeStruct.Type]',[NodeStruct.Id]','VariableNames',{'X1','Y1','Type','ID'});
%Actually creates a graph
G=graph(EdgeTable,NodeTable);

%plot(G,'XData',nX,'YData',nY);

%% shortest path and distance matrix
ToNode = 1:height(G.Nodes);
for i=1:height(G.Nodes)
    [TRR,DD] = shortestpathtree(G,i,ToNode,'OutputForm','cell');
    SDist(i,:) = DD;
    SPath(i,:)= TRR;
end

for i=1:height(G.Nodes)
    for j = 1:height(G.Nodes)
        if SDist(i,j)== Inf
            SDist(i,j)=999999;
        end
    end
end

%% pop nad betweenness

pop = [NodeStruct.CSDpop2016];
SDist_sn = [NodeStruct.Dis_To_Hwy];
fire = [NodeStruct.SUM_Fire];

%bet = [NodeStruct.Between];

%%
%sensitivity
% sen = zeros(20,10);
% for R = 1:20
%     for F=1:length(ShelterNodes)
%
F = 10;
RU = 150;
RL = 50;
% potential facility nodes, check if the facility is located within certain distance
DistCheck =(SDist<=RU & SDist> RL);

% subset of possible facility nodes for each community, find which nodes are within certain distance of each community, not used in further calculation
N = cell(height(G.Nodes),max(sum(DistCheck')));

for n =1:height(G.Nodes)
    index = find(DistCheck(n,:)==1);
    for m = 1:length(index)
        N{n,m} = index(m);
    end
end

%Community distance matrix
CDistCheck= DistCheck(CommunityNodes,ShelterNodes);

% Optimization
%
w1 = 1;                %   Population coverage weightage
w2 = 1;                %   Distance factor weightage
w3 = -1*1;                %   Fire coverage weitage
%
f     = -1.*[(w1.*(pop(CommunityNodes))+w3.*(fire(CommunityNodes)))'.*ones(length(CommunityNodes),1);...
    (w2.*SDist_sn(ShelterNodes)'.*ones(length(ShelterNodes),1))];
Aineq = [eye(length(CommunityNodes)),(-1*CDistCheck)] ;
bineq = zeros(length(CommunityNodes),1);
Aeq   = [zeros(length(CommunityNodes),1);ones(length(ShelterNodes),1)]';
beq   = F;

options = cplexoptimset;
options.Display = 'on';

[x, fval, eXintflag, output] = cplexbilp (f, Aineq, bineq, Aeq, beq, ...
    [ ], options);

fprintf ('\nSolution status = %s\n', output.cplexstatusstring);
fprintf ('Solution value = %d\n', fval);
disp ('Values = ');
disp (x');
disp ('Communities covered = ');
disp(CommunityNodes((x(1:length(CommunityNodes))'==1)));
disp ('# Communities covered = ');
disp (sum(((x(1:length(CommunityNodes))'==1))));
disp ('Communities not covered = ');
disp(CommunityNodes((x(1:length(CommunityNodes))'==0)));
disp ('#Communities not covered = ');
disp(sum((x(1:length(CommunityNodes))'==0)));
disp ('Shelter Nodes = ');
disp(ShelterNodes(x(end-(length(ShelterNodes))+1:end)'==1));
%
Cover=CommunityNodes((x(1:length(CommunityNodes))'==1));
NotCover= CommunityNodes((x(1:length(CommunityNodes))'==0));
Shelter= ShelterNodes(x(end-(length(ShelterNodes))+1:end)'==1);

CoverPop=0;
NotCoverPop=0;
for i=Cover
    CoverPop =CoverPop+pop(i);
end
for i=NotCover
    NotCoverPop =NotCoverPop+pop(i);
end
disp ('Covered Population = ');
disp(CoverPop);
disp ('Not Covered Population = ');
disp(NotCoverPop);
disp ('% Covered Population = ');
disp(CoverPop/(sum(pop)/2));

disp("------------%______________%---------------%______________%------------");
%disp([w1.*(pop(CommunityNodes).*ones(length(CommunityNodes),1));zeros(length(ShelterNodes),1)]*x);
%    sen(R,F)=CoverPop/(sum(pop)/2);
%
%    F=F+1;
%     end
%     R=R+1;
% end
