%clear all
start = tic;

%% TODO:
% - Does defining Z lead to improved computational time?
% - EdgeSet by scenario
% - Create scenario data



%% Problem Data
%Input files
fileName = 'H:\PhD\Programming\Stochastic\Stochastic-Programming\DemoData\StochasticDemonstrationData.csv';
VehiclePropFile = 'H:\PhD\Programming\Stochastic\Stochastic-Programming\DemoData\VehicleProperties.csv';
TravelTimeFile = 'H:\PhD\Programming\Stochastic\Stochastic-Programming\DemoData\TTFinal.csv';

%fileName = 'C:\Users\User\Documents\GitHub\Stochastic-Programming\DemoData\StochasticDemonstrationData.csv';
%VehiclePropFile = 'C:\Users\User\Documents\GitHub\Stochastic-Programming\DemoData\VehicleProperties.csv';
%TravelTimeFile = 'C:\Users\User\Documents\GitHub\Stochastic-Programming\DemoData\TTFinal.csv';

%numStageAssets = [8,15,15,15];
numStageAssets = [11,42,42,42];
scenProb = [1,0.33,0.33,0.34];
ScenarioSet = 1:4;
numScenarios = 4;

%% Comp testing
%fileName = 'test.csv';
%VehiclePropFile = 'VehicleProperties.csv';
%TravelTimeFile = 'TravelTime.csv';

%% Code starts 
DataFile = fileName;
display(DataFile);

%Read Data
DATA = csvread(DataFile);
TIME = csvread(TravelTimeFile);
VEHICLEPROP = csvread(VehiclePropFile);

%Parse data
N = DATA(1,1); %Total number of vertices in the graph representation of the problem
numDepots = DATA(1,2);
numStaging = DATA(1,3);
stageEndTime = DATA(1,4); %The time at which the second stage scenarios commence
Tmax = DATA(2,7)

service = DATA(2:N+1, 4);
value = DATA(2:N+1, 5);
openT = DATA(2:N+1, 6);
closeT = DATA(2:N+1, 7);

M = max(closeT) + max(max(TIME)) + max(service) - min(openT);


%Vehicle Data
numU = VEHICLEPROP(1,2);
TotalP = VEHICLEPROP(1,1);
numQ = VEHICLEPROP(1,3);
numP = VEHICLEPROP(1,4:3+numQ);
r(1:N,:) = DATA(2:N+1,8:7+numU);

%Define the sets
DepotSet = 1:numDepots;
counter = numDepots + 1;
StagingSet = counter:counter+numStaging-1;
counter = counter+numStaging;
for l = ScenarioSet
AssetSet{l} = counter:counter + numStageAssets(l)-1;
counter = counter + numStageAssets(l); 
end


index = 1;
for q = 1:numQ
    cap(q,:) = VEHICLEPROP(2:numU+1,index);
    index = index + numP(q);
end

counter = 1;
startVec = zeros(length(DepotSet),numQ);
for q = 1:numQ
    for i = 1:numP(q)
    k = VEHICLEPROP(numU+2,counter); 
    startVec(k,q) = startVec(k,q) + 1; 
    counter = counter+1;
    end
end

 
%EdgeSets for first stage 
for l = 1
    %Construct edge set
    k = 1;   
    newEdgeSet = [0 ,0];
    for i = [DepotSet, AssetSet{l}]
        for j = [AssetSet{l}, StagingSet]
            if i~=j && openT(i) + service(i)+ TIME(i,j) <= closeT(j) 
                newEdgeSet(k,:) = [i,j];
                k = k + 1;
            end
        end
    end
    EdgeSet{l} = newEdgeSet;
    %EdgeSetMatrix
    
    newEdgeSetMatrix = zeros(N,N);
    for i = [DepotSet, AssetSet{l}]
        for j = [AssetSet{l}, StagingSet]
            if i~=j && openT(i) + service(i)+ TIME(i,j) <= closeT(j) && j ~= 1 && ~ismember(j,DepotSet)
                newEdgeSetMatrix(i,j) = 1;
            end
        end
    end
    EdgeSetMatrix{l} = newEdgeSetMatrix;
end

%EdgeSets for second stage
for l = ScenarioSet(2:end)
    %Construct edge set
    k = 1;   
    newEdgeSet = [0 ,0];
    for i = [AssetSet{l}, StagingSet]
        for j = [AssetSet{l}, N]
            if i~=j && openT(i) + service(i)+ TIME(i,j) <= closeT(j) && j ~= 1 && ~ismember(j,DepotSet)
                newEdgeSet(k,:) = [i,j];
                k = k + 1;
            end
        end
    end
    EdgeSet{l} = newEdgeSet;
    %EdgeSetMatrix
    
    newEdgeSetMatrix = zeros(N,N);
    for i = [AssetSet{l}, StagingSet]
        for j = [AssetSet{l}, N]
            if i~=j && openT(i) + service(i)+ TIME(i,j) <= closeT(j) && j ~= 1 && ~ismember(j,DepotSet)
                newEdgeSetMatrix(i,j) = 1;
            end
        end
    end
    EdgeSetMatrix{l} = newEdgeSetMatrix;
end


clear newEdgeSetMatrix newEdgeSet

display('Building model...');
toc(start);
%% Problem variables
variableIndex = 1;

display('S...');
%S(i,l):real for all i in AssetSet(l), l in ScenarioSet;
S = zeros(N,numScenarios);
for l = ScenarioSet
    for i = AssetSet{l}
        S(i,l) = variableIndex;
        ctype(variableIndex) = 'C';
        colNames{variableIndex} = ['S[',num2str(i),',',num2str(l),']'];
        lb(variableIndex) = 0;
        if l == 1
            ub(variableIndex) = stageEndTime+30; 
        else
            ub(variableIndex) = Tmax;
        end
        variableIndex = variableIndex + 1;
    end
end

display('X...');
%X(EdgeSet(l),l):integer;
X = zeros(N,N,numQ,numScenarios);
for q = 1:numQ
    for l = ScenarioSet
        for k = 1:size(EdgeSet{l},1)
            i = EdgeSet{l}(k,1);
            j = EdgeSet{l}(k,2);
            X(i,j,q,l) = variableIndex;
            ctype(variableIndex) = 'I';
            colNames{variableIndex} = ['X[',num2str(i),',',num2str(j),',' ,num2str(q),',',num2str(l),']'];
            lb(variableIndex) = 0;        
            ub(variableIndex) = numP(q);
            variableIndex = variableIndex + 1;
        end
    end
end

display('Z...');
%Z(EdgeSet(l),l):binary;
Z = zeros(N,N,numQ,numScenarios);
for q = 1:numQ
    for l = ScenarioSet
        for k = 1:size(EdgeSet{l},1)
            i = EdgeSet{l}(k,1);
            j = EdgeSet{l}(k,2);
            Z(i,j,q,l) = variableIndex;
            ctype(variableIndex) = 'B';
            colNames{variableIndex} = ['Z[',num2str(i),',',num2str(j),',',num2str(q),',',num2str(l),']'];            
            lb(variableIndex) = 0;        
            ub(variableIndex) = 1;
            variableIndex = variableIndex + 1;
        end
    end
end

display('Y...');
%Y(AssetSet(l),l):integer;
Y = zeros(N,numScenarios);
for l = ScenarioSet
    for i = AssetSet{l}
        Y(i,l) = variableIndex;
        ctype(variableIndex) = 'B';
        colNames{variableIndex} = ['Y[',num2str(i),',',num2str(l),']'];        
        lb(variableIndex) = 0;        
        ub(variableIndex) = 1;
        variableIndex = variableIndex + 1;
    end
end

numVariables = variableIndex - 1;

display('Objective...');
toc(start)
%% Objectives
f = zeros(1,numVariables);
for l = ScenarioSet
    for i = AssetSet{l}
        f(Y(i,l)) = scenProb(l) * value(i);
    end
end
%f = -f; %Maximise the objective


%% Constraints
Aeq = []; % zeros(1,numVariables); %equality constraints.
beq = [];
eqi = 1;
Aineq = []; %zeros(1,numVariables); %equality constraints.
bineq = [];
ineqi = 1;

display('Start locations...');
toc(start)
%% Start locations: sum X^0_kjq <= start_kq forall k in DepotSet and q in TypeSet
for q = 1:numQ
    for k = DepotSet
        row = zeros(1,numVariables);
        for j = find(EdgeSetMatrix{1}(k,:)) 
            row(X(k,j,q,1)) = 1;
        end
        Aineq(ineqi,:) = row;
        bineq(ineqi) = startVec(k,q);
        ineqi = ineqi+1;
    end
end
 
display('Staging constriants...');
toc(start)
%% Staging constraints
% \sum_{k \in \mathcal{A}_d} \sum_{ j \in \Delta^+_q(k) }  X^0_{kjq} = \sum_{k \in \mathcal{A}_s}\sum_{ i \in \Delta^-_q(k) }  X^0_{ikq} \hspace{12pt} \forall\ q \in \mathcal{Q};
for q = 1:numQ
    row = zeros(1,numVariables);
    for k = DepotSet
        for j = find(EdgeSetMatrix{1}(k,:))
            row(X(k,j,q,1)) = 1;
        end
    end
    
    for k = StagingSet
        for i = find(EdgeSetMatrix{1}(:,k))
            row(X(i,k,q,1)) = -1;
        end
    end
       
    Aeq(eqi,:) = row;
    beq(eqi) = 0;
    eqi = eqi+1;
end

display('Flow balancing...');
toc(start)
%% Balance flow
%Aeq = [Aeq; zeros(length(AssetSet),numVariables)];
%beq = [beq, zeros(1,length(AssetSet))];
for q = 1:numQ
    for l = ScenarioSet
        for k = AssetSet{l}
            row = zeros(1,numVariables);
            for i = find(EdgeSetMatrix{l}(:,k))                
                    row(X(i,k,q,l)) = 1;                
            end
            
            for j = find(EdgeSetMatrix{l}(k,:))                
                    row(X(k,j,q,l)) = -1;
            end            
            Aeq(eqi,:) = row;
            beq(eqi) = 0;
            eqi = eqi+1;
        end
    end
end

display('Protection requirements...');
toc(start)
%% Protection requirements
% r_ku Y_k  - sum sum cap_qu X_ikq <= 0 for each u, k, and l
for u = 1:numU
    for l = ScenarioSet
        for k = AssetSet{l}        
            row = zeros(1,numVariables);                        
            row(Y(k,l)) = r(k,u);
            for q = 1:numQ
                for i = find(EdgeSetMatrix{l}(:,k))
                    row(X(i,k,q,l)) = -cap(q,u);
                end
            end            
            Aineq(ineqi,:) = row;
            bineq(ineqi) = 0;
            ineqi = ineqi+1;
        end
    end
end

display('Linking contstraints...');
toc(start)
%% Linking constraints 
% sum_i X(i,k,q,0) - sum_j X(k,j,q,l) for all l in ScenarioSet and q in Q,
% k in StagingSet

for k = StagingSet
    for q = 1:numQ
        for l = ScenarioSet(2:end)
            row = zeros(1,numVariables);

            for i = [find(EdgeSetMatrix{1}(:,k))]'
                row(X(i,k,q,1)) = 1;
            end


            for j = find(EdgeSetMatrix{l}(k,:))
                row(X(k,j,q,l)) = -1;
            end

            Aeq(eqi,:) = row;
            beq(eqi) = 0;
            eqi = eqi+1;        
        end
    end
end

display('Sink...');
toc(start)
%% Sink constraint
for q = 1:numQ
    for l = ScenarioSet(2:end)
        row = zeros(1,numVariables);
        
        for k = StagingSet
            for j = find(EdgeSetMatrix{l}(k,:))
                row(X(k,j,q,l)) = 1;
            end
        end
        
        for i = find(EdgeSetMatrix{l}(:,N))
            row(X(i,N,q,l)) = -1;
        end
        
        Aeq(eqi,:) = row;
        beq(eqi) = 0;
        eqi = eqi+1;
    end
end

display('Timing...');
toc(start)
%% Timing constraints 
% X_ijql - p_q Z_ijql <= 0 for all q in Q, l in ScenarioSet 
aVariable = 0;
for l = ScenarioSet
    aVariable = aVariable + size(EdgeSet{l},1);
end

aVariable = numQ * aVariable;

%Aineq = [Aineq; zeros(aVariable,numVariables)];
bineq = [bineq, zeros(1,aVariable)];
Aineq = sparse(Aineq);


for q = 1:numQ
    for l = ScenarioSet
        for k = 1:size(EdgeSet{l},1)
            row = zeros(1,numVariables);
            i = EdgeSet{l}(k,1);
            j = EdgeSet{l}(k,2);
            row(X(i,j,q,l)) = 1;
            row(Z(i,j,q,l)) = - numP(q);
            Aineq(ineqi,:) = row;
            bineq(ineqi) = 0;
            ineqi = ineqi + 1;
        end        
    end
end

%Si - Sj + M*Zijp  <= M -  tijp - ai
for q = 1:numQ
    for l = ScenarioSet
        Aineq = [Aineq; zeros(size(EdgeSet{l},1),numVariables)];
        bineq = [bineq, zeros(1,size(EdgeSet{l},1))];
        for k = 1:size(EdgeSet{l},1)
            endVertexCoeff = 0;
            stagingStart = 0;
            row = zeros(1,numVariables);
            i = EdgeSet{l}(k,1);
            j = EdgeSet{l}(k,2);
            
            if ismember(i,DepotSet) %The start time is 0
                %Leave this blank -> S_i = 0;
            elseif ismember(i,StagingSet) %The start time is t, must be subtracted from RHS
                stagingStart = stageEndTime;
            elseif ismember(i, AssetSet{l})
                row(S(i,l)) = 1;
            else
                error('Something went wrong in the timing constraints.')
            end
            
            if j == N
                endVertexCoeff = Tmax;
            elseif ismember(j,StagingSet)
                endVertexCoeff = stageEndTime;
            else
                row(S(j,l)) = -1;
            end
            row(Z(i,j,q,l)) = M;
            
            Aineq(ineqi,:) = row;
            bineq(ineqi) = M - TIME(i,j)-service(i) +  endVertexCoeff - stagingStart;
            ineqi = ineqi + 1;
        end
    end
end

display('Time windows...');
toc(start)
%% Time windows
% S_i <= c_i
% -S_i <=  -o_i

%Aineq = [Aineq; zeros(2*N-4 ,numVariables)];
%bineq = [bineq, zeros(1,2*N-4)];
for l = ScenarioSet
    for i = AssetSet{l}
        row = zeros(1,numVariables);
        row(S(i,l)) = -1;
        Aineq(ineqi,:) = row;
        bineq(ineqi) = -openT(i);
        ineqi = ineqi + 1;
        
        row = zeros(1,numVariables);
        row(S(i,l)) = 1;
        Aineq(ineqi,:) = row;
        bineq(ineqi) = closeT(i);
        ineqi = ineqi + 1;
    end
end

%% Bounds on variables


%% Solve using cplex
% options = cplexoptimset('cplex');
% options.Display = 'on';
% options.workmem = 2^16;
% options.workdir = '/lustre/pRMIT0107/WorkDirTOPTW';
% options.mip.strategy.variableselect = 3; 
% [solution, fval,~, output] = cplexmilp(f,Aineq,bineq',Aeq,beq',  [], [], [ ], lb, ub, ctype, [ ], options);
% disp(fval);
 
display('Model build done, defining cplex object...');
cplex = Cplex('');
cplex.Model.sense = 'maximize';  
cplex.Model.obj   = f';
cplex.Model.ctype = ctype;
cplex.Model.A =  [Aeq ; Aineq ];
cplex.Model.lhs = [ beq'; -inf( length(bineq),1)];
cplex.Model.rhs = [ beq';bineq' ];
      
cplex.Model.lb = lb; %[0*ones(size(f))];
cplex.Model.ub = ub; %[inf*ones(size(f))];
cplex.Model.colname = colNames;
cplex.Param.timelimit.Cur = 600;
elapsedTime = toc(start);

display('Solving...');
output = cplex.solve;
totalTime = toc(start);




