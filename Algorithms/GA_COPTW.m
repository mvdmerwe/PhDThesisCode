%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GA algorithm for solving cooperative orienteering problem.    %                                                               %
% Author: Martijn van der Merwe                                 %
%                                                               %
% Date: November 2014                                           %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [bestFitness, globalOptSol, pop,exeTime] = GA_COPTW(DataFile,TravelTimeFile,Nm)
global DATA N TravelTime Local_Search

coder.extrinsic('tic')

start = tic; %start timing
disp('GA Started')

%%% GA parameters %%%
Pop_Size = 100;     %
Cross_Prop = 0.7;   %
Mut_Prop = 0.4;     %
beta = 50;          %
Local_Search = 10;  %
initDensity = 0.4;  %
%%%%%%%%%%%%%%%%%%%%%

%Read data and calculate model parrameters.
TravelTime = csvread(TravelTimeFile);
DATA = csvread(DataFile);
N = int16(DATA(1,1));                           %Number of vertices
mu = Pop_Size * Cross_Prop;               %Number of offspring each generation
%protReq = DATA(2:N+1,8);


%% Initialise population
generation = 1;             %GA iteration counter
%globalOptVal(1) = 0;        %Stores the fitness value of the best solution found so far
globalOptSol = [];          %The best solution found
pop = zeros(Pop_Size,N*Nm); %Population of solutions

for i = 1:Pop_Size %Select routes at random
    
    vertexOrder = randperm(N-2)+1; %Determines the order in which assets will be visited    
    startPos = 1;
    for j = 1:Nm
        temp = vertexOrder;
        for k = 1:N-2 %Randomly select candidates for each vehicle path
            if rand < initDensity
                temp(k) = 0;
            end
        end
        indices = temp>0;
        content = temp(indices);
        endPos = startPos + length(content)+ 1;
        
        pop(i,startPos) = 1;
        pop(i,startPos+1:endPos-1) = content;
        pop(i,endPos) = N;
        startPos = endPos + 1;
    end
    
    if sum(pop(i,:) == N) < Nm    %-------------------------------
        disp('error!!!')          %An error check (can be removed)
        return                    %-------------------------------
    end                           %-------------------------------
end



%% Evaluate the population's fitness.
bestFitness(generation) = 0;
popFitness = zeros(1,Pop_Size);
popFeas = popFitness;
solution = zeros(1,N);
for i = 1:Pop_Size
    [popFitness(i), popFeas(i)] = evaluate(pop(i,:),Nm);
    if popFitness(i) > bestFitness(generation)
        bestFitness(generation) = popFitness(i);
        solution = pop(i,:);
    end
end
globalOptVal(generation) = bestFitness(generation);
globalOptSol  = solution;
disp(['Val: ', num2str(globalOptVal(generation))]);


%% Breed
disp('Breeding...')
unchanged = 0; %Counts the number of generation with no change in the best solution found.

%canPopFit = zeros(1,Pop_Size) - 1;
offspringFitness = zeros(1,mu) - 1;
while unchanged <= beta %&& toc < 20
    
    %Crossover
    offspringPop = zeros(mu,N*Nm);
    for i = 1:mu
        selected = randsample(Pop_Size,4);
        selFit = popFitness(selected);
        if selFit(1) > selFit(2)
            parent1 = pop(selected(1),:);
        else
            parent1 = pop(selected(2),:);
        end
        
        if selFit(3) > selFit(4)
            parent2 = pop(selected(3),:);
        else
            parent2 = pop(selected(4),:);
        end
        offspringPop(i,:)  = crossover(parent1,parent2,N,Nm);
    end 
    
    %Mutate
    selected = randsample(mu,mu*Mut_Prop); %Select individuals for mutation
    for i = selected
        %1-add, 2-remove, 3-replace
        [offspringPop(i,:),offspringFitness(i)] = mutate(offspringPop(i,:),randi(2),randi(Nm));
        
        %[test, ~] = evaluate(offspringPop(i,:),Nm);
        %if test ~= offspringFitness(i)
        %    display(['error, fitness miss match ', num2str(test),', ' ,num2str(offspringFitness(i))])            
        %end
    end
    
    mutateParentPop = true;
    if mutateParentPop
        selected2 = randsample(Pop_Size,Pop_Size*Mut_Prop); %Select individuals for mutation
        for i = selected2
            %1-add, 2-remove, 3-replace
            [pop(i,:), popFitness(i)] = mutate(pop(i,:),randi(2),randi(Nm));
        end
    end
    
    for i = 1:mu
        [offspringFitness(i), ~] = evaluate(offspringPop(i,:),Nm);
    end
    
        
    %Select individuals for next generation
    if Cross_Prop < 1
        newPop = zeros(Pop_Size,N*Nm);
        canPop = pop;
        newPopFit = zeros(1,Pop_Size);                
        canPopFit = popFitness;
        
        
        for i = 1:Pop_Size-mu
            selected = randsample(length(canPopFit),2);
            fit1 = canPopFit(selected(1));
            fit2 = canPopFit(selected(2));
            if fit1 > fit2
                win = selected(1);
            else
                win = selected(2);
            end
            
            newPop(i,:) = canPop(win,:);
            newPopFit(i,:) = canPopFit(win);
            canPop(win,:) = [];
            canPopFit(win) = [];
            
        end
        
        pop(1:mu,:) = offspringPop;
        pop(mu+1:Pop_Size,:) = newPop;
        popFitness(1:mu) = offspringFitness;
        popFitness(mu+1:Pop_Size) = newPopFit;
        
    else
        pop = offspringPop;
        popFitness = offspringFitness;
    end
    
    
    %Evaulate the new generations fitness
    generation = generation + 1;
    bestFitness(generation) = 0;
    for i = 1:Pop_Size
        %[popFitness(i), popFeas(i)] = evaluate(pop(i,:),Nm);
        if popFitness(i)  > bestFitness(generation)
            bestFitness(generation) = popFitness(i);
            solution = pop(i,:);
        end
    end
         
    if bestFitness(generation) > globalOptVal(generation-1)
        globalOptVal(generation) = bestFitness(generation);
        globalOptSol = solution;
        disp(['Val: ', num2str(globalOptVal(generation))]);
        unchanged = 0;
    else
        globalOptVal(generation) = globalOptVal(generation-1);
        unchanged = unchanged + 1;
    end
    
end

[val, feas] = evaluate(globalOptSol,Nm);
disp(['Solution:' ,num2str(val), ', ',num2str(feas) ]);
disp(['Generations: ',num2str(generation)]);
exeTime = toc(start)
end


%% Functions
function [offspring] = crossover(parent1, parent2,N,Nm)

cutPath = randi(Nm);                %Choose which path will be cut for the crossover
offspring = zeros(1,N*Nm);          %Initialise the offspring vector
pathEnds1 = find(parent1 == N);     %Find the indices of the endpoints of each path
pathEnds2 = find(parent2 == N);     %Find the indices of the endpoints of each path

maxP1 = pathEnds1(cutPath)-1;       %The maximum index value where the path 1 may be cut
maxP2 = pathEnds2(cutPath)-1;       %The maximum index value where the path 2 may be cut

if cutPath > 1
    minP1 = pathEnds1(cutPath-1)+1;
    minP2 = pathEnds2(cutPath-1)+1;
else
    minP1 = 1;
    minP2 = 1;
end

if maxP1 > minP1
    crossP1 = randi(maxP1-minP1+1) + minP1-1;
else
    crossP1 = minP1;
end

if maxP2 > minP2
    crossP2 = randi(maxP2-minP2+1) + minP2-1;
else
    crossP2 = minP2;
end

offspring(1:crossP1) = parent1(1:crossP1);
offspring(crossP1+1:pathEnds2(Nm)-(crossP2+1)+crossP1+1) = parent2(crossP2+1:pathEnds2(Nm));

pathStarts = find(offspring == 1);
pathEnds = find(offspring == N);

%Remove repeated entries
tempOffspring = offspring;
offspring = zeros(1,N*Nm);



uList = [];
if cutPath > 1;
    offspring(1:pathEnds(cutPath-1)) = tempOffspring(1:pathEnds(cutPath-1));
end
j  = pathStarts(cutPath);
for i = pathStarts(cutPath):pathEnds(cutPath)
    if tempOffspring(i) == 1 || tempOffspring(i) == N
        offspring(j) = tempOffspring(i);
        j = j + 1;
    else
        %check if unique
        if sum(uList == tempOffspring(i)) == 0;
            offspring(j) = tempOffspring(i);
            j = j + 1;
            uList = [uList,tempOffspring(i)];
        end
    end
end
if cutPath < Nm
    for i = pathStarts(cutPath+1):pathEnds(Nm)
        offspring(j) = tempOffspring(i);
        j = j+1;
    end
end
end


function [offspring, fitness] = mutate(parent,operation,pathN)
global N Local_Search

startP = find(parent == 1);
endP = find(parent == N);
Nm = length(endP);

pLength = sum(parent > 0);
subparent = parent(startP(pathN):endP(pathN));
spLength = sum(subparent > 0);

missing = setdiff(1:N,subparent);

if spLength <= 2
    operation = 1;
end

switch operation
    case 1 %%Add a point
        %missing = setdiff(1:N,parent);
        if ~isempty(missing)
            candidate = zeros(Local_Search,N*Nm);
            canFit = zeros(1,10);
            bestFit = -inf;
            bestFiti = 1;
            %pLength = sum(parent > 0);
            %i= 1;
            for i = 1:Local_Search
                subparent = parent(startP(pathN):endP(pathN));
                spLength = sum(subparent > 0);
                candidate(i,:) = parent;
                addPoint = missing(randi(length(missing)));
                if spLength-2 > 0
                    insPos = startP(pathN) + randi(spLength-2)-1;
                else
                    insPos = startP(pathN);
                end
                candidate(i,1:insPos) = parent(1:insPos);
                candidate(i,insPos+1) = addPoint;
                candidate(i,insPos+2:pLength+1) = parent(insPos+1:pLength);
                
                [canFit(i), ~ ] =  evaluate(candidate(i,:),Nm);
                if canFit(i) >= bestFit
                    bestFit = canFit(i);
                    bestFiti = i;
                end
            end
            
            offspring = candidate(bestFiti,:);
        else
            offspring = parent;
        end
        
    case 2 %%Omit a point if length is great than 2
        if spLength > 2
            candidate = zeros(10,N*Nm);
            canFit = zeros(1,10);
            bestFit = -inf;
            bestFiti = 1;
            for i = 1:10
                remPos = randi(spLength - 2)+ startP(pathN);
                candidate(i,1:pLength-1) = parent([1:remPos-1,remPos+1:pLength]);
                [canFit(i), ~ ] =  evaluate(candidate(i,:),Nm);
                if canFit(i) >= bestFit
                    bestFit = canFit(i);
                    bestFiti = i;
                end
            end
            offspring = candidate(bestFiti,:);
        else
            offspring = parent;
        end
        

        
end

fitness = bestFit;
end


function [fitness, feasible] = evaluate(PATH,Nm)

global N DATA %TravelTime %alpha %Tmax

%Read data
values = DATA(2:N+1,5);
O = DATA(2:N+1,6);
C = DATA(2:N+1,7);
protReq = DATA(2:N+1,8);


%Calculate service times
[S, ~] = BellmanFord(PATH);

for i = 1:length(S)
    S(i) = round(S(i) * 10) / 10;
end

%Calculate fitness
TimeLimitEx = 1;
i = 1;
pathEnds = find(PATH == N);
visited = zeros(1,N);
infeasible = 0;
while i <= pathEnds(Nm);
    if S(PATH(i)) > C(PATH(i))
        % fitness = fitness -1;
        infeasible = infeasible + 1;
        if PATH(i) == N
            %infeasible = infeasible + 100;
            TimeLimitEx = 0.9;
        end
    elseif S(PATH(i)) < O(PATH(i))
        disp('error')
    else
        visited(PATH(i)) =  visited(PATH(i))+1;
        
    end
    i = i + 1;
end

fitness = TimeLimitEx * sum(values(visited >= protReq')) -infeasible;

if S(N) > C(N)
    feasible = -1;
else
    feasible = 1;
end
    
    
end


function [weight, predecessor] = BellmanFord(solution)

global N TravelTime DATA

% This implementation takes in a graph, represented as
% lists of vertices and edges, and fills two arrays
% (weight and predecessor) with shortest-path
% (less cost/weight/metric) information

% Preprocessing
% graph -> list vertices, list edges, vertex source
vertices = unique(solution(solution>0));
source = 1;
ServiceTime = DATA(2:N+1,4);
O = DATA(2:N+1,6);


% Step 1: initialize graph
templ = zeros(1,length(vertices));
weight = templ;
predecessor = templ;
for v = vertices
    if v == source
        weight(v) = 0;
    else
        weight(v) = -1;
        predecessor(v) = NaN;
    end %if
end %for

% Step 2: relax edges repeatedly
for i = 1:20
    for j = 1:length(solution(solution > 0))
        if solution(j) ~= N
            u = solution(j);
            v = solution(j + 1);
            w = max(O(v), weight(u) + ServiceTime(u) + TravelTime(u,v));
            if  w > weight(v)
                weight(v) = w;
                predecessor(v) = u;
            end %if
        end
    end %for
end %for

% Step 3: check for negative-weight cycles
% for j = 1:length(solution(solution > 0)) %each edge (u, v) with weight w in edges:
%     if solution(j) ~= N
%         u = solution(j);
%         v = solution(j+1);
%        w = max(O(v), weight(u) + ServiceTime(u) + TravelTime(u,v));
%        if  w > weight(v)
%             disp('Graph contains a negative-weight cycle');
%         end
%     end
% end

end






