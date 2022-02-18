function [bestSoFarFit ,bestSoFarSolution ...
    ]=simpleEA( ...  % name of your simple EA function
    fitFunc, ... % name of objective/fitness function
    T, ... % total number of evaluations
    input) % replace it by your input arguments

% Check the inputs
if isempty(fitFunc)
  warning(['Objective function not specified, ''' objFunc ''' used']);
  fitFunc = 'objFunc';
end
if ~ischar(fitFunc)
  error('Argument FITFUNC must be a string');
end
if isempty(T)
  warning(['Budget not specified. 1000000 used']);
  T = '1000000';
end
eval(sprintf('objective=@%s;',fitFunc));
% Initialise variables
nbGen = 0; % generation counter
nbEval = 0; % evaluation counter
bestSoFarFit = 0; % best-so-far fitness value
bestSoFarSolution = NaN; % best-so-far solution
%recorders
fitness_gen=[]; % record the best fitness so far
solution_gen=[];% record the best phenotype of each generation
fitness_pop=[];% record the best fitness in current population 
%% Below starting your code
nbPop = 4;
lb = 0;
ub = 31;
lenGene = 5;
% Initialise a population
%% TODO
pop = randi([lb,ub],nbPop,1);
genotype = dec2bin(pop,lenGene);

% Evaluate the initial population
%% TODO
fitness = objective(pop);
[m,idx]=max(fitness);
if m > bestSoFarFit
    bestSoFarFit = m;
    bestSoFarSolution = pop(idx);
end
fitness_gen = [fitness_gen,bestSoFarFit];
solution_gen = [solution_gen,bestSoFarSolution];
fitness_pop = [fitness_pop,m];

nbEval = nbEval + nbPop;
nbGen =  nbGen +1;
% Start the loop
while (nbEval<T)
% Reproduction (selection, crossver)
%% TODO
selectProb = fitness./sum(fitness);
newpop_genotype = [];
for i=1:nbPop/2
    r1 = rand();
    r2 = rand();
    p1_idx = 0;
    p2_idx = 0;
    s = 0;
    for index = 1:nbPop
        s = s+selectProb(index);
        if p1_idx==0 && r1<s
            p1_idx = index;
        end
        if p2_idx==0 && r2<s
            p2_idx = index;
        end
    end
    crosspoint = randi(lenGene-1);
    c1 = [genotype(p1_idx,1:crosspoint),genotype(p2_idx,crosspoint+1:end)];
    c2 = [genotype(p2_idx,1:crosspoint),genotype(p1_idx,crosspoint+1:end)];
    newpop_genotype = [newpop_genotype;c1;c2];
end

% Mutation
%% TODO
mutationProb = 1/lenGene;
for i=1:nbPop
    for j=1:lenGene
        if rand()<mutationProb
            if newpop_genotype(i,j) == '0'
                newpop_genotype(i,j) = '1';
            else
                newpop_genotype(i,j) = '0';
            end
        end
    end
end
newpop = bin2dec(newpop_genotype);

pop = newpop;
genotype = newpop_genotype;

fitness = objective(pop);
[m,idx]=max(fitness);
if m > bestSoFarFit
    bestSoFarFit = m;
    bestSoFarSolution = pop(idx);
end
fitness_gen = [fitness_gen,bestSoFarFit];
solution_gen = [solution_gen,bestSoFarSolution];
fitness_pop = [fitness_pop,m];

nbEval = nbEval + nbPop;
nbGen =  nbGen +1;

end

bestSoFarFit
bestSoFarSolution

figure,plot(1:nbGen,fitness_gen,'b') 
title('Fitness\_Gen')

figure,plot(1:nbGen,solution_gen,'b') 
title('Solution\_Gen')

figure,plot(1:nbGen,fitness_pop,'b') 
title('Fitness\_Pop')
