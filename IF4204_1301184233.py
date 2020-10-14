import numpy as np
import random

# Desain Kromosom dan Metode Pendekodean

## Binary Encoding using 3 bits (3 Gens)
def createPopulation(jumlahKromosom,jumlahIndividu):
	pop = []
	i = 0
	while i < jumlahIndividu:
		gen = []
		j = 0
		while j<jumlahKromosom:
			gen.append(random.randint(0,1))
			j += 1
		pop.append(gen)
		i += 1
	return pop

## Decode Chromosome
def decodeChromosome(rminx1, rmax1, rminx2, rmax2,individu):
	x1 = 	rminx1 + (rmax1 - rminx1)/(pow(2,-1) + pow(2,-2) + pow(2,-3)) * (individu[0]*pow(2,-1) + individu[1]*pow(2,-2) + individu[2]*pow(2,-3))
	x2 = 	rminx2 + (rmax2 - rminx2)/(pow(2,-1) + pow(2,-2) + pow(2,-3)) * (individu[3]*pow(2,-1) + individu[4]*pow(2,-2) + individu[5]*pow(2,-3))
	return x1,x2

## Generate Chromosome Lists
def decodeChromosomeList(population,rminx1, rmax1, rminx2, rmax2):
	cx1 = []
	cx2 = []
	for individu in population:
		x1,x2 = decodeChromosome(rminx1,rmax1,rminx2,rmax2,individu)
		cx1.append(x1)
		cx2.append(x2)
	return cx1, cx2

# Fitness Calculation with Minimizing Function
def fitnessCalc(x1,x2):
	fit = 1/(((np.cos(x1))*(np.sin(x2)))-(x1/(pow(x2,2)+1))+0.0001)
	return fit

def allFitness(jIndividu, cx1, cx2):
	fitnessAll = []
	i = 0
	while i < jIndividu:
		fitnessAll.append(fitnessCalc(cx1[i],cx2[i]))
		i += 1
	return fitnessAll

# Fitness Probability
def fitnessCalcProb(fitness):
	fitProb = []
	i = 0
	while i < len(fitness):
			fitProb.append(fitness[i]/np.sum(fitness))
			i += 1
	return fitProb

def fitnessCalcProbAll(fitness):
	fitProb = fitnessCalcProb(fitness)
	fitProbList = []
	fp = 0
	i = 0
	while i < len(fitness):
		fp += fitProb[i]
		fitProbList.append(fp)
		i += 1
	return fitProbList

def randomFitness(fitness):
	rnd = []
	i = 0
	while i < len(fitness):
		rnd.append(random.random())
		i += 1
	return rnd

# Parent Selection using Roulette Wheel
def rouletteWheel(fitness, population, random_):
	fp = fitnessCalcProb(fitness)
	fpAll = fitnessCalcProbAll(fitness)
	i = 0
	while i < len(population):
		if i > 0 and random_ > fpAll[i-1] and random_ < fpAll[i]:
			parents = population[i]
			break
		elif random_ < fpAll[0]:
			parents = population[i]
			break
		i += 1
	return parents

# Crossover
def crossover(parentsAll,pc):
        childAll = []
        j = 1
        i = 0
        r = random.random()
        if r<pc:
                while i < len(parentsAll):
                        pointC = random.randint(1,jKromosom-1)
                        child1 = parentsAll[i]
                        child2 = parentsAll[j]
                        childAll.append(child1[:pointC] + child2[pointC:])
                        childAll.append(child2[:pointC] + child1[pointC:])
                        i += 2
                        j += 2
        return childAll

# Mutation
def mutation(childAll,pm):
        mutationResult = childAll
        i = 0
        while i < len(mutationResult):
                if random.random() < pm:
                        rando = random.randint(1,jKromosom-1)
                        if mutationResult[i][rando] == 1: mutationResult[i][rando] = 0
                        else: mutationResult[i][rando] = 1
                i += 1
        return mutationResult 

def finalResult(maxim):
        k = 0
        while k < len(fitnessAll):
                if (fitnessAll[k] > maxim):
                        maxim = fitnessAll[k]
                        hasil = k
                k += 1
        if hasil:
                print("<<<<< RESULT >>>>>")
                print("Best Chromosome\t= ",population[hasil])
                print("Best Fitness\t= ",fitnessAll[hasil])
                print("Best x1\t\t= ",x1[hasil],"\t\tBest x2\t= ", x2[hasil])
        else:
                print("Please Try Again")
# MAIN
jKromosom = 6 ## Jumlah Kromosom
jIndividu = 40 ## Jumlah Individu
generations = 115 ## Jumlah Generasi
## Probability
pm = 0.9
pc = 0.69
## Nilai min dan max dari x1 dan x2
rminx1 = -1 
rmax1 = 2
rminx2 = -1
rmax2 = 1
# Create Population
populesyen = createPopulation(jKromosom, jIndividu)
# Decode Chromosome
x1, x2 = decodeChromosomeList(populesyen,rminx1, rmax1, rminx2, rmax2);
# Create list for all fitness value
fitnessAll = []

for i in range(generations): # Stopping Criteria (Max Generations)
	parents = [] # Create Parents List
	fitnessNew = [] # Create list for Final Result of Fitness
	population = populesyen # Assigned Population
	fitnessAll = allFitness(jIndividu, x1, x2) # Assigned all Fitness value
	rand = randomFitness(fitnessAll) # Random variable for Roulette Wheel
	j = 0
	while j<jIndividu: # Parent Selection
		parents.append(rouletteWheel(fitnessAll, population,rand[j]))
		j += 1
	childAll = crossover(parents,pc) # Generate Child
	mutationL = mutation(childAll, pm) # Mutation Start
	x1, x2 = decodeChromosomeList(population,rminx1, rmax1, rminx2, rmax2)
	fitnessNew = allFitness(len(mutationL),x1,x2)

maxim = fitnessAll[0]
finalResult(maxim)
