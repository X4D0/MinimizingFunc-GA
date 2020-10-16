import numpy as np
import random

# Desain Kromosom dan Metode Pendekodean

## Binary Encoding using 3 bits (3 Gens)
def createPopulation(jumlahGen,jumlahKromosom):
        pop = []
        i = 0
        while i < jumlahKromosom:
                gen = []
                j = 0
                while j<jumlahGen:
                        gen.append(random.randint(0,1))
                        j += 1
                pop.append(gen)
                i += 1
        return pop

## Decode Chromosome
def decodeChromosome(rminx1, rmax1, rminx2, rmax2,kromosom):
        x1 = rminx1 + (rmax1 - rminx1)/(pow(2,-1) + pow(2,-2) + pow(2,-3)) * (kromosom[0]*pow(2,-1) + kromosom[1]*pow(2,-2) + kromosom[2]*pow(2,-3))
        x2 = rminx2 + (rmax2 - rminx2)/(pow(2,-1) + pow(2,-2) + pow(2,-3)) * (kromosom[3]*pow(2,-1) + kromosom[4]*pow(2,-2) + kromosom[5]*pow(2,-3))
        return x1,x2

## Generate Chromosome Lists
def decodeChromosomeList(population,rminx1, rmax1, rminx2, rmax2):
        cx1 = []
        cx2 = []
        for kromosom in population:
                x1,x2 = decodeChromosome(rminx1,rmax1,rminx2,rmax2,kromosom)
                cx1.append(x1)
                cx2.append(x2)
        return cx1, cx2

# Fitness Calculation with Minimizing Function
def fitnessCalc(x1,x2):
        fit = 1/(((np.cos(x1))*(np.sin(x2)))-(x1/(pow(x2,2)+1))+0.0001)
        return fit

def allFitness(jKromosom, cx1, cx2):
        fitnessAll = []
        i = 0
        while i < jKromosom:
                fitnessAll.append(fitnessCalc(cx1[i],cx2[i]))
                i += 1
        return fitnessAll

# Parent Selection using Roulette Wheel
def rouletteWheel(fitness, population, random_):
        fitProb = []
        fp = 0
        i = 0
        while i < len(fitness):
                fp += fitness[i]/np.sum(fitness)
                fitProb.append(fp)
                i += 1
        r = 0
        while r < len(population):
                if r > 0 and fitProb[r-1] < random_ and fitProb[r] > random_:
                        indiv = population[r]
                        break
                elif random_ < fitProb[0]:
                        indiv = population[r]
                        break
                r += 1
        return indiv

# Crossover
def crossover(parentsAll,pc):
        childAll = []
        j = 1
        i = 0
        while i < len(parentsAll):
                if np.random.uniform(0,1)<pc:
                        titik = random.randint(1,jGens-1)
                        child1 = parentsAll[i]
                        child2 = parentsAll[j]
                        childAll.append(child1[:titik] + child2[titik:])
                        childAll.append(child2[:titik] + child1[titik:])
                        i += 2
                        j += 2
        return childAll

# Mutation
def mutation(childAll,pm):
        mutationResult = childAll
        i = 0
        while i < len(mutationResult):
                if random.random() < pm:
                        rando = random.randint(1,jGens-1)
                        if mutationResult[i][rando] == 1: mutationResult[i][rando] = 0
                        else: mutationResult[i][rando] = 1
                i += 1
        return mutationResult 

# Steady-State Fitness-Based Procedures
def steadyState(fitpop,fitmut):
        baru = []
        i = 0
        a = min(fitpop)
        b = min(fitmut)
        for fitness in fitpop:
                if a<b:
                        baru.append(b)
                elif a>b:
                        baru.append(a)
                else:
                        baru.append(fitness)
        return baru

def finalResult(maxim):
        k = 0
        while k < len(fitnessAll):
                if (fitnessAll[k] > maxim):
                        maxim = fitnessAll[k]
                        hasil = k
                k += 1
        res = hasil
        print("\nResult :\n")
        print("Best Chromosome\t= ",population[res])
        print("Best Fitness\t= ",fitnessAll[res])
        print("Best x1\t\t= ",x1[res],"\t\tBest x2\t= ", x2[res])
        
# MAIN

jGens = 6 ## Jumlah Gen
jKromosom = 40 ## Jumlah kromosom
generations = 115 ## Jumlah Generasi
## Probability
pm = 0.9
pc = 0.69
## Nilai min dan max dari x1 dan x2
rminx1 = -1 
rmax1 = 2
rminx2 = -1
rmax2 = 1
# Create First Population
populesyen = createPopulation(jGens, jKromosom)
# First Decode Chromosome
x1, x2 = decodeChromosomeList(populesyen,rminx1, rmax1, rminx2, rmax2)
# Create list for all fitness value
fitnessAll = []

for i in range(generations): # Stopping Criteria (Max Generations)
        parents = [] # Create Parents List
        fitnessNew = [] # Create list for Final Result of Fitness
        population = populesyen # Assigned Population
        fitnessAll = allFitness(jKromosom, x1, x2) # Assigned all Fitness value
        j = 0
        while j<jKromosom: # Parent Selection
                parents.append(rouletteWheel(fitnessAll,population,random.random()))
                j += 1
        childAll = crossover(parents,pc) # Do Crossover
        mutasi = mutation(childAll, pm) # Do Mutation
        x1,x2 = decodeChromosomeList(population,rminx1, rmax1, rminx2, rmax2)
        fitnessNew = allFitness(len(mutasi),x1,x2)
        fitnessAll = steadyState(fitnessAll,fitnessNew)
        print("Generation ",i+1," DONE!")

maxim = fitnessAll[0]
finalResult(maxim)
