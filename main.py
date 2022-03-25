from bitarray import bitarray
import random as rd
import math
import matplotlib.pyplot as plt1
from pylab import *
import colorsys

# rd.seed(566)
# 445
# 50
numeracids =12
punishment = 1
bits = 6*numeracids
viruses = 8

islands = 2
mutationRate = .3
maxMuation = 10
initPop = 50
make = 50
keepBest=10
keepnext =50
shareInfo = False
best = True
acug = {'a':bitarray('00'),
            'c':bitarray('01'),
            'u':bitarray('10'),
            'g':bitarray('11')}
acugbit = { bitarray('00').to01():'a',
            bitarray('01').to01():'c',
            bitarray('10').to01():'u',
            bitarray('11').to01():'g'}

ala = {'gcu', 'gcc', 'gca', 'gcg'}
arg = {'cgu', 'cgc', 'cga', 'cgg', 'aga', 'agg'}
asn = {'aau', 'aac'}
asp = {'gau', 'gac'}
cys = {'ugu', 'ugc'}
gln = {'caa', 'cag'}
glu = {'gaa', 'gag'}
gly = {'ggu', 'ggc', 'gga', 'ggg'}
his = {'cau', 'cac'}
lle = {'auu', 'auc', 'aua'}
start = {'aug'}
leu = {'uua', 'uug', 'cuu', 'cuc', 'cua', 'cug'}
lys = {'aaa', 'aag'}
met = {'aug'}
phe = {'uuu', 'uuc'}
pro = {'ccu', 'ccc', 'cca', 'ccg'}
ser = {'ucu', 'ucc', 'uca', 'ucg', 'agu', 'agc'}
thr = {'acu', 'acc', 'aca', 'acg'}
trp = {'ugg'}
tyr = {'uau', 'uac'}
val = {'guu', 'guc', 'gua', 'gug'}
stop = {'uag', 'uga', 'uaa'}

all = [ala,arg,asn,asp,cys,gln,glu,gly,his,lle,start,leu,lys,met,phe,pro,ser,thr,trp,tyr,val,stop]
mapacug= {}

for x in all:
    for y in x:
        mapacug[y]=x

def test111(a, b):
    if (not ((a) in mapacug[b])):
        return min(
            list(map(lambda z: len(list(filter(lambda x: not (x[0] == x[1]), zip(list(a), list(z))))), mapacug[b])))
    else:
        return 0



class bitArray:
    counter = 0

    def __init__(self, size, fitness=100000000000000000000, valid=False, array=None, iteration=0, parent=0, id=0,
                 parentid=0):

        bitArray.counter += 1
        self.fitness = fitness
        self.valid = valid
        self.size = size
        self.iteration = iteration

        self.parentid = parentid
        if (id == 0):
            self.id = bitArray.counter
        else:
            self.id = id

        if (parent == 0):
            self.parent = bitArray.counter
        else:
            self.parent = parent

        if array == None:
            self.array = self.makeArray(size)
        else:
            self.array = array

        # print(self.array)

    def makeArray(self, size):
        zero_count = rd.randint(0, size)
        one_count = size - zero_count
        my_list = [0] * zero_count + [1] * one_count
        rd.shuffle(my_list)
        return bitarray(''.join(map(str, my_list)))

    def flip(self, iteration):
        self.iteration = iteration
        if self.valid:
            bitArray.counter += 1
            self.parentid = self.id
            self.id = bitArray.counter

        self.valid = False
        self.fitness = 10000000000000000
        loction = rd.randint(0, self.size / 2 - 1)
        self.array[loction * 2:loction * 2 + 2] = self.makeArray(2)

        return self

    def calcFitness(self, others):
        return self.calcFitness1(others)

    def calcFitness1(self, others):

        self.valid = True
        y = [acugbit[self.array[i:i + 2].to01()] for i in range(0, bits, 2)]
        z = [''.join(y[i:i + 3]) for i in range(0,  len(y), 3)]
        print(list(map(lambda x: sum(list(map(lambda q: test111((q[1]), q[0]), zip(x, z)))), others)))
        arrayMin = min(
            zip(list(map(lambda x: sum(list(map(lambda q: test111((q[1]), q[0]), zip(x, z)))), others)),
                range(viruses)),
            key=lambda x: x[0])
        self.fitness = arrayMin[0]
        return self

    def clone(self):
        return bitArray(self.size, fitness=self.fitness, valid=self.valid, array=self.array.copy(),
                        iteration=self.iteration, parent=self.parent, id=self.id, parentid=self.parentid)

    def __str__(self):
        return self.array


def makePopulation(amount, size):
    return [bitArray(size) for _ in range(amount)]


virus = makePopulation(viruses, bits)
bitVirus = list(map(lambda y:  [ ''.join( y[i:i+3]) for i in range(0, len(y), 3)],  map(lambda x: [acugbit[x.array[i:i + 2].to01()] for i in range(0, bits, 2)], virus)))

populations = [makePopulation(initPop, bits) for _ in range(islands)]

# avgList =list()
avgList = [list() for _ in range(islands)]
actions = [set() for _ in range(islands)]
infolist = [list() for _ in range(islands)]


def pprint(pop,island):
    par = len(set([ind.parent for ind in pop]))
    fits = [ind.fitness for ind in pop]
    length = len(pop)
    mean = sum(fits) / length
    sum2 = sum(x * x for x in fits)
    std = abs(sum2 / length - mean ** 2) ** 0.5
    min_ = min(fits)
    max_ =max(fits)
    print("  Min %s" % min_)
    print("  Max %s" % max_)
    print("  Avg %s" % mean)
    print("  Std %s" % std)
    print("  unique blood lines %s" % par)
    infolist[island].append((min_,max_,mean))
    return mean


def tick(pop, count, popNumber):
    newClones = [ind.clone() for ind in pop]

    buildinglist = list()
    buildinglist+=newClones[0:keepBest]
    while len(buildinglist) < make:
        buildinglist.append(rd.choice(newClones))

    for mutant in buildinglist:
        counter = 0
        while (rd.random() < mutationRate) and counter < maxMuation:
            mutant.flip(count)
            counter += 1
    # bcell = 0


    for i in list(filter(lambda x: not x.valid, buildinglist)):
        i.calcFitness(bitVirus)
        # bcell += 1
        actions[popNumber].add((i.parentid, i.id, count, i.fitness, i.parent))

    mean = pprint(buildinglist,popNumber)

    avgList[popNumber].append(mean)
    # offspring = sorted((filter(lambda x: x.fitness <= mean, buildinglist)), key=lambda x: x.fitness)
    offspring = sorted(buildinglist, key=lambda x: x.fitness)
    return offspring[0:keepnext]


it = 0
fitsMin = 1000000000

while it < 42 and fitsMin != 0:
    populations = [tick(pop, it, popNumber) for pop, popNumber in zip(populations, range(islands))]
    it += 1
    # if(shareInfo):
    #     populations = [tick(pop, it) for pop in populations]



# print(bloodlined)
#
# print(actions[0])

# fig =plt1.figure()

for island in range(islands):
    print(infolist[island])
    r  = range(len( infolist[island]))
    color= colorsys.hsv_to_rgb((island+1)/islands,(island+1)/islands,
                                          .5)
    plt.plot(r,list(map(lambda  x : x[0], infolist[island])),color=color,linestyle='--' )
    plt.plot(r,list(map(lambda  x : x[1], infolist[island])),color=color,linestyle='--' )

    plt.fill_between(r,list(map(lambda  x : x[0], infolist[island])), list(map(lambda  x : x[1], infolist[island])),color=color,alpha=0.05)
    # ax.fill_between(
    #     time, income, expenses, where=(income <= expenses),
    #     interpolate=True, color="red", alpha=0.25,
    #     label="Negative"
    # )
    plt.plot(r,list(map(lambda  x : x[2], infolist[island])),color=color  )
plt.show()




fig = plt.figure()
ax = fig.gca(projection='3d')
for island in range(islands):
    best_ind = min(list(filter(lambda y: y[2] == it - 1, actions[island])), key=lambda x: x[3])
    # avgList.append(best_ind[4])
    bloodlined = list()
    last = best_ind
    bloodlined.append(last)
    for line in reversed(sorted(actions[island])):
        if (last[0] == line[1]):
            bloodlined.append(line)
            last = line
    yy = set(filter(lambda y: best_ind[4] == y[4], actions[island]))

    if(best):
        for x in actions[island]:
            for y in actions[island]:
                if  not( y in yy) and  not( x in yy):
                    if 0 == y[0]:
                        ax.scatter(y[3], avgList[island][y[2]], y[2])
                    if y[1] == x[0]:
                        if y[3] > x[3]:
                            ax.plot3D([y[3], x[3]],
                                      [avgList[island][y[2]], avgList[island][x[2]]],
                                      [y[2], x[2]],
                                      color=(.5, 0, .5), linewidth=1)
                        else:
                            ax.plot3D([y[3], x[3]],
                                      [avgList[0][y[2]], avgList[0][x[2]]],
                                      [y[2], x[2]],
                                      color=(0, .5, .5), linewidth=1)
        for x in yy:
            for y in yy:
                # last[0] == line[1]
                if y[0] == x[1]:
                    if y[3] > x[3]:
                        ax.plot3D([y[3], x[3]],
                                  [avgList[island][y[2]], avgList[island][x[2]]],
                                  [y[2], x[2]],
                                  color=(1, 0, 1), linewidth=.5)
                    else:
                        ax.plot3D([y[3], x[3]],
                                  [avgList[island][y[2]], avgList[island][x[2]]],
                                  [y[2], x[2]],
                                  color=(0, 1, 1), linewidth=.5)

    # print(len(avgList))
    for x in yy:
        for y in yy:
            if y[0] == x[1]:
                if (bloodlined.__contains__(x) and bloodlined.__contains__(y)):
                    ax.plot3D([y[3], x[3]],
                              [avgList[island][y[2]], avgList[island][x[2]]],
                              [y[2], x[2]],
                              color=(0, 0, 1))

ax.set_ylabel('Avg fitness')
ax.set_zlabel('iterations')
ax.set_xlabel('Real fitness')

ax.invert_zaxis()
plt.show()