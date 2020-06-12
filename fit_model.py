import numpy
from random import random
from sys import argv

from utils import today, day_zero

### EVALUATION FUNCTIONS
def propagate(days, N, I0, beta0, ommega, gammaD, gammaR):
    susceptible = [N]
    active = [I0]
    deaths = [0]
    recoveries = [0]
    dailyCases = []
    dailyDeaths = []
    beta = beta0
    for d in range(1,days+1):
        dsdt = -beta * susceptible[d-1] * active[d-1] / N
        didt = beta * susceptible[d-1] * active[d-1] / N - gammaD * active[d-1] - gammaR * active[d-1]
        dddt = gammaD * active[d-1]
        drdt = gammaR * active[d-1]
        beta = beta0 * numpy.exp(-ommega * d)
        susceptible.append(susceptible[-1] + dsdt)
        active.append(active[-1] + didt)
        deaths.append(deaths[-1] + dddt)
        recoveries.append(recoveries[-1] + drdt)
        dailyCases.append(-dsdt)
        dailyDeaths.append(dddt)
    susceptible = numpy.array(susceptible)
    active = numpy.array(active)
    deaths = numpy.array(deaths)
    recoveries = numpy.array(recoveries)
    return [susceptible, active, deaths, recoveries, dailyCases, dailyDeaths]

def error(data_active, data_deaths, data_recoveries, _active, _deaths, _recoveries):
    # Crop future data if there is
    tiltoday = len(data_active)
    if len(_active) > tiltoday:
        _active = _active[:tiltoday]
        _deaths = _deaths[:tiltoday]
        _recoveries = _recoveries[:tiltoday]
    SEC = sum(numpy.square(data_active-_active))
    SEC += sum(numpy.square(data_deaths - _deaths))
    SEC += sum(numpy.square(data_recoveries - _recoveries))
    SEC /= len(_active)
    return SEC
    
### INITIAL SETUP
popSize = 5000
generations = 100
new_population = []

## Create initial population
print("Creating initial population of {}".format(popSize))
for i in range(popSize):
    N = 500000 + int(random() * 500000)
    I0 = int(random() * 5000)
    beta0 = random()
    ommega = random()
    gammaD = random()
    gammaR = random()
    new_population.append([N, I0, beta0, ommega, gammaD, gammaR])

## Read data
print("Reading data from file {}".format(argv[1]))
dataA = []
dataD = []
dataR = []
dataDC = []
dataDD = []
dataTotal = []
fopen = open(argv[1])
days = 0
for line in fopen:
    days, total, active, deaths, recoveries, dailyC, dailyD, dailyR = list(map(float,line.strip().split(',')))
    dataA.append(active)
    dataD.append(deaths)
    dataR.append(recoveries)
    dataDC.append(dailyC)
    dataDD.append(dailyD)
    dataTotal.append(total)
    ## Remove last days of data
    if days == today() - day_zero - 12: break
days = int(days)
fopen.close()
dataA = numpy.array(dataA)
dataD = numpy.array(dataD)
dataR = numpy.array(dataR)


### EVOLVE
print("Evolving...")
best = []
try:
    for gen in range(generations):
        ## Evaluate
        fitness = numpy.zeros(popSize)
        population = list(new_population)
        for i in range(len(population)):
            specimen = population[i]
            N = specimen[0]
            I0 = specimen[1]
            beta0 = specimen[2]
            ommega = specimen[3]
            gammaD = specimen[4]
            gammaR = specimen[5]
            susceptible, active, deaths, recoveries, dc, dd = propagate(days, N, I0, beta0, ommega, gammaR, gammaD)
            fitness[i] = error(dataA, dataD, dataR, active, deaths, recoveries)

        ## Select best 5
        fitness = list(fitness)
        fittest = []
        top_specimen = numpy.format_float_scientific(min(fitness), precision=4)
        for i in range(5):
            fit_idx = fitness.index(min(fitness))
            fittest.append(fit_idx)
            fitness.remove(min(fitness))
        print("Generation {} fitness (error) {}".format(gen,top_specimen))
        # Save best fit from generation to animate evolution later
        best.append(population[fittest[0]])
        ## Procreate
        new_population = []
        for i in range(popSize):
            new_specimen = []
            # Inherit trait j from either 5 top or random mutation based on the best one
            for j in range(6):
                chances = random()
                if chances < 1/6:
                    new_specimen.append(population[fittest[0]][j])
                elif chances < 2/6:
                    new_specimen.append(population[fittest[1]][j])
                elif chances < 3/6:
                    new_specimen.append(population[fittest[2]][j])
                elif chances < 4/6:
                    new_specimen.append(population[fittest[3]][j])
                elif chances < 5/6:
                    new_specimen.append(population[fittest[4]][j])
                else:
                    new_specimen.append(population[fittest[0]][j] * (0.5 + random()))
            new_population.append(new_specimen)
except KeyboardInterrupt:
    print("\nInterrupting evolution process")

specimen = best[-1]
N = specimen[0]
I0 = specimen[1]
beta0 = specimen[2]
ommega = specimen[3]
gammaD = specimen[4]
gammaR = specimen[5]
susceptible, active, deaths, recoveries, dc, dd= propagate(days, N, I0, beta0, ommega, gammaR, gammaD)

print("Fittest data:")
print(specimen)
print("Plotting")
import matplotlib.pyplot as plot
tim = list(range(days+1))
plot.plot(tim,active)
plot.plot(tim,dataA)
plot.grid(which='major',alpha=0.5)
plot.grid(which='minor',alpha=0.2)
plot.legend(("Model","Data"))
plot.show()

print("Forecast")
extra_days = 300
susceptible, active, deaths, recoveries, dc, dd = propagate(extra_days, N, I0, beta0, ommega, gammaR, gammaD)
tim_full = list(range(extra_days+1))

plot.plot(tim_full,active)
plot.plot(tim, dataA)
plot.plot(tim_full, deaths)
#plot.plot(tim_full, recoveries)
plot.grid(which='major',alpha=0.5)
plot.grid(which='minor',alpha=0.2)
plot.legend(("Prediccion de infecciones","Datos","Muertes"))
plot.plot(tim, dataD, color='r')
plot.plot(tim_full, N - susceptible)
plot.plot(tim, dataTotal, color='r')
plot.xlabel("Dias desde Feb 22")
# Draw a line to mark today
today_x = today() - day_zero
top_y = active[today_x]
plot.plot([today_x,today_x], [0,top_y], color='g')
plot.show()
plot.scatter(tim, dataDC, color='r')
plot.scatter(tim, dataDD, color='k')
plot.plot(tim_full, [0]+dc, color='r')
plot.plot(tim_full, [0]+dd, color='k')
plot.legend(("Casos","","Muertes",""))
top_y = dc[today_x]
plot.plot([today_x,today_x], [0,top_y], color='g')
plot.grid(which='major',alpha=0.5)
plot.grid(which='minor',alpha=0.2)
plot.show()
mulplt = plot.figure(1)
anmplt = plot.figure(2)
print("Generation animation frames of evolution process")
for i in range(len(best)):
    specimen = best[i]
    N = specimen[0]
    I0 = specimen[1]
    beta0 = specimen[2]
    ommega = specimen[3]
    gammaD = specimen[4]
    gammaR = specimen[5]
    susceptible, active, deaths, recoveries, dc, dd= propagate(extra_days, N, I0, beta0, ommega, gammaR, gammaD)
    # Save plot image
    plot.figure(1)
    plot.plot(tim_full,active)
    plot.plot(tim, dataA)
    plot.plot(tim_full, deaths)
    plot.grid()
#    plot.legend(("Prediccion de infecciones","Datos","Muertes"))
    plot.xlabel("Dias desde Feb 22")
    plot.savefig("imgs/best{}".format(i),format='png')
    plot.cla()
    err = error(dataA, dataD, dataR, active, deaths, recoveries)
    if err < 1e7:
        plot.figure(2)
        plot.plot(tim_full,active, alpha=0.033, color='b')
        plot.plot(tim_full, deaths, alpha=0.033, color='k')
plot.figure(2)
plot.plot(tim, dataA, color='r')
top_y = active[today_x]
plot.plot([today_x,today_x], [0,top_y], color='g')
plot.plot(tim, dataD)
plot.grid(which='major',alpha=0.5)
plot.grid(which='minor',alpha=0.2)
plot.show()
