from sys import argv
from numpy import zeros, array, histogram

from utils import *

recoveryThreshold = 16

# We need a csv file with the following data:
# + Days from start of epidemic (2020-02-22)
# + Total cases
# + Total deaths
# + Total recoveries (Estimated with first symptoms date + 14 days)
# + Daily new cases
# + Daily new deaths
# + Daily new recoveries


fname = argv[1]

fopen = open(fname, encoding= 'unicode_escape')


fields = fopen.readline().strip().replace('"','').split(',')

total_days = totalDays()

totalCases = zeros(total_days)
totalDeaths = zeros(total_days)
totalRecoveries = zeros(total_days)

dailyCases = zeros(total_days)
dailyDeaths = zeros(total_days+1)
dailyRecoveries = zeros(total_days)

daysUntilDeath = []

uciDeaths = 0
uciRecoveries = 0
totalUci = 0

diabetesDeaths = 0
diabetesRecoveries = 0
totalDiabetes = 0

epocDeaths = 0
epocRecoveries = 0
totalEpoc = 0

inmunoDeaths = 0
inmunoRecoveries = 0
totalInmuno = 0

hipertDeaths = 0
hipertRecoveries = 0
totalHipert = 0

obesDeaths = 0
obesRecoveries = 0
totalObes = 0

tabacoDeaths = 0
tabacoRecoveries = 0
totalTabaco = 0

_today = today()

#states = ['09','15'] # Take only CDMX('09') and EDOMEX('15') ( use 'all' for every state)
states = ['all']

print("Processing data...\n")
print("Taking data from states:")
for state in states: print(" - "+state)

print("Using recovery threshold of {} days".format(recoveryThreshold))


for line in fopen:
    l = line.strip().replace('"','').split(',')
    
    if l[30] != '1': continue # Not confirmed CoVid. Next line
    if l[6] not in states and states[0] != 'all': continue
    
    daySymptoms = dayOffset(dayOfYear(readDate(l[11])))
    
    if daySymptoms < 0: continue # Symptoms before Feb 22? Skip

    ## Has other comorbidities?
    uci = l[34] == '1'
    diabetes = l[19] == '1'
    epoc = l[20] == '1'
    inmuno = l[22] == '1'
    hipert = l[23] == '1'
    obes = l[26] == '1'
    tabaco = l[28] == '1'


    if l[12] != "9999-99-99": # Died
        dayDeath = dayOffset(dayOfYear(readDate(l[12])))
        lasted = dayDeath - daySymptoms
        if lasted < 0: # Died before showing symptoms? Ignore
            #print(l)
            continue
        daysUntilDeath.append(lasted)
        try:dailyDeaths[dayDeath] += 1
        except:continue
        if uci:
            uciDeaths += 1
        if diabetes:
            diabetesDeaths += 1
        if epoc:
            epocDeaths += 1
        if inmuno:
            inmunoDeaths += 1
        if hipert:
            hipertDeaths += 1
        if obes:
            obesDeaths += 1
        if tabaco:
            tabacoDeaths += 1

    elif _today - daySymptoms > recoveryThreshold: # Days until recovery passed
        try:
            dailyRecoveries[daySymptoms + recoveryThreshold] += 1
            if uci:
                uciRecoveries += 1
            if diabetes:
                diabetesRecoveries += 1
            if epoc:
                epocRecoveries += 1
            if inmuno:
                inmunoRecoveries += 1
            if hipert:
                hipertRecoveries += 1
            if obes:
                obesRecoveries += 1
            if tabaco:
                tabacoRecoveries += 1
        except:
            pass
    dailyCases[daySymptoms] += 1
    if uci:
        totalUci += 1
    if diabetes:
        totalDiabetes += 1
    if epoc:
        totalEpoc += 1
    if inmuno:
        totalInmuno += 1
    if hipert:
        totalHipert += 1
    if obes:
        totalObes += 1
    if tabaco:
        totalTabaco += 1
fopen.close()

# Now accumulate daily data
for i in range(len(dailyCases)):
    totalCases[i] = totalCases[i-1] + dailyCases[i]
    totalDeaths[i] = totalDeaths[i-1] + dailyDeaths[i]
    totalRecoveries[i] = totalRecoveries[i-1] + dailyRecoveries[i]

# Export data
print("Exporting data to transformed_data.csv\n\n")
output = open("transformed_data.csv",'w')
output.write("")
for i in range(total_days):
    activeCases = totalCases[i] - totalDeaths[i] - totalRecoveries[i]
    output.write("{},{},{},{},{},{},{},{}\n".format(i,totalCases[i],activeCases,totalDeaths[i],totalRecoveries[i],\
                                                dailyCases[i],dailyDeaths[i],dailyRecoveries[i]))
output.close()
#
print("Days until death")
daysUntilDeath_hist = histogram(array(daysUntilDeath),70)
print("Min: {}".format(min(daysUntilDeath)))
print("Max: {}".format(max(daysUntilDeath)))
print("Avg: {:.1f}".format(sum(daysUntilDeath) / len(daysUntilDeath)))
print ("Total deaths: {}".format(len(daysUntilDeath)))
print()
# UCI
print("UCI")
print("Total admitted: {}".format(totalUci))
print("Died at UCI: {:0.1f}%".format(100*uciDeaths/totalUci))
print("Recovered after UCI: {:0.1f}%".format(100*uciRecoveries/totalUci))
# DIABETES
print("\nDIABETES")
print("Total admitted: {}".format(totalDiabetes))
print("Died with DIABETES: {:0.1f}%".format(100*diabetesDeaths/totalDiabetes))
print("Recovered with DIABETES: {:0.1f}%".format(100*diabetesRecoveries/totalDiabetes))
# EPOC
print("\nEPOC")
print("Total admitted: {}".format(totalEpoc))
print("Died with EPOC: {:0.1f}%".format(100*epocDeaths/totalEpoc))
print("Recovered with EPOC: {:0.1f}%".format(100*epocRecoveries/totalEpoc))
# INMUNO
print("\nINMUNO")
print("Total admitted: {}".format(totalInmuno))
print("Died with INMUNO: {:0.1f}%".format(100*inmunoDeaths/totalInmuno))
print("Recovered with INMUNO: {:0.1f}%".format(100*inmunoRecoveries/totalInmuno))
# HIPERTENSION
print("\nHIPERTENSION")
print("Total admitted: {}".format(totalHipert))
print("Died with HIPERTENSION: {:0.1f}%".format(100*hipertDeaths/totalHipert))
print("Recovered with HIPERTENSION: {:0.1f}%".format(100*hipertRecoveries/totalHipert))
# OBESIDAD
print("\nOBESIDAD")
print("Total admitted: {}".format(totalObes))
print("Died with OBESIDAD: {:0.1f}%".format(100*obesDeaths/totalObes))
print("Recovered with OBESIDAD: {:0.1f}%".format(100*obesRecoveries/totalObes))
# TABACO
print("\nTABACO")
print("Total admitted: {}".format(totalTabaco))
print("Died with TABACO: {:0.1f}%".format(100*tabacoDeaths/totalTabaco))
print("Recovered with TABACO: {:0.1f}%".format(100*tabacoRecoveries/totalTabaco))
# Plot
print("Plotting...")
import matplotlib.pyplot as plot
# Days until death
#plot.hist(daysUntilDeath, 50)
plot.bar(daysUntilDeath_hist[1][:-1], 100*daysUntilDeath_hist[0]/len(daysUntilDeath))
plot.grid()
plot.title("Días desde primeros síntomas hasta muerte")
plot.xlabel("Días")
plot.ylabel("Porcentaje (%)")
plot.show()
# Total cases/deaths/recoveries
plot.plot(range(total_days), totalCases)
plot.plot(range(total_days), totalDeaths)
plot.plot(range(total_days), totalRecoveries)
plot.legend(("Cases","Deaths","Recoveries"))
plot.grid()
plot.title("Total cases/deaths/recoveries")
plot.xlabel("Days from Feb 22")
plot.ylabel("Number of cases")
plot.show()
# Daily cases/deaths/recoveries
plot.plot(range(total_days), dailyCases)
plot.plot(range(total_days+1), dailyDeaths)
plot.plot(range(total_days), dailyRecoveries)
plot.legend(("Cases","Deaths","Recoveries"))
plot.grid()
plot.title("Daily cases/deaths/recoveries")
plot.xlabel("Days from Feb 22")
plot.ylabel("Number of cases")
plot.show()
