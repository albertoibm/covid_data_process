from sys import argv
from random import random

from numpy import array, log, exp, zeros, abs
from scipy import optimize

from utils import *

print("CoVid19 mortality predictor")

recoveryThreshold = 16

fname = argv[1]

fopen = open(fname, encoding= 'unicode_escape')

fields = fopen.readline().strip().replace('"','').split(',')

# x vector
# x         =   [1, sexo, edad, diabetes, epoc, asma, inmuno, hipert, otra, cardio, obes, renal, tabaco, otro]

X = []
y = []

_today = today() - day_zero
print("Reading data...")
for line in fopen:
    l = line.strip().replace('"','').split(',')
    # Not confirmed CoVid. Next line
    if l[30] != '1': continue     
    # Murio (1=Si, 0=No)
    _y = (l[12] != "9999-99-99") * 1   
    # Drop surviving cases randomly in order to have a better dataset
    #if _y == 0:
    #    if random() < 2.1/6: continue
    # Day of first symptoms (Is recovered?)
    daySymptoms = dayOffset(dayOfYear(readDate(l[11])))
    # Alive and not recovered, skip
    if _y == 0 and _today - daySymptoms < recoveryThreshold: 
        continue
    # Else continue and collect data
    _x = [1]
    _x.append( (l[5] == '2') * 1.) # Sexo (1=H, 0=M)
    _x.append( int(l[15]) ) # Edad
    _x.append( (l[19] == '1') * 1.) # Diabetes (1=Si, 0=No)
    _x.append( (l[20] == '1') * 1.) # EPOC
    _x.append( (l[21] == '1') * 1.) # Asma
    _x.append( (l[22] == '1') * 1.) # Inmunosupresion
    _x.append( (l[23] == '1') * 1.) # Hipertension
    _x.append( (l[24] == '1') * 1.) # Otra complicacion
    _x.append( (l[25] == '1') * 1.) # Cardiovascular
    _x.append( (l[26] == '1') * 1.) # Obesidad
    _x.append( (l[27] == '1') * 1.) # Renal Cronica
    _x.append( (l[28] == '1') * 1.) # Tabaco
    _x.append( (l[29] == '1') * 1.) # Otro caso
    X.append(_x)
    y.append(_y) 

fopen.close()
## Normalize
X = array(X)
y = array(y)
_mean = X.mean(axis=0)
_mean[0] = 0 # Keep our first element of X equal to 1 as is the bias multiplier
_std = X.std(axis=0)
_std[0] = 1 # Same
X -= _mean
X /= _std

print("Dataset has a {:.2f} recovered per dead ratio".format((len(y)-y.sum())/y.sum()))
## Start learning
print("Learning from data...")
## Prepare learning and testing sets
## Use 70% for learning, 30% for testing
limit = int(len(X) * 0.7)
X_learn = array(list(X[:limit]))
y_learn = array(list(y[:limit]))

X_test = array(list(X[limit:]))
y_test = array(list(y[limit:]))

print("Total cases: {}".format(len(y)))
print("Size of training set: {}".format(len(y_learn)))
print("Size of validation set: {}".format(len(y_test)))

## Minimize cost function
theta = array([ -2, .1, 1.1, .2, .2, .2, .1, .2, .1, .1, .1, .1, .2, .1])
       
from sklearn.linear_model import LogisticRegression

predictor = LogisticRegression(max_iter=5000,solver='liblinear',tol=0.000001)
predictor.fit(X_learn, y_learn)


print("Coeficients:")
print(predictor.intercept_)
print(predictor.coef_)

print("\nScore on test set:")
print(predictor.score(X_test,y_test))


## Calculate AUROC (Area Under the Receiver Operating Characteristics)
## Predict with test set
print("Calculating AUROC")
res_test = predictor.predict_proba(X_test)[:,1]
thresh = []
fpr = []
fnr = []
tpr = []
tp = []
tn = []
fp = []
fn = []
thresh_n = 10000
for i in range(thresh_n):
    thresh.append(i / thresh_n)
    res = array(((res_test > thresh[-1]) * 1))
    tp.append( (res * y_test).sum() )# true positives
    tn.append( ((1 - res) * (1 - y_test)).sum()) # true negatives
    fp.append( (res * (1 - y_test)).sum() )
    fn.append( (y_test * (1 - res)).sum())
    tpr.append( tp[-1] / (tp[-1] + fn[-1]) )
    fnr.append( 1 - tpr[-1])
    fpr.append( 1 - tn[-1] / (tn[-1] + fp[-1]) )

## Integrate area under curve
A = 0
for i in range(1,len(fpr)):
    A += abs(fpr[i] - fpr[i-1]) * (tpr[i] + tpr[i-1]) / 2
print("AUC = {}".format(A))

import matplotlib.pyplot as plot
plot.title("Todos los casos")
plot.plot( thresh, tp )
plot.plot( thresh, tn )
plot.plot( thresh, fp )
plot.plot( thresh, fn )
plot.grid()
plot.legend(("Verdaderos positivos","Verdaderos negativos","Falsos positivos","Falsos negativos"))
plot.xlabel("Umbral")
plot.show()
plot.title("Tasas de falsos positivos y negativos")
plot.plot( thresh, fpr )
plot.plot( thresh, fnr )
plot.legend(("Tasa de falsos positivos","Tasa de falsos negativos"))
plot.xlabel("Umbral")
plot.grid()
plot.show()
plot.title("Curva ROC")
plot.plot(fpr, tpr)
plot.fill_between(fpr, tpr)
plot.plot([0,1],[0,1],alpha=0.2,linestyle='--',color='k')
plot.xlabel("1 - Especificidad")
plot.ylabel("Sensitividad")
plot.figtext(.60,.30,"Area: {:.2f}".format(A),color='w')
plot.grid()
plot.show()


thresh = list(array(fn) > array(fp)).index(True) / thresh_n

res = array(((res_test > thresh) * 1))
score = 1 - abs(y_test - res).sum() / len(y_test)
print("Score with threshold = {:.6f} : {:.3f}".format(thresh, score))


tp = (res * y_test).sum() # true positives
tn = ((1 - res) * (1 - y_test)).sum() # true negatives
fp = (res * (1 - y_test)).sum() 
fn = (y_test * (1 - res)).sum()

sensitivity = tp / (tp + fn)
specificity = tn / (tn + fp)


print("\nPorcentaje de positivos en datos:")
print(y_test.sum() / len(res_test))
print("\nPorcentaje de positivos en prediccion:")
print(res.sum() / len(res_test))
print("\nFalse positives:")
print(fp / (fp + tn))
print("\nFalse negatives:")
print(fn / (fn + tp))
print("\nTotal:")
print(fp / (fp + tn) + fn / (fn + tp))
print("\nSensitivity:\n{:.3f}".format(sensitivity))
print("\nSpecificity:\n{:.3f}".format(specificity))

#print("\nMean:")
#print(_mean)
#print("\nStandard dev:")
#print(_std)


# Predict me
print("\nPredicting me:")
me = array([1.,1.,30.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]).reshape(1, -1)
me -= _mean
me /= _std
print(me)
print(predictor.predict_proba(me)[0][1])

# Predict abue
print("\nPredicting abue:")
abue = array([1.,0.,84.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.]).reshape(1, -1)
abue -= _mean
abue /= _std
print(predictor.predict_proba(abue)[0][1])
# Predict mama
print("\nPredicting mama:")
mama = array([1.,0.,61.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]).reshape(1, -1)
mama -= _mean
mama /= _std
print(predictor.predict_proba(mama)[0][1])
# Predict eva
print("\nPredicting eva:")
eva = array([1.,0.,27.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]).reshape(1, -1)
eva -= _mean
eva /= _std
print(eva)
print(predictor.predict_proba(eva)[0][1])

# Predict bad
print("\nPredicting bad:")
bad = array([1.,1.,50.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]).reshape(1, -1)
bad -= _mean
bad /= _std
print(bad)
print(predictor.predict_proba(bad)[0][1])

print("If probability of death > {}, then counted as will die".format(thresh))
