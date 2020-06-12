from datetime import datetime
from numpy import array, exp, log
## Day 0 : 2020-02-22 (53)
day_zero = 53

def readDate(strDatetime):
    return datetime.strptime(strDatetime,"%Y-%m-%d")

def dayOfYear(objDatetime):
    return objDatetime.timetuple().tm_yday

def dayOffset(yDay):
    return yDay - day_zero

def today():
    return datetime.today().timetuple().tm_yday

def totalDays():
    return today() - day_zero

def cost(theta, X, y):
    m = len(X)
    alfa = 0.8
    res = logisticFunc(theta, X)
    cost = 1 / m * (- y.transpose()*log(res) -(1-y).transpose()*log(1-res)).sum()
    grad = alfa/m * (X.transpose() * (res - y)).sum(axis=1)
    return [cost, grad]

def logisticFunc(X, theta):
    z = (theta * X).sum(axis = 1)
    res = 1/(1+exp(-z))
    return res

