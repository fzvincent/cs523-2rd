import numpy as np
import math

mydata=np.load("mydata.npz")
mutationRateList=mydata["mutationRateList"]
eliteList=mydata["eliteList"]
islandsList=mydata["islandsList"]
a=mydata["a"]
b=mydata["b"]
c=mydata["c"]

mTime=10
eliteTimes=6


mutationRateList=np.logspace(math.log10(0.01),math.log10(0.90),mTime,endpoint=True)
eliteList=np.logspace(math.log10(0.10),math.log10(0.90),eliteTimes,endpoint=True)