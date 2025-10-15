from gurobipy import *

m = read('plant_07f_.mps')

m.update()

m.params.LogFile = 'tune_log.txt'

m.params.TuneTrials = 6
m.params.OutputFlag = 1
m.params.TimeLimit = 300
m.params.Threads=16
m.params.IntegralityFocus = 1
m.params.TuneTimeLimit =60*60*24 #24 hours of tuning 


m.tune()

for i in range(m.tuneResultCount):
	m.getTuneResult(i)
	m.write('tune'+str(i)+'.prm')