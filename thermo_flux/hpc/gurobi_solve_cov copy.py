from gurobipy import *

file_in = 'double_reg'

m = read(file_in+'.mps')

log_name = file_in+_.txt'
sol_name = file_in+'_.sol'
sol_write = file_in

m.params.LogFile = log_name

m.params.OutputFlag = 1
m.params.TimeLimit = 60*60*16
m.params.NonConvex = 2
m.params.Threads=8
m.params.Seed = seed
m.params.SolFiles = sol_write
m.params.IntegralityFocus = 1


m.optimize()

m.write(sol_name)


