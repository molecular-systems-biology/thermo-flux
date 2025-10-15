import sys
from gurobipy import *

print(sys.argv)

seed_in = str(sys.argv[1])
seed = int(seed_in)

file_in = 'double_reg'

m = read(file_in+'.mps')

log_name = file_in+seed_in+'_.txt'
sol_name = file_in+seed_in+'_.sol'
mps_name = file_in+seed_in+'_out_.mps'
sol_write = file_in+seed_in

m.params.LogFile = log_name

m.params.OutputFlag = 1
m.params.TimeLimit = 60*60*16
m.params.NonConvex = 2
m.params.Threads=4
m.params.Seed = seed
m.params.SolFiles = sol_write
m.params.IntegralityFocus = 1


def mycallback(model, where):
    if where == GRB.Callback.MIPSOL:
        print('solution found', model.cbGet(GRB.Callback.MIPSOL_OBJ))
        with open('res_'+file_in+".txt", "a") as myfile:
            time = model.cbGet(GRB.Callback.RUNTIME)
            objval = model.cbGet(GRB.Callback.MIPSOL_OBJ)
            myfile.write('\n Seed_'+seed_in+'_'+str(time)+'_'+str(objval))

m.optimize(mycallback)

m.write(sol_name)


