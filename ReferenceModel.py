from pyomo.environ import *

# Generic linear program with inequalities and equations

model = AbstractModel()

model.m = Param(within=NonNegativeIntegers)
model.n = Param(within=NonNegativeIntegers)
model.p = Param(within=NonNegativeIntegers)



model.J = RangeSet(1, model.n)

model.pi = Var(model.J, within=NonNegativeReals)

model.v = Param(model.J)
model.Ineqs = RangeSet(1, model.m)
model.Eqs = RangeSet(1, model.p)

model.bIneq = Param(model.Ineqs)
model.bEq = Param(model.Eqs)

model.AIneq = Param(model.Ineqs, model.J)
model.AEq = Param(model.Eqs, model.J)

def inequalities_constraint(model, i):
    return sum(model.AIneq[i, j] * model.pi[j] for j in model.J) >= model.bIneq[i]

def equalities_constraint(model, i):
    return sum(model.AEq[i, j] * model.pi[j] for j in model.J) == model.bEq[i]

model.inequality_constraint = Constraint(model.Ineqs, rule=inequalities_constraint)
model.equality_constraint = Constraint(model.Eqs, rule=equalities_constraint)

def obj(model):
    return summation(model.v, model.pi)

model.Obj = Objective(rule=obj, sense=minimize)


