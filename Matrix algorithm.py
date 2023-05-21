import numpy as np
from pyomo.environ import ConcreteModel, Var, Objective, SolverFactory, ConstraintList, RangeSet, Constraint
from pyomo.environ import PositiveIntegers,Binary, Integers

# a=np.array([[0,2],
#             [1,1]])

# b=np.array([[-5,2],
#             [1,6]])

# c=np.dot(a,b)
# print(c)


def basic_multiplication_alogrithm(n):
    a=np.zeros((n,n))
    c=np.zeros((n**2,n**2,n**2))
    indice=0
    for i in range(n):
        for j in range(n):
            indice+=1
            a[i,j]=indice
    b=np.copy(a)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                indice_a=a[i,k]
                indice_b=b[k,j]
                c[int(a[i,j])-1, int(indice_a)-1, int(indice_b)-1] = 1
    return c


#c=basic_multiplication_alogrithm(2)
#print(c)

def resolution(n):

    model = ConcreteModel()

    K=n**3
    model.I=RangeSet(1,n**2)
    model.M=RangeSet(1,K)
    model.u=Var(model.I, model.M,within=Binary)
    model.v=Var(model.I, model.M,within=Binary)
    model.w=Var(model.I, model.M,within=Binary)
    model.ab=Var(model.I,model.I,model.M,within=Integers)
    #model.c=Var(model.I,model.I,model.I,within=Binary)
    c=basic_multiplication_alogrithm(n)
    """ #if i convert c in a variable it will get modified to full 0
    for i in range(n**2):
        for j in range(n**2):
            for k in range(n**2):
                model.c[j+1,k+1,i+1] = c[i,j,k] #pyomo : line, column, 3D ; numpy : 3D, line, column
    #display
    for k in range(n**2):
        for i in range(n**2):
            l=[]
            for j in range(n**2):
                l.append(c[i+1,j+1,k+1]())
            print(l) """
    
    #we want to reduce the number of additions for a given K ==> reduce non-zero terms
    #absolute value for -1, 0, 1 : can be cell squared
    #for binary
    model.obj = Objective(expr=sum(sum(model.u[i,m] for m in model.M) for i in model.I) + sum(sum(model.v[i,m] for m in model.M) for i in model.I) + sum(sum(model.w[i,m] for m in model.M) for i in model.I))
    #for -1,0,1 :
    #model.obj = Objective(expr=sum(sum(model.u[i,m]*model.u[i,m] for m in model.M) for i in model.I) + sum(sum(model.v[i,m]*model.v[i,m] for m in model.M) for i in model.I) + sum(sum(model.w[i,m]*model.w[i,m] for m in model.M) for i in model.I))
    
    model.sum=ConstraintList()
    for h in model.I:
        for i in model.I:
            for j in model.I:
                #model.sum.add(expr=sum(model.w[h,k] * model.u[i,k] * model.v[j,k] for k in model.M) == model.c[i,j,h])
                model.sum.add(expr=sum(model.w[h,k] * model.ab[i,j,k] for k in model.M) == c[h-1,i-1,j-1])

    model.product=ConstraintList()
    for i in model.I:
        for j in model.I:
            for k in model.M:
                model.product.add(expr=model.ab[i,j,k] == model.u[i,k] * model.v[j,k])


    # Write the LP model in standard format
    #model.write("matrix_{}.lp".format(n))

    # Solve the model
    sol = SolverFactory('gurobi').solve(model, tee=True)
    #sol = SolverFactory('gurobi').solve(model, tee=True, options_string="NonConvex=2")

    return model.u, model.v, model.w

n=2
u,v,w=resolution(n)
for i in range(1,5):
    print([int(u[i,m]()) for m in range(1,9)])
print("\n")
for i in range(1,5):
    print([int(v[i,m]()) for m in range(1,9)])
print("\n")
for i in range(1,5):
    print([int(w[i,m]()) for m in range(1,9)])
print("\n")
