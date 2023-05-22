import numpy as np
from pyomo.environ import ConcreteModel, Var, Objective, SolverFactory, ConstraintList, RangeSet, Constraint
from pyomo.environ import PositiveIntegers,Binary, Integers, NonNegativeIntegers

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
    #the first dimension of c is the number of the cell : 1,2,...n**2
    for i in range(n):
        for j in range(n):
            for k in range(n):
                indice_a=a[i,k]
                indice_b=b[k,j]
                c[int(a[i,j])-1, int(indice_a)-1, int(indice_b)-1] = 1
    return c


#c=basic_multiplication_alogrithm(2)
#print(c)

def resolution(n,k):

    model = ConcreteModel()

    model.I=RangeSet(1,n**2)
    model.M=RangeSet(1,k)
    model.u=Var(model.I, model.M,within={-1,0,1})
    model.v=Var(model.I, model.M,within={-1,0,1})
    model.w=Var(model.I, model.M,within={-1,0,1})
    model.ab=Var(model.I,model.I,model.M,within={-1,0,1})
    #model.m=Var(model.M,within=NonNegativeIntegers)
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
    #1-norm for binary
    #model.obj = Objective(expr=sum(sum(model.u[i,m] for m in model.M) for i in model.I) + sum(sum(model.v[i,m] for m in model.M) for i in model.I) + sum(sum(model.w[i,m] for m in model.M) for i in model.I))
    #1-norm for -1,0,1 :
    model.obj = Objective(expr=sum(sum(model.u[i,m]*model.u[i,m] for m in model.M) for i in model.I) + sum(sum(model.v[i,m]*model.v[i,m] for m in model.M) for i in model.I) + sum(sum(model.w[i,m]*model.w[i,m] for m in model.M) for i in model.I))
    
    #new formula must be equivalent to the basic formula with n**3 multiplication
    model.sum=ConstraintList()
    for h in model.I:
        for i in model.I:
            for j in model.I:
                #model.sum.add(expr=sum(model.w[h,k] * model.u[i,k] * model.v[j,k] for k in model.M) == model.c[i,j,h]) #body with non linear terms
                model.sum.add(expr=sum(model.w[h,k] * model.ab[i,j,k] for k in model.M) == c[h-1,i-1,j-1])

    #in order to reduce the above constraint in a quadratic form, we describe here the multiplication of each aibj for each k
    model.product=ConstraintList()
    for i in model.I:
        for j in model.I:
            for k in model.M:
                model.product.add(expr=model.ab[i,j,k] == model.u[i,k] * model.v[j,k])


    #this represents the number of multiplication : each non zero cell counts for one multiplication
    # model.multipl=ConstraintList()
    # for k in model.M:
    #     model.multpli.add(expr=model.m[k] == sum(w[i,k]*w[i,k] for i in model.I))

    # Write the LP model in standard format
    #model.write("matrix_{}.lp".format(n))

    # Solve the model
    #sol = SolverFactory('gurobi').solve(model, tee=True)
    sol = SolverFactory('gurobi').solve(model, tee=True, options={"NonConvex":2,"TimeLimit":45*60})
    #the solver found a feasible solution if there is an "incumbent" value

    return model.u, model.v, model.w

def algo_verification(u,v,w,n):
    #verifies if the algorithm works
    for test in range(10):
        a=np.random.randint(-10,10,(n,n))
        b=np.random.randint(-10,10,(n,n))
        ab=np.dot(a,b)
        a_line=[]
        b_line=[]
        ab_line=[]
        for i in range(n):
            for j in range(n):
                a_line.append(a[i,j])
                b_line.append(b[i,j])
                ab_line.append(ab[i,j])
        for indice in range(n):
            if ab_line[indice] != sum(w[indice+1,k]() * sum(u[i+1,k]() * a_line[i] for i in range(n**2)) * sum(v[i+1,k]() * b_line[i] for i in range(n**2)) for k in range(1,n**3+1)):
                return False
    return True



n=2
k=n**3-1
u,v,w=resolution(n,k)
I=RangeSet(1,n**2)
M=RangeSet(1,k)
for i in I:
    print([int(u[i,m]()) for m in M])
print("\n")
for i in I:
    print([int(v[i,m]()) for m in M])
print("\n")
for i in I:
    print([int(w[i,m]()) for m in M])
print("\n")

if algo_verification(u,v,w,n):
    print("Functional")
else:
    print("Error, not functional")
