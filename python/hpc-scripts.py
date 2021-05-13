#the modules, do not modify
import math
import cmath
import copy
# im too lazy to solve a quadratic then round the sols, so this is the result
def quad_formula(a,b,c):
    delta = b*b - 4 * a * c
    if (delta >= 0): return {(-b + math.sqrt(delta)) / (2 * a), (-b - math.sqrt(delta)) / (2 * a)}
    return {(-b + math.sqrt(-delta)*1j) / (2 * a), (-b - math.sqrt(-delta)*1j) / (2 * a)}
# the repeated problem of tension in which a weight is suspended by two strings/ropes/whatever
def tension(weight, Langle,Rangle):
    '''angles in degrees'''
    x=weight/(math.tan(Langle)+math.tan(Rangle))
    return [x/math.cos(Langle),x/math.cos(Rangle)]
#triangle completion, also assigned ad nauseam, thank god these exist
#the weird '^\circ' exists to allow for direct copypasting into LaTeX compilers
def trianglesolveAAS(A,B,a):
    C=math.pi-A-B
    b= a*math.sin(B)/math.sin(A)
    c= a*math.sin(C)/math.sin(A)
    return [b,c,str(180*C/math.pi)+'^{\circ}']
def trianglesolveSSS(a,b,c):
    A=math.acos((b**2+c**2-a**2)/(2*b*c))
    B=math.acos((a**2+c**2-b**2)/(2*a*c))
    C=math.acos((a**2+b**2-c**2)/(2*a*b))
    return [str(A*180/math.pi)+'^{\circ}',str(B*180/math.pi)+'^{\circ}',str(C*180/math.pi)+'^{\circ}']
def trianglesolveSAS(a,b,C):
    c=math.sqrt(a**2+b**2-2*a*b*math.cos(C))
    L=trianglesolveSSS(a,b,c)
    return [L[0],L[1],c]
def trianglesolveSSA(A,a,b):
    B1=math.asin(b*math.sin(A)/a)
    B2=math.pi-B1
    C1=B2-A
    C2=B1-A
    c1=trianglesolveAAS(A,B1,a)[1]
    c2=trianglesolveAAS(A,B2,a)[1]
    return [[c1,str(B1*180/math.pi)+'^{\circ}',str(C1*180/math.pi)+'^{\circ}'],[c2,str(B2*180/math.pi)+'^{\circ}',str(C2*180/math.pi)+'^{\circ}']]
def heron(a,b,c): # (duh)
    s= (a+b+c)/2
    return math.sqrt(s*(s-a)*(s-b)*(s-c))
#
# how common core teaches very very basic linear alg:
# scripts for matrix operations, theyre too tedious to do manually more than once, no numpy bc im lazy
def minor(a:int,b:int,M:list):
    l=len(M)
    T=[]
    for i in range(l):
        T.append(tuple(M[i]))
    T=tuple(T)
    L=[]
    for j in range(l):
        L.append(list(T[j]))
    L.pop(a)
    for k in range(l-1):
        L[k].pop(b)
    return L
def det(M:list):
    total=0
    l=len(M)
    if l==1:
        return M[0][0]
    minor_orig=[]
    for i in range(l):
        minor_orig.append(M[i][1:])
    minor_orig=tuple(minor_orig)
    for i in range(l):
        minor=list(minor_orig)
        minor.pop(i)
        total+=((-1)**(i%2))*M[i][0]*det(minor)
    return total
def inverse(M):
    l=len(M)
    T=[]
    for i in range(l):
        T.append(tuple(M[i]))
    T=tuple(T)
    L=[]
    for m in range(l):
        L.append(list(T[m]))
    for p in range(l):
        for q in range(l):
            L[p][q]=det(minor(q,p,M))*(-1)**((p+q)%2)/det(M)
    return L
def multiply(M1:list,M2:list):
    rows1=len(M1)
    cols1=len(M1[0])
    cols2=len(M2[0])
    product=[]
    for s in range(rows1):
        element=[]
        for t in range(cols2):
            s1=0
            for u in range(cols1):
                s1+=M1[s][u]*M2[u][t]
            element.append(s1)
        product.append(element)
    return product
def equationsolve(M:list):
    l=len(M)
    rhs_vec=[]
    for h in range(l):
        rhs_vec.append([M[h][-1]])
    print(f"{rhs_vec=}")
    lhs_matrix=[]
    for g in range(l):
        lhs_matrix.append(M[g][:-1])
    print(f"{lhs_matrix=}")
    if det(lhs_matrix)==0:
        return 'No unique solution'
    inv=inverse(lhs_matrix)
    print(f"{inv=}")
    print(f"{det(lhs_matrix)=}")
    return multiply(inv,rhs_vec)
def scale(r,M:list):
    length=len(M[0])
    height=len(M)
    M1=copy.deepcopy(M)
    for i in range(height):
        for j in range(length):
            M1[i][j]*=r
    return M1
def add(M1,M2):
    ans=copy.deepcopy(M1)
    length1=len(M1[0])
    height1=len(M1)
    for i in range(height1):
        for j in range(length1):
            ans[i][j]+=M2[i][j]
    return ans
def shoelace(vertices):
    '''vertices contains ordered pairs(list)'''
    l=len(vertices)
    middle=copy.deepcopy(vertices)
    middle.append(middle[0])
    left=0
    right=0
    for i in range(l):
        left+=middle[i][1]*middle[i+1][0]
        right+=middle[i][0]*middle[i+1][1]
    return (left-right)/2
# harmonic numbers because they appeared multiple times here and there
def harmonic1(n):
    '''naive'''
    ans=0
    for i in range(n):
        ans+=(1/(i+1))
    return ans
def harmonic2(n):
    gamma=0.577215664901532 #euler-mascheroni const
    return math.log(n)+gamma+1/(2*n)-1/(12*n*n)+1/(120*(n**4)) #first few terms of series but alr rlly accurate
#
#self-explanatory, for use on the stupid riemann sums
def foo(x): 
    return x*x*x/2  #change this to your discretion, to be referenced in 'fakeintegral' functions
def lowfakeintegral(a,b,n):
    '''a,b are endpts, n is number of sections'''
    ans,intervalsize=0,(b-a)/n
    for i in range(n):
        ans+=foo(a+i*intervalsize)
    return ans*intervalsize
def highfakeintegral(a,b,n):
    '''a,b are endpts, n is number of sections'''
    ans,intervalsize=0,(b-a)/n
    for i in range(n):
        ans+=foo(a+(i+1)*intervalsize)
    return ans*intervalsize
highfakeintegral(0,2,4)
#more will be added as more repetitive tasks appear, but hopefully i saved you some time lel
