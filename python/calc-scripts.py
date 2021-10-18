import math as m
def f(x): return x*x-4 # function to be referenced in all below funct's
def derivative(r,n):
    '''1/n = horizontal length of secant segment'''
    return n*(f(r+1/n)-f(r))
def lowfakeintegral(a,b,n):
    '''a,b are endpts, n is number of sections'''
    ans,intervalsize=0,(b-a)/n
    for i in range(n):
        ans+=f(a+i*intervalsize)
    return ans*intervalsize
def highfakeintegral(a,b,n):
    '''a,b are endpts, n is number of sections'''
    ans,intervalsize=0,(b-a)/n
    for i in range(n):
        ans+=f(a+(i+1)*intervalsize)
    return ans*intervalsize
def trapezoid(a,b,n):
    ans,intervalsize=0,(b-a)/n
    for i in range(n):
        ans+=f(a+i*intervalsize)a+f(a+(i+1)*intervalsize)
    return ans*intervalsize/2
def simpson(a,b,n):
    '''n MUST be even for simpson to work!'''
    ans,intervalsize=f(a)-f(b),(b-a)/n
    for i in range(n//2):
        ans+=4*f(a+(2*i+1)*intervalsize) + 2*f(a+(2*i+2)*intervalsize)
    return ans*intervalsize/3
def newton(g,m):
    '''m iterations, any int m geq 0'''
    N,guess=10**15,g
    for i in range(m):
        line= [derivative(guess,N),guess*derivative(guess,N)-f(guess)]
        # tangent line: y=l[0]*x-l[1]
        guess=line[1]/line[0]
    return guess
