# scripts used during ap bc >.<
import math as m
def f1(x): return x*x-4 #change this to your discretion, to be referenced in 'fakeintegral' functions;
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
    '''m iterations'''
    N=10**15
    guess=g
    for i in range(m):
        line= [derivative(guess,N),guess*derivative(guess,N)-f(guess)] # tangent line: y=l[0]*x-l[1]
        guess=line[1]/line[0]
    return guess
def f2(x,y): return m.e**(x*y)
def euler(x,y,n,step):
    x0,y0=x,y
    for i in range(n):
        y0+=f(x0,y0)*step
        x0+=step
    return y0
# more scripts coming soon >:)
