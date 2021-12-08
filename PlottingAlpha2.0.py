import numpy as np
import numpy.linalg as nl
import matplotlib.pyplot as plt
import math as ma
  
fig = plt.figure(figsize=(50,40))
plt.rcParams.update({'font.size': 30})
plt.style.context('light_background')
  
 
with open('x.txt', 'r') as fx:
    x = np.array([float(i) for i in fx.read().split()])
with open('y.txt', 'r') as fy:
    y = np.array([float(i) for i in fy.read().split()])
with open('delx.txt', 'r') as fdelx:
    delx = np.array([float(i) for i in fdelx.read().split()])
with open('dely.txt', 'r') as fdely:
    dely = np.array([float(i) for i in fdely.read().split()])

n = int(x.shape[0])
print('\nNumber of x-y: {}'.format(n))
if y.shape[0] != n:
    print('len(y) != len(x)')
    input('Press Enter to continue...')
if delx.shape[0] != n and delx.shape[0] != 0:
    print('len(delx) != len(x)')
    input('Press Enter to continue...')
if dely.shape[0] != n and dely.shape[0] != 0:
    print('len(dely) != len(x)')
    input('Press Enter to continue...')

St = np.sum((y-np.full(shape=y.shape,fill_value=np.mean(y),dtype=np.float))**2)

typ = int(input('\nWhat is your function type: \n (1) y = k(n).x^n + ... + k(m).x^m \n (2) y = k(n).x^1/n + ... + k(0) + ... + k(m).x^-1/m \n (3) y = A(0) + ... + A(n).cos(n.x) + B(n).sin(n.x) \n (4) Plot only \nEnter your type number: ' ))
nex1 = 2j
ney1 = 2j
nex2 = 2j
ney2 = 2j
cox = -1
if typ in [1]:
    nex2 = complex(input("\nWhat is your n (If you dont know, type 1j): " ))
    nex1 = complex(input("What is your m (If you dont know, type 1j): "))
if typ in [2]:
    ney2 = complex(input("\nWhat is your n (If you dont know, type 1j): " ))
    ney1 = complex(input("What is your m (If you dont know, type 1j): "))
if typ in [3]:
    cox = int(input("\nWhat is your n (If you dont know, type 0): " ))
if nex1.imag == 0:
    nex1 = int(nex1.real)
if nex2.imag == 0:
    nex2 = int(nex2.real)
if ney1.imag == 0:
    ney1 = int(ney1.real)
if ney2.imag == 0:
    ney2 = int(ney2.real)
if typ == 1 and nex2.imag == 0 and nex1.imag == 0 and nex2.real - nex1.real + 2 > n:
    print('Insufficient data')
    nex1 = 2j
    nex2 = 2j
if typ == 2 and ney2.imag == 0 and ney1.imag == 0 and ney2.real - ney1.real + 2 > n:
    print('Insufficient data')
    ney1 = 2j
    ney2 = 2j
if typ == 3 and cox != 0 and 2*cox + 2 > n:
    print('Insufficient data')
    cox = -1
    
if cox == 0:
    ctemp = int(input('What is the range (type 0 if you dont know): '))
        
con = 0
exln = 0
if typ in [1,2]:
    con = int(input('\nDo you want to convert x or y: \n (0) No \n (1) convert x \n (2) convert y \n (3) convert both x and y \nEnter your choice: '))
if con in [1,2,3]:
    exln = int(input('\nWhat do you want to convert to: \n (1) a into exp(a) \n (2) a in to ln(a) \n (3) a into b^a \n (4) a in to logb(a) \nEnter your choice: '))
if exln in [3,4]:
    mun = float(input('\nEnter your b: '))

if con in [1,3]:
    if exln == 1:
        for i in range(n):
            x[i] = ma.exp(x[i])
    if exln == 2:
        for i in range(n):
            x[i] = ma.log(x[i])
    if exln == 3:
        for i in range(n):
            x[i] = mun**(x[i])
    if exln == 4:
        for i in range(n):
            x[i] = ma.log(x[i],mun)
if con in [2,3]:
    if exln == 1:
        for i in range(n):
            y[i] = ma.exp(y[i])
    if exln == 2:
        for i in range(n):
            y[i] = ma.log(y[i])
    if exln == 3:
        for i in range(n):
            y[i] = mun**(y[i])
    if exln == 4:
        for i in range(n):
            y[i] = ma.log(y[i],mun)

if typ == 4:
    plt.scatter(x,y,s=50,c='r',marker='o')
    plt.plot(x,y, 'y')
    if len(delx) != 0:
        plt.errorbar(x,y,xerr=delx,ecolor='b')    
    if len(dely) != 0:
        plt.errorbar(x,y,yerr=dely,ecolor='b')
    plt.savefig('graph.png', bbox_inches='tight')
    plt.close(fig)
if typ in [1,2,3]:
    plt.scatter(x,y,s=50,c='r',marker='o')
    if len(delx) != 0:
        plt.errorbar(x,y,xerr=delx,ecolor='b') 
    if len(dely) != 0:
        plt.errorbar(x,y,yerr=dely,ecolor='b')

    
# Type 1
  
def k_matrix_1(nu,mu):
    global x,y
    y_matrix = np.array([np.sum(y*x**(nu-i)) for i in range(abs(nu+1-mu))])
    
    x_matrix = np.array([[np.sum(x**(nu-i)*x**(nu-j)) for i in range(abs(nu+1-mu))] for j in range(abs(nu+1-mu))])
        
    if nl.det(x_matrix) == 0:
        k_matrix = np.zeros_like(y_matrix)
    else:
        k_matrix = nl.inv(x_matrix).dot(y_matrix)
    return k_matrix

def r2_1(nu,mu,k):
    global x,y,n
    if np.array_equal(k,np.zeros_like(k)):
        Sr = St
    else:
        the_x = np.array([x**(nu-i) for i in range(abs(nu+1-mu))])
        the_y = np.array([the_x[:,i].dot(k) for i in range(n)])
        Sr = np.sum((y-the_y)**2)
    return (St - Sr)/St

def res_1(k,nu,mu,r):
    f = open('result.txt', 'w+')
    print('\nFunction y = k({}).x^{} + ... + k({}).x^{} has: \n'.format(nu,nu,mu,mu))
    f.write('Function y = k({}).x^{} + ... + k({}).x^{} has: \n'.format(nu,nu,mu,mu))
    for i in range(abs(nu+1-mu)):
        print(' k({}) = {}'.format(nu-i,k[i]))
        f.write('\n k({}) = {}'.format(nu-i,k[i]))
    print('\nr^2 = {}'.format(r))
    f.write('\n \nr^2 = {}'.format(r))
    f.close()

def lin_1(nu,mu,k):
    global x
    tim = (np.max(x) - np.min(x)) / 200
    xtem = np.arange(np.min(x)-30*tim,np.max(x)+30*tim,tim)
    the_x = np.array([xtem**(nu-i) for i in range(abs(nu+1-mu))])
    ytem = np.array([the_x[:,i].dot(k) for i in range(xtem.shape[0])])
    plt.plot(xtem,ytem, 'y')
    plt.savefig('graph.png', bbox_inches='tight')

def find_1(nu,mu,k):
    with open('findx.txt' , 'r') as fdx:
        a = np.array([float(i) for i in fdx.read().split()])
    if len(a) != 0:
        the_x = np.array([a**(nu-i) for i in range(abs(nu+1-mu))])
        yr = [the_x[:,i].dot(k) for i in range(the_x.shape[1])]
        fdy = open('y_result.txt', 'w+')
        for j in range(len(yr)):
            fdy.write('\nx = {} => y = {:.6f}'.format(a[j],yr[j]))
            print('\nx = {} => y = {:.6f}'.format(a[j],yr[j]))
        fdy.close()
        plt.scatter(a,yr,c='g')
        

if nex1.imag == 0 and nex2.imag == 0:
    k_matrix = k_matrix_1(nex1,nex2)
    r2 = r2_1(nex1,nex2,k_matrix)       
    res_1(k_matrix,nex1,nex2,r2)
    lin_1(nex1,nex2,k_matrix)
    find_1(nex1,nex2,k_matrix)
                    
if nex1.imag == 0 and nex2.imag == 1:
    k = [k_matrix_1(nex1,nex1-(n-2)+i) for i in range(n-1)]
    r = np.array([r2_1(nex1,nex1-(n-2)+i,k[i]) for i in range(n-1)])
    ind = np.argmax(r)
    nex2 = nex1-n+2+ind
    k_matrix = k[ind]
    r2 = r[ind]
    res_1(k_matrix,nex1,nex2,r2)
    lin_1(nex1,nex2,k_matrix)
    find_1(nex1,nex2,k_matrix)
                   
if nex1.imag == 1 and nex2.imag == 0:
    k = [k_matrix_1(nex2+(n-2)-i,nex2) for i in range(n-1)]
    r = np.array([r2_1(nex2+(n-2)-i,nex2,k[i]) for i in range(n-1)])
    ind = np.argmax(r)
    nex1 = nex2 + n - 2 - ind
    k_matrix = k[ind]
    r2 = r[ind]
    res_1(k_matrix,nex1,nex2,r2)
    lin_1(nex1,nex2,k_matrix)
    find_1(nex1,nex2,k_matrix)
    
if nex1.imag == 1 and nex2.imag == 1:
    k = np.array([[k_matrix_1(i-15,i-15-j) for j in range(n-1)] for i in range(31)])
    r = np.array([[r2_1(i-15,i-15-j,k[i,j]) for j in range(n-1)] for i in range(31)])
    i,j = np.unravel_index(r.argmax(),r.shape)
    nex1 = i-15
    nex2 = i-15-j
    k_matrix = k[i,j]
    r2 = r[i,j]
    res_1(k_matrix,nex1,nex2,r2)
    lin_1(nex1,nex2,k_matrix)
    find_1(nex1,nex2,k_matrix)
                    
                    
# Type 2
                    
def k_matrix_2_1(nu,mu):
    global x,y
    y_matrix = np.array([np.sum(y*x**(1/(nu - i))) for i in range(abs(nu+1-mu))])

    x_matrix = np.array([[np.sum(x**(1/(nu - i))*x**(1/(nu - j))) for j in range(abs(nu + 1 - mu))] for i in range(abs(nu + 1 - mu))])
    
    if nl.det(x_matrix) == 0:
        k_matrix = np.zeros_like(y_matrix)
    else:
        k_matrix = nl.inv(x_matrix).dot(y_matrix)
    return k_matrix

def k_matrix_2_down(mu):
    global x,y,n
    y_matrix = np.array([np.sum(y)] + [np.sum(y*(x**(1/(-1-i)))) for i in range(abs(mu))])
            
    x_matrix = np.array([[n] + [np.sum(x**(1/(- 1 - i))) for i in range(abs(mu))]] \
                        + [[np.sum(x**(1/(- 1 - i)))] + [np.sum(x**(1/(- 1 - i))*x**(1/(- 1 - j))) for j in range(abs(mu))] for i in range(abs(mu))])
        
    if nl.det(x_matrix) == 0:
        k_matrix = np.zeros_like(y_matrix)
    else:
        k_matrix = nl.inv(x_matrix).dot(y_matrix)
    return k_matrix
        
def k_matrix_2_up(nu):
    global x,y,n
    y_matrix = np.array([np.sum(y*x**(1/(nu - i))) for i in range(nu)] + [np.sum(y)])
                   
    x_matrix = np.array([[np.sum(x**(1/(nu - i))*x**(1/(nu - j))) for j in range(nu)] + [np.sum(x**(1/(nu - i)))] for i in range(nu)] \
                        + [[np.sum(x**(1/(nu - i))) for i in range(nu)] + [n]])

    if nl.det(x_matrix) == 0:
        k_matrix = np.zeros_like(y_matrix)
    else:
        k_matrix = nl.inv(x_matrix).dot(y_matrix)
    return k_matrix
    
def k_matrix_2_2(nu,mu):
    global x,y,n
    y_matrix = np.array([np.sum(y*x**(1/(nu - i))) for i in range(nu)] + [np.sum(y)] + [np.sum((1/(- 1 - i))) for i in range(abs(mu))])
            
    x_matrix = np.array([[np.sum(x**(1/(nu - j))*x**(1/(nu - i))) for j in range(nu)] + [np.sum(x**(1/(nu - i)))] + [np.sum(x**(1/(- 1 - j))*x**(1/(nu - i))) for j in range(abs(mu))] for i in range(nu)] \
                        + [[np.sum(x**(1/(nu - i))) for i in range(nu)] + [n] + [np.sum(x**(1/(- 1 - i))) for i in range(abs(mu))]] \
                        + [[np.sum(x**(1/(- 1 - i))*x**(1/(nu - j))) for j in range(nu)] + [np.sum(x**(1/(- 1 - i)))] + [np.sum(x**(1/(- 1 - i))*x**(1/(- 1 - j))) for j in range(abs(mu))] for i in range(abs(mu))])

    if nl.det(x_matrix) == 0:
        k_matrix = np.zeros_like(y_matrix)
    else:
        k_matrix = nl.inv(x_matrix).dot(y_matrix)
    return k_matrix

def r2_2_1(nu,mu,k):
    global x,y,n
    if np.array_equal(k,np.zeros_like(k)):
        Sr = St
    else:
        the_x = np.array([x**(1/(nu - i)) for i in range(abs(nu+1-mu))])
        the_y = np.array([the_x[:,i].dot(k) for i in range(n)])
        Sr = np.sum((y-the_y)**2)
    return (St - Sr)/St

def r2_2_down(mu,k):
    global x,y,n
    if np.array_equal(k,np.zeros_like(k)):
        Sr = St
    else:
        the_x = np.append(np.array([[1 for i in range(n)]]),np.array([x**(1/(- 1 - i)) for i in range(abs(mu))])).reshape((abs(mu)+1,n))
        the_y = np.array([the_x[:,i].dot(k) for i in range(n)])
        Sr = np.sum((y-the_y)**2)
    return (St - Sr)/St

def r2_2_up(nu,k):
    global x,y,n
    if np.array_equal(k,np.zeros_like(k)):
        Sr = St
    else:
        the_x = np.append(np.array([x**(1/(nu - i)) for i in range(nu)]),np.array([[1 for i in range(n)]])).reshape((nu+1,n))
        the_y = np.array([the_x[:,i].dot(k) for i in range(n)])
        Sr = np.sum((y-the_y)**2)
    return (St - Sr)/St

def r2_2_2(nu,mu,k):
    global x,y,n
    if np.array_equal(k,np.zeros_like(k)):
        Sr = St
    else:
        the_x = np.append(np.array([x**(1/(nu - i)) for i in range(nu)]),np.array([[1 for i in range(n)]]))
        the_x = np.append(the_x,np.array([x**(1/(- 1 - i)) for i in range(abs(mu))])).reshape(nu+1-mu,n)
        the_y = np.array([the_x[:,i].dot(k) for i in range(n)])
        Sr = np.sum((y-the_y)**2)
    return (St - Sr)/St

def res_2(k,nu,mu,r):
    f = open('result.txt', "w+")
    print('\nFunction y = k({}).x^1/{} + ... + k({}).x^1/{} has: \n'.format(nu,nu,mu,mu))
    f.write('Function y = k({}).x^1/{} + ... + k({}).x^1/{} has: \n'.format(nu,nu,mu,mu))
    for i in range(abs(nu+1-mu)):
        print(' k({}) = {:.6f}'.format(nu-i,k[i]))
        f.write('\n k({}) = {:.6f}'.format(nu-i,k[i]))
    print('\nr^2 = {:.6f}'.format(r))
    f.write('\n \nr^2 = {:.6f}'.format(r))
    f.close()

def lin_2_1(nu,mu,k):
    global x
    tim = (np.max(x) - np.min(x)) / 200
    xtem = np.arange(np.min(x)-2*tim,np.max(x)+3*tim,tim)
    the_x = np.array([xtem**(1/(nu-i)) for i in range(abs(nu+1-mu))])
    ytem = np.array([the_x[:,i].dot(k) for i in range(xtem.shape[0])])
    plt.plot(xtem,ytem, 'y')
    plt.savefig('graph.png', bbox_inches='tight')
                
def lin_2_down(mu,k):
    global x
    tim = (np.max(x) - np.min(x)) / 200
    xtem = np.arange(np.min(x)-2*tim,np.max(x)+3*tim,tim)
    the_x = np.append(np.array([[1 for i in range(xtem.shape[0])]]),np.array([xtem**(1/(- 1 - i)) for i in range(abs(mu)+1)])).reshape((abs(mu)+1,xtem.shape[0]+1))
    ytem = np.array([the_x[:,i].dot(k) for i in range(xtem.shape[0])])
    plt.plot(xtem,ytem, 'y')
    plt.savefig('graph.png', bbox_inches='tight')
                
def lin_2_up(nu,k):
    global x
    tim = (np.max(x) - np.min(x)) / 200
    xtem = np.arange(np.min(x)-2*tim,np.max(x)+3*tim,tim)
    the_x = np.append(np.array([xtem**(1/(nu - i)) for i in range(nu)]),np.array([[1 for i in range(xtem.shape[0])]])).reshape((nu+1,xtem.shape[0]))
    ytem = np.array([the_x[:,i].dot(k) for i in range(xtem.shape[0])])
    plt.plot(xtem,ytem, 'y')
    plt.savefig('graph.png', bbox_inches='tight')
                
def lin_2_2(nu,mu,k):
    global x
    tim = (np.max(x) - np.min(x)) / 200
    xtem = np.arange(np.min(x)-2*tim,np.max(x)+3*tim,tim)
    the_x = np.append(np.array([xtem**(1/(nu - i)) for i in range(nu)]),np.array([[1 for i in range(xtem.shape[0])]]))
    the_x = np.append(the_x,np.array([xtem**(1/(- 1 - i)) for i in range(abs(mu))])).reshape(nu+1-mu,xtem.shape[0])
    ytem = np.array([the_x[:,i].dot(k) for i in range(xtem.shape[0])])
    plt.plot(xtem,ytem, 'y')
    plt.savefig('graph.png', bbox_inches='tight')

def find_2_1(nu,mu,k):
    with open('findx.txt' , 'r') as fdx:
        a = np.array([float(i) for i in fdx.read().split()])
    if len(a) != 0:
        the_x = np.array([a**1(nu-i) for i in range(abs(nu+1-mu))])
        yr = [the_x[:,i].dot(k) for i in range(the_x.shape[1])]
        fdy = open('y_result.txt', 'w+')
        for j in range(len(yr)):
            fdy.write('\nx = {} => y = {:.6f}'.format(a[j],yr[j]))
            print('\nx = {} => y = {:.6f}'.format(a[j],yr[j]))
        fdy.close()
        plt.scatter(a,yr,c='g')
        
def find_2_down(mu,k):
    with open('findx.txt' , 'r') as fdx:
        a = np.array([float(i) for i in fdx.read().split()])
    if len(a) != 0:
        the_x = np.append(np.array([[1 for i in range(a.shape[0])]]),np.array([a**(1/(- 1 - i)) for i in range(abs(mu))])).reshape((abs(mu)+1,a.shape[0]))
        yr = [the_x[:,i].dot(k) for i in range(the_x.shape[1])]
        fdy = open('y_result.txt', 'w+')
        for j in range(len(yr)):
            fdy.write('\nx = {} => y = {:.6f}'.format(a[j],yr[j]))
            print('\nx = {} => y = {:.6f}'.format(a[j],yr[j]))
        fdy.close()
        plt.scatter(a,yr,c='g')

def find_2_up(nu,k):
    with open('findx.txt' , 'r') as fdx:
        a = np.array([float(i) for i in fdx.read().split()])
    if len(a) != 0:
        the_x = np.append(np.array([a**(1/(nu - i)) for i in range(nu)]),np.array([[1 for i in range(a.shape[0])]])).reshape((nu+1,a.shape[0]))
        yr = [the_x[:,i].dot(k) for i in range(the_x.shape[1])]
        fdy = open('y_result.txt', 'w+')
        for j in range(len(yr)):
            fdy.write('\nx = {} => y = {:.6f}'.format(a[j],yr[j]))
            print('\nx = {} => y = {:.6f}'.format(a[j],yr[j]))
        fdy.close()
        plt.scatter(a,yr,c='g')
        
def find_2_2(nu,mu,k):
    with open('findx.txt' , 'r') as fdx:
        a = np.array([float(i) for i in fdx.read().split()])
    if len(a) != 0:
        the_x = np.append(np.array([a**(1/(nu - i)) for i in range(nu)]),np.array([[1 for i in range(a.shape[0])]]))
        the_x = np.append(the_x,np.array([x**(1/(- 1 - i)) for i in range(abs(mu))])).reshape(nu+1-mu,a.shape[0])
        yr = [the_x[:,i].dot(k) for i in range(the_x.shape[1])]
        fdy = open('y_result.txt', 'w+')
        for j in range(len(yr)):
            fdy.write('\nx = {} => y = {:.6f}'.format(a[j],yr[j]))
            print('\nx = {} => y = {:.6f}'.format(a[j],yr[j]))
        fdy.close()
        plt.scatter(a,yr,c='g')


if ney1.imag == 0 and ney2.imag == 0:
    if ney1 * ney2 > 0:
        k_matrix = k_matrix_2_1(ney1,ney2)   
        r2 = r2_2_1(ney1,ney2,k_matrix)    
        res_2(k_matrix,ney1,ney2,r2)    
        lin_2_1(ney1,ney2,k_matrix)
        find_2_1(ney1,ney2,k_matrix)
    if ney1 == 0:
        k_matrix = k_matrix_2_down(ney2)        
        r2 = r2_2_down(ney2,k_matrix)        
        res_2(k_matrix,0,ney2,r2)   
        lin_2_down(ney2,k_matrix)
        find_2_down(ney2,k_matrix)
    if ney2 == 0:        
        k_matrix = k_matrix_2_up(ney1)        
        r2 = r2_2_up(ney1,k_matrix)        
        res_2(k_matrix,ney1,0,r2)    
        lin_2_up(ney1,k_matrix)
        find_2_up(ney1,k_matrix)            
    if ney1.real * ney2.real < 0:        
        k_matrix = k_matrix_2_2(ney1,ney2)        
        r2 = r2_2_2(ney1,ney2,k_matrix)        
        res_2(k_matrix,ney1,ney2,r2)        
        lin_2_2(ney1,ney2,k_matrix)                    
        find_2_2(ney1,ney2,k_matrix)           
        
if ney1.imag == 0 and ney2.imag == 1:   
    if ney1 - n + 2 > 0 or ney1 < 0:                   
        k = [k_matrix_2_1(ney1,ney1-(n-2)+i) for i in range(n-1)]
        r = np.array([r2_2_1(ney1,ney1-(n-2)+i,k[i]) for i in range(n-1)])
        ind = np.argmax(r)  
        ney2 = ney1 - n + 2 + ind
        k_matrix = k[ind]    
        r2 = r[ind]
        res_2(k_matrix,ney1,ney2,r2)
        lin_2_1(ney1,ney2,k_matrix)
        find_2_1(ney1,ney2,k_matrix)
    else:
        k = [k_matrix_2_1(ney1,ney1-i) for i in range(ney1)] + [k_matrix_2_up(ney1)] + [k_matrix_2_2(ney1,-1-i) for i in range(abs(ney1-(n-2)))]
        r = np.array([r2_2_1(ney1,ney1-i,k[i]) for i in range(ney1)] + [r2_2_up(ney1,k[ney1])] + [r2_2_2(ney1,-1-i,k[ney1+1+i]) for i in range(abs(ney1-(n-2)))])
        ind = np.argmax(r)
        ney2 = ney1 - ind
        k_matrix = k[ind]
        r2 = r[ind]
        res_2(k_matrix,ney1,ney2,r2)
        if ney2 > 0:
            lin_2_1(ney1,ney2,k_matrix)
            find_2_1(ney1,ney2,k_matrix)
        if ney2 == 0:
            lin_2_up(ney1,k_matrix)
            find_2_up(ney1,k_matrix)
        if ney2 < 0:
            lin_2_2(ney1,ney2,k_matrix)
            find_2_2(ney1,ney2,k_matrix)
                                                   
if ney1.imag == 1 and ney2.imag == 0: 
    if ney2 + n - 2 < 0 or ney2 > 0:    
        k = [k_matrix_2_1(ney2+i,ney2) for i in range(n-1)]
        r = np.array([r2_2_1(ney2+i,ney2,k[i]) for i in range(n-1)])
        ind = np.argmax(r)
        ney1 = ney2 + ind  
        k_matrix = k[ind]    
        r2 = r[ind]
        res_2(k_matrix,ney1,ney2,r2)    
        lin_2_1(ney1,ney2,k_matrix)
        find_2_1(ney1,ney2,k_matrix)  
    else:      
        k = [k_matrix_2_2(ney2+(n-2)-i,ney2) for i in range(ney2+(n-2))] + [k_matrix_2_down(ney2)] + [k_matrix_2_1(-1-i,ney2) for i in range(abs(ney2))]
        r = np.array([r2_2_2(ney2+(n-2)-i,ney2,k[i]) for i in range(ney2+(n-2))] + [r2_2_down(ney2,k[ney2+(n-2)])] + [r2_2_1(-1-i,ney2,k[ney2+(n-2)+1+i]) for i in range(abs(ney2))])
        ind = np.argmax(r)
        ney1 = ney2 + (n-2) - ind
        k_matrix = k[ind]
        r2 = r[ind]
        res_2(k_matrix,ney1,ney2,r2)
        if ney1 > 0:
            lin_2_2(ney1,ney2,k_matrix)
            find_2_2(ney1,ney2,k_matrix)
        if ney1 == 0:
            lin_2_down(ney2,k_matrix)
            find_2_down(ney2,k_matrix)
        if ney1 < 0:
            lin_2_1(ney1,ney2,k_matrix)
            find_2_1(ney1,ney2,k_matrix)
             
if ney1.imag == 1 and ney2.imag == 1:
    k = np.array([[k_matrix_2_1(i-15,i-15-j) for j in range(n-1)] for i in range(15)] +\
                   [[np.zeros(n+1)]+[k_matrix_2_down(1-j) for j in range(n-2)]] +\
                   [[k_matrix_2_1(i+1,i+1-j) for j in range(i+1)] + [k_matrix_2_up(i+1)] + [k_matrix_2_2(i+1,-1-j) for j in range(abs(i+1-(n-2)))] for i in range(n-2)] +\
                   [[k_matrix_2_1(n-1+i,n-1+i-j) for j in range(n-1)] for i in range(15-(n-2))]).reshape((31,n-1))
    r = np.array([[r2_2_1(i-15,i-15-j,k[i,j]) for j in range(n-1)] for i in range(15)]+\
                   [[0]+[r2_2_down(1-j,k[15,1+j]) for j in range(n-2)]] +\
                   [[r2_2_1(i+1,i+1-j,k[16+i,j]) for j in range(i+1)] + [r2_2_up(i+1,k[16+i,i+1])] + [r2_2_2(i+1,-1-j,k[16+i,i+2+j]) for j in range(abs(i+1-(n-2)))] for i in range(n-2)] +\
                   [[r2_2_1(n-1+i,n-1+i-j,k[16+n-2+i,j]) for j in range(n-1)] for i in range(15-(n-2))]).reshape((31,n-1))
    i,j = np.unravel_index(r.argmax(),r.shape)
    ney1 = i-15
    ney2 = ney1 - j
    k_matrix = k[i,j]    
    r2 = r[i,j]
    res_2(k_matrix,ney1,ney2,r2)
    if ney1 * ney2 > 0:
        lin_2_1(ney1,ney2,k_matrix)
        find_2_1(ney1,ney2,k_matrix)
    if ney1 == 0:
        lin_2_down(ney2,k_matrix)
        find_2_down(ney2,k_matrix)
    if ney2 == 0:
        lin_2_up(ney1,k_matrix)
        find_2_up(ney1,k_matrix)
    if ney1 * ney2 < 0:
        lin_2_2(ney1,ney2,k_matrix)
        find_2_2(ney1,ney2,k_matrix)
            
            
# Type 3
    
def k_matrix_3(c):
    global x,y,n
    y_matrix = np.append(np.array([np.sum(y)]),np.array([[np.sum(y*np.cos((i+1)*x))] + [np.sum(y*np.sin((i+1)*x))] for i in range(c)]))
        
    x_matrix = np.append(np.append(np.array([n]),np.array([[np.sum(np.cos((i+1)*x))] + [np.sum(np.sin((i+1)*x))] for i in range(c)]).reshape(-1)),
                         np.stack(tuple([[np.append(np.sum(np.cos((j+1)*x)),np.array([[np.sum(np.cos((i+1)*x)*np.cos((j+1)*x))] + [np.sum(np.sin((i+1)*x)*np.cos((j+1)*x))] for i in range(c)]).reshape(-1))] +\
                                         [np.append(np.sum(np.sin((j+1)*x)),np.array([[np.sum(np.cos((i+1)*x)*np.sin((j+1)*x))] + [np.sum(np.sin((i+1)*x)*np.sin((j+1)*x))] for i in range(c)]).reshape(-1))] for j in range(c)])).reshape(-1)).reshape((2*c+1,2*c+1))
    print(x_matrix.shape)
        
    if nl.det(x_matrix) == 0:
        k_matrix = np.zeros_like(y_matrix)
    else:
        k_matrix = nl.inv(x_matrix).dot(y_matrix)
    return k_matrix
    
def r2_3(c,k):
    global x,y,n
    if np.array_equal(k,np.zeros_like(k)):
        Sr = St
    else:
        the_x = np.append(np.array([1 for i in range(n)]),np.stack(tuple([[np.cos((i+1)*x)] + [np.sin((i+1)*x)] for i in range(c)])).reshape((c*2,n))).reshape((c*2+1,n))
        the_y = np.array([the_x[:,i].dot(k) for i in range(n)])
        Sr = np.sum((the_y-y)**2)
    return (St - Sr)/St
 
def res_3(k,c,r):
    f = open('result.txt', "w+")
    print('\nFunction y = A(0) + ... + A({}).cos({}.x) + B({}).sin({}.x) has: \n'.format(c,c,c,c))
    f.write('Function y = A(0) + ... + A({}).cos({}.x) + B({}).sin({}.x) has: \n'.format(c,c,c,c))
    print(' A(0) = {:.6f}'.format(k[0]))
    f.write('\n A(0) = {:.6f}'.format(k[0]))
    for i in range(cox):
        print(' A({}) = {:.6f}'.format(i+1,k[2*i + 1]))
        f.write('\n A({}) = {:.6f}'.format(i+1,k[2*i + 1]))
        print(' B({}) = {:.6f}'.format(i+1,k[2*i + 2]))
        f.write('\n B({}) = {:.6f}'.format(i+1,k[2*i + 2]))
    print('\nr^2 = {:.6f}'.format(r))
    f.write('\n \nr^2 = {:.6f}'.format(r))
    f.close()
        
def lin_3(c,k):
    global x
    tim = (np.max(x) - np.min(x)) / 400
    xtem = np.arange(np.min(x)-4*tim,np.max(x)+5*tim,tim)
    the_x = np.append(np.array([1 for i in range(xtem.shape[0])]),np.stack(tuple([[np.cos((i+1)*xtem)] + [np.sin((i+1)*xtem)] for i in range(c)])).reshape((c*2,xtem.shape[0]))).reshape((c*2+1,xtem.shape[0]))
    ytem = np.array([the_x[:,i].dot(k) for i in range(xtem.shape[0])])
    plt.plot(xtem,ytem, 'y')
    plt.savefig('graph.png', bbox_inches='tight')
    
def find_3(c,k):
    with open('findx.txt' , 'r') as fdx:
        a = np.array([float(i) for i in fdx.read().split()])
    if len(a) != 0:
        the_x = np.append(np.array([1 for i in range(a.shape[0])]),np.stack(tuple([[np.cos((i+1)*a)] + [np.sin((i+1)*a)] for i in range(c)])).reshape((c*2,a.shape[0]))).reshape((c*2+1,a.shape[0]))
        yr = [the_x[:,i].dot(k) for i in range(the_x.shape[1])]
        fdy = open('y_result.txt', 'w+')
        for j in range(len(yr)):
            fdy.write('\nx = {} => y = {:.6f}'.format(a[j],yr[j]))
            print('\nx = {} => y = {:.6f}'.format(a[j],yr[j]))
        fdy.close()
        plt.scatter(a,yr,c='g')
                
                
if cox > 0:
    k_matrix = k_matrix_3(cox)
    r2 = r2_3(cox,k_matrix)
    res_3(k_matrix,cox,r2)
    lin_3(cox,k_matrix)
    find_3(cox,k_matrix)
    
if cox == 0:
    if ctemp != 0:
        k = [k_matrix_3(1+i) for i in range(ctemp)]
        r = np.array([r2_3(1+i,k[i]) for i in range(ctemp)])
    else:
        k = [k_matrix_3(1+i) for i in range(int((n-2)/2))]
        r = np.array([r2_3(1+i,k[i]) for i in range(int((n-2)/2))])     
    ind = np.argmax(r)
    cox = ind + 1
    k_matrix = k[ind]    
    r2 = r[ind] 
    res_3(k_matrix,cox,r2)
    lin_3(cox,k_matrix) 
    find_3(cox,k_matrix)
    
plt.close(fig)