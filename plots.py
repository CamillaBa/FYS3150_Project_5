import matplotlib.pyplot as plt
from   mpl_toolkits.mplot3d import Axes3D
from   matplotlib import cm
from   matplotlib.ticker import LinearLocator, FormatStrFormatter
from   matplotlib.colors import BoundaryNorm
from   matplotlib.ticker import MaxNLocator
import numpy as np

#=========================================================================
# Problem c)
#=========================================================================

def relerror(target, approximation):
    return np.linalg.norm(target-approximation)/np.linalg.norm(target)


def color_style_analytical(t):
    if t == "t1":
        return "lightpink"
    if t == "t2":
        return "deepskyblue"

def color_style(t):
    if t == "t1":
        return "darkred"
    if t == "t2":
        return "blue"

N = 100 # terms of analytical solution
A = lambda n: 2/(np.pi*n)*(-1)**n
An = np.array([A(n) for n in range(1,N+1)])

def analytical(x,t):
    """
    Takes an np.array x shaped as a vector together with a time t,
    and returns the analytical solution as a vector.
    """
    pi_x           = np.pi*x
    mat_1          = np.array([np.sin(n*pi_x) for n in range(1,N+1)])
    exponentials   = np.array([np.exp(-n**2*np.pi**2*t) for n in range(1,N+1)])
    mat_2          = An*exponentials
    multiplication = np.matmul(mat_2,mat_1)
    return  multiplication + x

timelist = {}
timelist["t1"]=0.05
timelist["t2"]=1

for dx in ["0.100000","0.010000"]:
    x = np.arange(0,1+float(dx),float(dx))
    for method in ["explicit_euler","implicit_euler","crank_nicolson"]:
        plt.figure(method+"_dx_"+dx)
        plt.title(method+"\n"+r" $dx=$"+dx+r" $dt=$"+"0.01"+r"$dx^2$", fontsize=14)
        for t in ["t1","t2"]:
            u = np.loadtxt(method+"_dx_"+dx+"_"+t+".txt", 
                           delimiter=',', 
                           unpack=True)
            plt.plot(x,u,label = "t="+str(timelist[t])+", rel err: " + str(relerror(analytical(x,timelist[t]),u)),
                     color = color_style(t))
            plt.plot(x,analytical(x,timelist[t]), 
                     color = color_style_analytical(t), 
                     label = "t="+str(timelist[t]) + "(analytical)",
                     linestyle = '--')
        plt.xlabel(r"$x$", fontsize=14)
        plt.ylabel(r"$u$", fontsize=14)
        plt.legend(loc="best")
        plt.show(block=False)

#=========================================================================
# Problem f)
#=========================================================================

N = 20 # terms of analytical solution
A = lambda m, n: 4/(np.pi*np.pi*m*n)*((-1)**m-(-1)**(m+n))

def analytical_2D(x,y,t):
    """ Takes as input three floats x, y and t,
    and returns the analytical value for u.
    """
    pi_x = np.pi*x
    pi_y = np.pi*y
    arr  = np.array([A(m,n)*np.sin(m*pi_x)*np.sin(n*pi_y)*np.exp(-(m**2+n**2)*np.pi**2*t) for m in range(1,N+1) for n in range(1,N+1) if m+n < N+1]) 
    return  np.sum(arr) + x


for h in ["0.100000", "0.010000"]:
    for t in ["t1","t2"]:
        filename = "two_dim_implicit_euler_h_"+h+"_"+t
        u = np.loadtxt(filename+".txt",
                       delimiter=',', 
                       unpack=True)
        m = np.size(u[:,0])
        n = np.size(u[0,:])
        
        u_analytical = np.zeros((m,n));
        for i, x in enumerate(np.linspace(0,1,m)):
            for j, y in enumerate(np.linspace(0,1,n)):
                u_analytical[j][i] = analytical_2D(x,y,timelist[t])

        x    = np.arange(0, 1 +float(h), float(h))
        y    = np.arange(0, 1 +float(h), float(h))

        plt.figure(filename)
        plt.title(r"$t=$"+str(timelist[t])+r", $h=$"+h+r" $dt=$"+"0.01"+r"$h^2$"+"\n"
                     + "rel err: " + str(relerror(u_analytical ,u))
                     , fontsize=14)
        plt.plot(x, u[int(0.5/float(h)),:], label = r"$y = 0.5$")
        plt.plot(x, u[int(0.1/float(h)),:], label = r"$y \in \{0.1, 0.9\}$")
        plt.plot(x,  u_analytical[0][:], label = "steady state (analytical)")
        plt.xlabel(r"$x$", fontsize=14)
        plt.ylabel(r"$u$", fontsize=14)
        plt.legend(loc="best")
        plt.show(block = False)

        plt.figure(filename+"_color")
        plt.pcolor(x, y, u, cmap='RdBu_r', vmin=np.min(u), vmax=np.max(u))
        plt.title(r"$t=$"+str(timelist[t])+r", $h=$"+h+r" $dt=$"+"0.01"+r"$h^2$"+"\n"
                     + "rel err: " + str(relerror(u_analytical ,u))
                     , fontsize=14)
        plt.colorbar().set_label(label = r"$u$",size = 14)
        plt.xlabel(r"$x$", fontsize=14)
        plt.ylabel(r"$y$", fontsize=14)
        plt.show(block=False)

#=========================================================================
# Problem g)
#=========================================================================

Q1 = 1.4;       Q2 = 0.35;     Q3 = 0.05
c1 = -71/180;   c2 = -1/1125
d1 = -229/900;  d2 = -23/2250
f1 = -157/900;  f2 = -47/2250

def steady_state_solution(x,y,t):
    output = 0
    if   x > 4./15: output = -Q3/2*x*x-f1*x-f2
    elif x < 2./15: output = -Q1/2*x*x-c1*x-c2
    else          : output = -Q2/2*x*x-d1*x-d2
    return output

h           = 0.01
x           = np.arange(0,0.8+h,h)
T_steady    = np.array([steady_state_solution(x_,0,0) for x_ in x])

T = np.loadtxt("post_enr_explicit_euler_h_0.010000.txt", 
               delimiter = ',',
               unpack    = True)

# scaled plots

#plt.figure("post_enr_explicit_euler")
#plt.title(r"After enrichment, $\overline{t}=1$")
#plt.axvline(x=2./15,linestyle = "dotted", color = "black")
#plt.axvline(x=4./15,linestyle = "dotted", color = "black")
#plt.plot(x, T[50,:], label = r"$\overline{y} = 0.5$")
#plt.plot(x, T[25,:], label = r"$\overline{y} \in \{ 0.25, 0.75\}$")
#plt.plot(x, T[10,:], label = r"$\overline{y} \in \{0.1, 0.9\}$")
#plt.text(0.2/15, 0.1, r"$\overline{Q}_1$")
#plt.text(3.0/15, 0.1, r"$\overline{Q}_2$")
#plt.text(10./15, 0.1, r"$\overline{Q}_3$")
#plt.plot(x, T_steady, label = "steady state (analytical)")
#plt.xlabel(r"$\overline{x}$")
#plt.ylabel(r"$\overline{T}$")
#plt.legend(loc="best")
#plt.show(block = False)

#X    = np.arange(0, 0.8+float(h), float(h))
#Y    = np.arange(0, 1  +float(h), float(h))

#plt.figure("post_enr_explicit_euler_color")
#plt.pcolor(X, Y, T, cmap='RdBu_r', vmin=np.min(T), vmax=np.max(T))
#plt.title(r"After enrichment, $\overline{t}=1$")
#plt.axvline(x=2./15,linestyle = "dotted", color = "black")
#plt.axvline(x=4./15,linestyle = "dotted", color = "black")
#plt.colorbar(label = r"$\overline{T}$")
#plt.xlabel(r"$\overline{x}$")
#plt.ylabel(r"$\overline{y}$")
#plt.show(block=False)

# unscaled plots

L  = 150  #km
T0 = 9000 # degrees celsius

plt.figure("post_enr_explicit_euler_unscaled")
plt.title(r"1 G yr after enrichment", size=14)
plt.axvline(x=2./15*L,linestyle = "dotted", color = "black")
plt.axvline(x=4./15*L,linestyle = "dotted", color = "black")
plt.plot(x*L, T[50,:]*T0, label = r"$y = 75$ km")
plt.plot(x*L, T[25,:]*T0, label = r"$y = 37.5$ km$, 112.5$ km")
plt.plot(x*L, T[10,:]*T0, label = r"$y = 15$ km$, 135$ km")
plt.text(0.2/15*L, 0.1*T0, r"$Q_1$")
plt.text(3.0/15*L, 0.1*T0, r"$Q_2$")
plt.text(10./15*L, 0.1*T0, r"$Q_3$")
plt.plot(x*L, T_steady*T0, label = "steady state (analytical)")
plt.xlabel(r"$x$ [km]", fontsize=14)
plt.ylabel(r"$T\,$ [$^\circ C$]", fontsize=14)
plt.legend(loc="best")
plt.show(block = False)

X    = np.arange(0, 0.8+float(h), float(h))
Y    = np.arange(0, 1  +float(h), float(h))

plt.figure("post_enr_explicit_euler_color_unscaled")
plt.pcolor(X*L, Y*L, T*T0, cmap='RdBu_r', vmin=np.min(T)*T0, vmax=np.max(T)*T0)
plt.title(r"1 G yr after enrichment", size=14)
plt.axvline(x=2./15*L,linestyle = "dotted", color = "black")
plt.axvline(x=4./15*L,linestyle = "dotted", color = "black")
plt.colorbar().set_label(label = r"$T\,$ [$^\circ C$]", size=14)
plt.xlabel(r"$x$ [km]", fontsize=14)
plt.ylabel(r"$y$ [km]", fontsize=14)
plt.show(block=False)

#=========================================================================
# Success message
#=========================================================================

print("Success!")
plt.show()


#back up code for 3D surface

#for h in ["0.100000", "0.010000"]:
#    for t in ["t1","t2"]:
#        filename = "two_dim_implicit_euler_h_"+h+"_"+t
#        u = np.loadtxt(filename+".txt",
#                       delimiter=',', 
#                       unpack=True)
#        m = np.size(u[:,0])
#        n = np.size(u[0,:])
        
#        u_analytical = np.zeros((m,n));
#        for i, x in enumerate(np.linspace(0,1,m)):
#            for j, y in enumerate(np.linspace(0,1,n)):
#                u_analytical[j][i] = analytical_2D(x,y,timelist[t])

#        fig = plt.figure(filename)
#        ax = fig.add_subplot(111, projection='3d')
#        ax.set_title(r"$t=$"+str(timelist[t])+r", $h=$"+h+r" $dt=$"+"0.01"+r"$h^2$"+"\n"
#                     + "rel err: " + str(relerror(u_analytical ,u))
#                     , fontsize=14)
#        X = np.arange(0, 1+float(h),float(h))
#        Y = np.arange(0, 1+float(h),float(h))
#        X, Y = np.meshgrid(X, Y)

#        # make 3D plot
#        surf = ax.plot_surface(X, Y, u, 
#                               cmap=cm.coolwarm,
#                               linewidth=0, antialiased=False)
#        ax.set_xlabel(r'$x$', fontsize=14)
#        ax.set_ylabel(r'$y$', fontsize=14)
#        ax.set_zlabel(r'$u$', fontsize=14)
#        ax.view_init(5, 244)






































#T_t1        = np.loadtxt("pre_enr_explicit_euler_dx_0.010000_t1.txt", 
#                         delimiter = ',', 
#                         unpack    = True)
#T_t2        = np.loadtxt("pre_enr_explicit_euler_dx_0.010000_t2.txt", 
#                         delimiter = ',', 
#                         unpack    = True)

#plt.figure("Steady state solution versus numerical simulation")
#plt.title("Steady state solution")
#plt.plot(x, T_steady, label = "steady state (analytical)")
#plt.plot(x, T_t1,     label = r"$\overline{t} = 0.05$")
#plt.plot(x, T_t2,     label = r"$\overline{t} = 1$")
#plt.xlabel(r"$\overline{x}$")
#plt.ylabel(r"$\overline{T}$")
#plt.legend(loc="best")
#plt.show(block=False)
#print("Relative error: ",  relerror(T_steady,T_t2))