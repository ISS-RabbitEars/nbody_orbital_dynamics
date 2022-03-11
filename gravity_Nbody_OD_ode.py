import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from matplotlib import animation
import sys

def a_alpha(rt, t, p):
	rvto_list=rt
	m,G=p

	r=[]
	v=[]
	th=[]
	om=[]
	for i in range(m.size):
		r.append(rvto_list[4*i])
		v.append(rvto_list[4*i+1])
		th.append(rvto_list[4*i+2])
		om.append(rvto_list[4*i+3])
	a=[]
	for i in range(m.size):
		sumrdd=0
		sumthdd=0
		for j in range(m.size):
			if i!=j:
				cos=np.cos(th[i]-th[j])
				sin=np.sin(th[i]-th[j])
				rij=(r[i]**2+r[j]**2-2*r[i]*r[j]*cos)**(-3/2)
				fac=-2*G*m[j]*rij
				sumrdd+=fac*(r[i]-r[j]*cos)
				sumthdd+=fac*sin*(r[j]/r[i])
		a.append(v[i])
		a.append(r[i]*om[i]**2+sumrdd)
		a.append(om[i])
		a.append(-2*(v[i]/r[i])*om[i]+sumthdd)
	return a		
	
rng=np.random.default_rng(92314311)

N=10
mass_max=3
mass_min=1
size_factor=1/100
radius_max=5
radius_min=1
velocity_max=0
velocity_min=0
angular_velocity_max=1
angular_velocity_min=-1

m=(mass_max-mass_min)*np.random.rand(N)+mass_min
ro=(radius_max-radius_min)*np.random.rand(N)+radius_min
vo=(velocity_max-velocity_min)*np.random.rand(N)+velocity_min
thetao=2*np.pi*np.random.rand(N)
omegao=(angular_velocity_max-angular_velocity_min)*np.random.rand(N)+angular_velocity_min
G=1

p=[m,G]
rt=[]
for i in range(N):
	rt.append(ro[i])
	rt.append(vo[i])
	rt.append(thetao[i])
	rt.append(omegao[i])

tf = 120
nfps = 60
nframes = tf * nfps
t = np.linspace(0, tf, nframes)

rth = odeint(a_alpha, rt, t, args = (p,))

r=[]
th=[]
for i in range(0,4*N,4):
	r.append(rth[:,i])
	th.append(rth[:,i+2])

x=r*np.cos(th)
y=r*np.sin(th)

xamax=[]
xamin=[]
yamax=[]
yamin=[]
for i in range(N):
	xamax.append(max(x[i]))
	xamin.append(min(x[i]))
	yamax.append(max(y[i]))
	yamin.append(min(y[i]))

dx=max(xamax)-min(xamin)
dy=max(yamax)-min(yamin)
dr=np.sqrt(dx**2+dy**2)
maxmr=size_factor*dr
mmax=m.max()
mmax=1/mmax
f=maxmr*mmax
mr=[]
for i in range(N):
	mr.append(f*m[i])
shift=max(mr)

xmax=max(xamax)+2*shift
xmin=min(xamin)-2*shift
ymax=max(yamax)+2*shift
ymin=min(yamin)-2*shift

v=[]
w=[]
for i in range(1,4*N,4):
	v.append(rth[:,i])
	w.append(rth[:,i+2])

keTot=[]
peTot=[]
for i in range(nframes):
	keTot.append(0)
	peTot.append(0)

ke=[]
for i in range(N):
	ke.append(0.5*m[i]*((v[i]**2)+((r[i]*w[i])**2)))
	keTot+=ke[i]

for i in range(N-1):
	for j in range(i+1,N):
		rij=np.sqrt((r[i]**2)+(r[j]**2)-2*r[i]*r[j]*np.cos(th[i]-th[j]))
		for ii in range(nframes):
			if rij[ii]==0:
				sys.exit(1)
		peTot-=2*G*m[i]*m[j]/rij

ETot=keTot+peTot
#Emax=abs(max(ETot))
#keTot/=Emax
#peTot/=Emax
#ETot/=Emax
#Emax=max(ETot)
#keTot-=Emax
#peTot+=Emax
#ETot-=Emax

fig, a=plt.subplots()
fig.tight_layout()

def run(frame):
	plt.clf()
	plt.subplot(211)
	for i in range(N):
		circle=plt.Circle((x[i][frame],y[i][frame]),radius=mr[i],fc='r')
		plt.gca().add_patch(circle)
	plt.title(str(N)+" Body Orbital Dynamics")
	ax=plt.gca()
	ax.set_aspect(1)
	plt.xlim([xmin,xmax])
	plt.ylim([ymin,ymax])
	ax.xaxis.set_ticklabels([])
	ax.yaxis.set_ticklabels([])
	ax.xaxis.set_ticks_position('none')
	ax.yaxis.set_ticks_position('none')
	ax.set_facecolor('xkcd:black')
	plt.subplot(212)
	plt.plot(t[0:frame],keTot[0:frame],'r',lw=1)
	plt.plot(t[0:frame],peTot[0:frame],'b',lw=1)
	plt.plot(t[0:frame],ETot[0:frame],'g',lw=1)
	plt.xlim([0,tf])
	plt.title("Energy")
	ax=plt.gca()
	ax.legend(['T','V','E'],labelcolor='w',frameon=False)
	ax.set_facecolor('xkcd:black')


ani=animation.FuncAnimation(fig,run,frames=nframes)
writervideo = animation.FFMpegWriter(fps=nfps)
ani.save('gravity_Nbody_ode_wgphs.mp4', writer=writervideo)

#plt.show()


