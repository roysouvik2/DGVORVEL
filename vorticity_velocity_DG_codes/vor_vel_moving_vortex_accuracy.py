from dolfin import *


nx = 64

# Creating mesh and defining function space
mesh = UnitCircle(nx)

W = FunctionSpace(mesh, "DG", 1)
V = VectorFunctionSpace(mesh, "CG", 1)
U = W * V
Z = FunctionSpace(mesh,"CG", 1)
R = FunctionSpace(mesh,"DG",0)

u = TrialFunction(V)
v = TestFunction(V)

up = Function(V)
vp = Function(W)

# Parameter values 
T = 1
CFL = 3 
h = 1.0/nx
dt = h*CFL
t = dt

outfile = File("Exact_velocity.pvd")
outfile2 = File("Exact_vorticity.pvd")

while t <= T:
	a = 0.5
	b = 0.0

	

	# Perelman Vortex: initial vorticity.
	class vortex(Expression):
		def __init__(self, t):
			self.t = t
		def eval(self, value, x):
			r = sqrt(pow((x[0]-a*self.t),2)+ pow((x[1]-b*self.t),2))
			if(r <= 0.5):
				value[0] = pow((1-r*r/0.25),7)
			else:
				value[0] = 0

	vorin = vortex(t)

	# Perelman vortex: steady velocity field.

	class velocity_exact(Expression):
		def __init__(self, t):
			self.t = t
		def eval(self, value, x):
			r = sqrt(pow((x[0]-a*self.t),2)+ pow((x[1]-b*self.t),2))
			if(r == 0):
				value[0]=0
				value[1]=0
			else: 
				if(r <= 0.5):
				   value[0] = -1/(16*r*r/0.25)*(1-pow((1-r*r/0.25),8))*x[1]
				   value[1] = 1/(16*r*r/0.25)*(1-pow((1-r*r/0.25),8))*x[0]
				else:
				   value[0] = -1/(16*r*r/0.25)*x[1]
				   value[1] = 1/(16*r*r/0.25)*x[0]
		def value_shape(self):
			   return (2,)

	vel_ex = velocity_exact(t)
	vel_ex = interpolate(vel_ex,V)
	up.assign(vel_ex)
	#plot(up)
	
	vor_ex = interpolate(vorin,W)
	vp.assign(vor_ex)
	#plot(vp)
	
	outfile << vel_ex
	outfile2 << vor_ex
	
	print "Time = ", t
	t += dt











