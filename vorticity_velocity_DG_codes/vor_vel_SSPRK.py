"""

Solving the vorticity velocity equation for Euler using explicit SSPRK2/SSPRK3
time stepping schemes and DG spacial discretization.

Authors: Praveen C., Souvik Roy.

"""


from dolfin import *

def wform(wtemp, wt, u0, n):
	un  = dot(u0,n)
	unp = 0.5*(un+abs(un))
	unm = 0.5*(un-abs(un))
	H   = unp('+')*wtemp('+') + unm('+')*wtemp('-')
	f = Constant(0.0)
	F1  = -wtemp*inner(u0, grad(wt))*dx + H*jump(wt)*dS + unp*wtemp*wt*ds + unm*f*wt*ds
	return F1 

def compute(nx,deg):
	if deg == 3:
		ark = [0.0, 3.0/4.0, 1.0/3.0]
		brk = [1.0, 1.0/4.0, 2.0/3.0]
	else:
		ark = [0.0, 0.5]
		brk = [1.0, 0.5]
	

	# Creating mesh and defining function space
	mesh = UnitCircle(nx)
	W = FunctionSpace(mesh, "DG", deg-1)
	V = VectorFunctionSpace(mesh, "CG", deg-1)
	U = W * V
	Z = FunctionSpace(mesh,"CG", 1)
	R = FunctionSpace(mesh,"DG",0)

	# Functions
	u = TrialFunction(V)
	v = TestFunction(V)
	w0 = Function(W)
	w1 = Function(W)
	wtemp = Function(W)
	wt = TestFunction(W)
	wob = Function(W)
	u0 = Function(V)
	u1 = Function(V)
	rhs_vec = Function(W)

	# Perelman Vortex: initial vorticity.
	class vortex(Expression):
		def eval(self, value, x):
			r = sqrt(pow((x[0]),2)+ pow((x[1]),2))
			if(r <= 0.5):
				value[0] = pow((1-r*r/0.25),7)
			else:
				value[0] = 0

	vorin = vortex()

	# Perelman vortex: steady velocity field.

	class velocity_exact(Expression):
		def eval(self, value, x):
			r = sqrt(pow((x[0]),2)+ pow((x[1]),2))
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

	vel_ex = velocity_exact()

	# Boundary condition.
	def Boundary(x, on_boundary):
		return on_boundary

	class gbound(Expression):
		 def eval(self, value, x):
				 r = sqrt(pow((x[0]),2)+ pow((x[1]),2))
				 value[0] = -1/(16*r*r/0.25)*x[1]
				 value[1] = 1/(16*r*r/0.25)*x[0]
		 def value_shape(self):
				 return (2,)

	g = gbound()
					  

	vbc = DirichletBC(V, g, Boundary)


	# Parameter values 
	T = 1
	CFL = 1.0/deg
	h = 1.0/nx
	dt = h*CFL

	# Penalty parameter values.
	alpha = 1
	beta = 1

	print "Nx= %d, h = %.2f, dt = %.2f, alpha = %.2f, beta = %.2f" % (nx,h,dt,alpha,beta)

	set_log_level(40)

	w0 = project(vorin, W)

	velocity_form = inner(grad(u), grad(v))*dx -w0*Dx(v[1],0)*dx + w0*Dx(v[0],1)*dx+\
	pow(alpha,2.0)*div(u)*div(v)*dx + pow(beta,2.0)*(Dx(u[1],0)-Dx(u[0],1)-w0)*(Dx(v[1],0)-Dx(v[0],1))*dx

	a1 = lhs(velocity_form)
	L1 = rhs(velocity_form)

	solve(a1==L1, u0,vbc)


	#---------------------------------------------------------------------- 
	# Now calculating the vorticity and velocity for subsequent time steps.
	wtemp.assign(w0)
	n = FacetNormal(mesh)


	w2 = TrialFunction(W)

	# rhs of the problem
	vorticity_form_rhs = -wform(wtemp,wt,u0,n)
	# Computing mass matrix
	m = w2*wt*dx 

	A = assemble(m)
	solver = LUSolver(A)
	solver.parameters['reuse_factorization']=True

	velocity_form = inner(grad(u), grad(v))*dx -wtemp*Dx(v[1],0)*dx + wtemp*Dx(v[0],1)*dx+\
	pow(alpha,2.0)*div(u)*div(v)*dx + pow(beta,2.0)*(Dx(u[1],0)-Dx(u[0],1)-wtemp)*(Dx(v[1],0)-Dx(v[0],1))*dx

	a1 = lhs(velocity_form)
	L1 = rhs(velocity_form)

	A1 = assemble(a1)
	vbc.apply(A1)
	solver_vel = LUSolver(A1)
	solver_vel.parameters['reuse_factorization']=True


	t = 0

	while t <= T:
		for i in range(deg):
			b = assemble(vorticity_form_rhs)
			solver.solve(rhs_vec.vector(),b)
			wtemp.vector()[:] = ark[i]*w0.vector()+ brk[i]*(wtemp.vector()+dt*rhs_vec.vector())
			
			b1 = assemble(L1)
			vbc.apply(b1)
			solver_vel.solve(u0.vector(),b1)
		w0.assign(wtemp)
		print "Time = ", t
		t += dt

	# Calculating error in velocity.
		 
	err_vel = errornorm(vel_ex, u0, 'L2')
	vel_ex=project(vel_ex,V)
	norm_exact_velocity = norm(vel_ex,'l2')

	# Calculating error in divergence of velocity.
	uob = Function(Z)
	uob = project(div(u0),Z)
	err_div = norm(uob,'l2')

	# Calculating error in vorticity

	err_vor = errornorm(vorin, w0, 'L2')
	vorin = project(vorin,W)
	norm_exact_vorticity = norm(vorin,'l2')

	# Collect error measures in a dictionary with self-explanatory keys
	errors = {	'Relative L2 error in velocity': err_vel/norm_exact_velocity,
				'Divergence error': err_div,
				'Relative L2 error in vorticity':err_vor/norm_exact_vorticity}

	return errors


# Implementing the numerical scheme for various mesh divisions.

Nx = [] # Mesh divisions
E = []  # Errors
for nx in [ 16, 32, 64]:
	Nx.append(nx)
	E.append(compute(nx,3))  

	# Computing convergence rates
	from math import log as ln  
	error_types = E[0].keys()
for error_type in sorted(error_types):
	Ei = E[0][error_type]
	print '\n', error_type
	print '\nNx \tRelative error \t\t Rate of convergence\n'
	print '%d \t %8.2E \t\t - '% (Nx[0], Ei)
	for i in range(1, len(E)):
		Ei   = E[i][error_type] 
		Eim1 = E[i-1][error_type]
		r = ln(Eim1/Ei)/ln(Nx[i]/Nx[i-1])	
		print '%d \t %8.2E \t\t %.2f' % (Nx[i], Ei, r)



    
    



  



