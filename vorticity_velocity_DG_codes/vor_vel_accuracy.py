from dolfin import *


def compute(nx):
    
   
    # Creating mesh and defining function space
    mesh = UnitCircle(nx)
    W = FunctionSpace(mesh, "DG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    U = W * V
    Z = FunctionSpace(mesh,"CG", 1)
    R = FunctionSpace(mesh,"DG",0)
    
    u = TrialFunction(V)
    v = TestFunction(V)

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
    T = 3
    CFL = 3 
    h = 1.0/nx
    dt = h*CFL

	# Penalty parameter values.
    alpha = 1
    beta = 1
	
    print "Nx= %d, h = %.2f, dt = %.2f, alpha = %.2f, beta = %.2f" % (nx,h,dt,alpha,beta)

	# Functions
    w0 = Function(W)
    w1 = Function(W)
    wob = Function(W)
    u0 = Function(V)
    u1 = Function(V)

    

    set_log_level(40)

	# Determining vorticity and velocity at the first time step using 
	# backward Euler for vorticity.

	# First we compute initial velocity using the velocity equation

    inivel = inner(grad(u), grad(v))*dx -vorin*Dx(v[1],0)*dx + vorin*Dx(v[0],1)*dx+\
			 pow(alpha,2.0)*div(u)*div(v)*dx  +\
			 pow(beta,2.0)*(Dx(u[1],0)-Dx(u[0],1)-vorin)*(Dx(v[1],0)-Dx(v[0],1))*dx
	 
    a = lhs(inivel)
    L = rhs(inivel)

    solve(a == L, u0 , vbc)
    w0 = project(vorin, W)
    
    outfile = File("pat2Velocitypen.pvd")
    outfile2 = File("pat2Vorticitypen.pvd")
    outfile4 = File("patdiv.pvd")

	# Now solving vorticity w1 using backward Euler.

    w2 = TrialFunction(W)
    wt = TestFunction(W)
    f = Constant(0.0)

    n   = FacetNormal(mesh)
    un  = dot(u0,n)
    unp = 0.5*(un+abs(un))
    unm = 0.5*(un-abs(un))
    H   = unp('+')*w2('+') + unm('+')*w2('-')


    vorticity_form = (1/dt)*(w2-w0)*wt*dx - w2*inner(u0, grad(wt))*dx + H*jump(wt)*dS +\
		   unp*w2*wt*ds + unm*f*wt*ds
    a = lhs(vorticity_form)
    L = rhs(vorticity_form)

    solve(a==L, w1)

	# Solving for second step velocity

    velocity_form = inner(grad(u), grad(v))*dx -w1*Dx(v[1],0)*dx + w1*Dx(v[0],1)*dx+\
			 pow(alpha,2.0)*div(u)*div(v)*dx +\
			 pow(beta,2.0)*(Dx(u[1],0)-Dx(u[0],1)-w1)*(Dx(v[1],0)-Dx(v[0],1))*dx
    
    a = lhs(velocity_form)
    L = rhs(velocity_form)

    solve(a == L, u1 , vbc)

	#---------------------------------------------------------------------- 
	# Now calculating the vorticity and velocity for subsequent time steps.
    umod = 2*u1-u0
    un  = dot(umod,n)
    unp = 0.5*(un+abs(un))
    unm = 0.5*(un-abs(un))
    H   = unp('+')*w2('+') + unm('+')*w2('-')

    vorticty_form  = (3.0/(2*dt))*(w2-(4.0/3.0)*w1+(1.0/3.0)*w0)*wt*dx -\
			w2*inner(umod, grad(wt))*dx + H*jump(wt)*dS + unp*w2*wt*ds + unm*f*wt*ds 

    a = lhs(vorticity_form)
    L = rhs(vorticity_form)

    a1 = lhs(velocity_form)
    L1 = rhs(velocity_form)

    t = 2*dt

    while t <= T:

		 solve(a == L, wob)
		 
		 w0.assign(w1)
		 w1.assign(wob)
		 u0.assign(u1)
		 
		 solve(a1 == L1, u1 , vbc)
		 print "Time = ", t
		 
		
		 t += dt

	# Calculating error in velocity.
		 
    err_vel = errornorm(vel_ex, u1, 'L2')
    vel_ex=project(vel_ex,V)
    norm_exact_velocity = norm(vel_ex,'l2')
	
	# Calculating error in divergence of velocity.
    uob = Function(Z)
    uob = project(div(u1),Z)
    err_div = norm(uob,'l2')
	
	# Calculating error in vorticity
	
    err_vor = errornorm(vorin, w1, 'L2')
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
for nx in [ 64]:
    Nx.append(nx)
    E.append(compute(nx))  

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



    
    



  



