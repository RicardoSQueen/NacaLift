from math import sin, cos, pi
import numpy as np
import matplotlib.pyplot as plt, mpld3
# import mpld3
from dolfin import *
from fenics import *
from mshr import *
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from matplotlib import ticker
# from plotter import *
from naca import return_naca_msh

class BladeMesh():
    def __init__(self, profNaca, n_of_points_naca, length, diameter,flow_speed, alpha = 5, delta=1, rho = 1.292, **kwargs):
        super().__init__(**kwargs)
        self.delta = delta # spacing between points
        self.length = length
        self.diameter = diameter
        self.profNaca = profNaca
        self.n_of_points = n_of_points_naca
        self.nx = int(length/delta)
        self.ny = int(diameter/delta)
        self.alpha = alpha*np.pi/180
        self.flow_speed = flow_speed
        self.rho = rho
        return None

    def create_mesh(self, mesh_resolution = 50):
        p0 = Point(np.array([0.0, 0.0]))
        p1 = Point(np.array([self.length, self.diameter]))
        X, Y = return_naca_msh(self.profNaca, self.n_of_points)
        aero_points = np.vstack((X, Y)).T
        baricenter = np.mean(aero_points, axis=0)
        rmatrix = np.matrix([[np.cos(self.alpha), np.sin(self.alpha)], [-np.sin(self.alpha), np.cos(self.alpha)]])
        # baricenter_rotated = np.dot(rmatrix, baricenter.T).T
        aero_points_rotated = [np.dot(rmatrix, (aeropoint.T)) for aeropoint in aero_points]
        aero_points_rotated = np.squeeze(np.asarray(aero_points_rotated))
        baricenter_rotated = np.mean(aero_points_rotated, axis=0)
        disloc_x = self.length/2 - baricenter_rotated[0]
        disloc_y = self.diameter/2 - baricenter_rotated[1]
        aerofoil = [Point(i[0] + disloc_x, i[1] + disloc_y) for i in aero_points_rotated]
        self.original_baricenter = [baricenter[0], baricenter[1]]
        self.baricenter = [disloc_x + baricenter_rotated[0], disloc_y + baricenter_rotated[1]]
        aerofoil = Polygon(aerofoil)
        domain = Rectangle(p0, p1) - aerofoil
        self.mesh = generate_mesh(domain, mesh_resolution)
        return (200)

    def create_formulation(self, circ_cc=None):
        #Creating Mesh according to geometry
        self.V = FunctionSpace(self.mesh, 'P', 1)
        # Define boundary condition expression
        expr_str = str(self.flow_speed)+'*x[1]'
        # print(expr_str)
        self.u_D = Expression(expr_str, degree=1)
        self.leading_edge_radius = 1.1019*(float(self.profNaca[-2:])/100)
        if circ_cc is None:
            self.circ_cc = 4*np.pi*self.flow_speed*self.leading_edge_radius*np.sin(self.alpha + abs(np.arcsin(self.original_baricenter[1]/self.leading_edge_radius)))*self.diameter
            self.circ_cc = abs((self.circ_cc/2*np.pi)*np.log(self.leading_edge_radius))/4
            print(self.circ_cc)
        else:
            self.circ_cc = circ_cc
        # Making a mark dictionary for bcs
        mark = {"generic": 0,"lower_wall": 1,"upper_wall": 2,"left": 3,"right": 4, "airfoil": 5 }
        self.subdomains = MeshFunction("size_t", self.mesh, 1)
        self.subdomains.set_all(mark["generic"])
        length = self.length
        diameter = self.diameter
        class Left(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[0], 0)
        class Right(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[0], length)
        class UpperWall(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[1], diameter)
        class LowerWall(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[1], 0)
        class Airfoil(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and not (near(x[1], 0) or near(x[1], diameter) or near(x[0], length) or near(x[0], 0))
        # Marking the self.subdomains
        left = Left()
        left.mark(self.subdomains, mark["left"])
        right = Right()
        right.mark(self.subdomains, mark["right"])
        upper_wall = UpperWall()
        upper_wall.mark(self.subdomains, mark["upper_wall"])
        lower_wall = LowerWall()
        lower_wall.mark(self.subdomains, mark["lower_wall"])
        airfoil = Airfoil()
        airfoil.mark(self.subdomains, mark["airfoil"])
        bc_lower_wall = DirichletBC(self.V, 0, self.subdomains, mark["lower_wall"])
        bc_upper_wall = DirichletBC(self.V, self.u_D((0,diameter)), self.subdomains, mark["upper_wall"])
        bc_left = DirichletBC(self.V, self.u_D, self.subdomains, mark["left"])
        bc_right = DirichletBC(self.V, self.u_D, self.subdomains, mark["right"])
        bc_airfoil = DirichletBC(self.V, self.circ_cc, self.subdomains, mark["airfoil"])
        self.bcs = [bc_lower_wall, bc_upper_wall, bc_left, bc_right, bc_airfoil]
        # Define variational problem
        self.u_trial = TrialFunction(self.V)
        self.v = TestFunction(self.V)
        self.f = Constant(0)
        self.a = dot(grad(self.u_trial), grad(self.v))*dx
        self.L = self.f*self.v*dx
        return (200)
    def compute_solution(self):
        # Compute solution
        self.u = Function(self.V)
        solve(self.a == self.L, self.u, self.bcs)
        self.W = VectorFunctionSpace(self.mesh, 'P', 1)
        self.vx = project(self.u.dx(1), self.V)
        self.vy = project(-self.u.dx(0), self.V)
        self.vvec = project(as_vector([self.vx, self.vy]),self.W)
        return (200)
    def calculate_lift(self):
        ds = Measure("ds", domain=self.mesh, subdomain_data=self.subdomains)
        GammaP = ds(5)
        n = FacetNormal(self.mesh)
        tangent = as_vector([n[1], -n[0]])
        self.circ = -assemble((dot(self.vvec, tangent))*GammaP)
        l1 = self.rho*self.vy*dot(self.vvec,n)*GammaP
        lift=assemble(l1)
        lift = self.rho*self.flow_speed*self.circ
        return self.circ, lift, l1
    def visualize_solutions(self, plot_type):
        if not self.u:
            return('Please compute_solution() first!')
        else:
            aspect = 20
            pad_fraction = 0.5
            plt.figure(figsize=(12,12))
            ax = plt.gca()
            plt.title('Função Corrente \n $\\psi = {}$'.format(self.u_D((0,self.diameter))))
            plt.ylabel('y \n $\\psi = {}y$'.format(self.u_D((0,self.diameter))))
            plt.xlabel('$\\psi = 0$ \n x')
            im = plot(self.u, cmap='jet')
            divider = make_axes_locatable(ax)
            width = axes_size.AxesY(ax, aspect=1./aspect)
            pad = axes_size.Fraction(pad_fraction, width)
            cax = divider.append_axes("right", size=width, pad=pad)
            plt.colorbar(im,cax=cax)
            fig1 = plt.gcf()
            #fig2
            plt.figure(figsize=(12,12))
            plt.title('Velocidade vetorial ($\\vec{v}$) $[m/s]$ \n')
            plt.xlabel('x $[m]$')
            plt.ylabel('y $[m]$')
            plt.xlim(self.baricenter[0] - 0.8, self.baricenter[0] + 0.8)
            plt.ylim(self.baricenter[1] - 0.8, self.baricenter[1] + 0.8)
            ax = plt.gca()
            im = plot(self.vvec)
            divider = make_axes_locatable(ax)
            width = axes_size.AxesY(ax, aspect=1./aspect)
            pad = axes_size.Fraction(pad_fraction, width)
            cax = divider.append_axes("right", size=width, pad=pad)
            cb = plt.colorbar(im, cax=cax)
            tick_locator = ticker.MaxNLocator(nbins=7)
            cb.locator = tick_locator
            cb.update_ticks()
            # plt.savefig('./Resultados/correnteaerovx.png', dpi=500, bbox_inches='tight')
            fig2 = plt.gcf()
            plt.figure(figsize=(12,12))
            plt.title('Velocidade em x ($v_x$) \n')
            plt.xlabel('x $[m]$')
            plt.ylabel('y $[m]$')
            ax = plt.gca()
            im = plot(self.u.dx(1))
            divider = make_axes_locatable(ax)
            width = axes_size.AxesY(ax, aspect=1./aspect)
            pad = axes_size.Fraction(pad_fraction, width)
            cax = divider.append_axes("right", size=width, pad=pad)
            plt.colorbar(im, cax=cax)
            fig3 = plt.gcf()
            
            plt.figure(figsize=(12,12))
            plt.title('Velocidade em y ($v_y$) \n')
            plt.xlabel('x $[m]$')
            plt.ylabel('y $[m]$')
            ax = plt.gca()
            im = plot(-self.u.dx(0))
            divider = make_axes_locatable(ax)
            width = axes_size.AxesY(ax, aspect=1./aspect)
            pad = axes_size.Fraction(pad_fraction, width)
            cax = divider.append_axes("right", size=width, pad=pad)
            plt.colorbar(im, cax=cax)
            fig4 = plt.gcf()
            
            return fig1, fig2, fig3, fig4



# naca2412 = BladeMesh('2412', 240, 0.1, 4, 1,100)
# naca2412.create_mesh(mesh_resolution=50)
# naca2412.create_formulation()
# naca2412.compute_solution()
# naca2412.visualize_solutions(plot_type='matplotlib')

# naca01420 = BladeMesh('01420', 240, 0.1, 4, 1,100)
# naca01420.create_mesh(mesh_resolution=50)
# naca01420.create_formulation()
# naca01420.compute_solution()
# naca01420.visualize_solutions()
