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
    def __init__(self, profNaca, n_of_points_naca, length, diameter,bc_upperwall_value, alpha = 5, delta=1, **kwargs):
        super().__init__(**kwargs)
        self.delta = delta # spacing between points
        self.length = length
        self.diameter = diameter
        self.profNaca = profNaca
        self.n_of_points = n_of_points_naca
        self.nx = int(length/delta)
        self.ny = int(diameter/delta)
        self.alpha = alpha*np.pi/180
        self.bc_upperwall_value = bc_upperwall_value
        return None

    def create_mesh(self, mesh_resolution = 50):
        parameters['allow_extrapolation'] = True
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
        self.baricenter = [disloc_x + baricenter_rotated[0], disloc_y + baricenter_rotated[1]]
        aerofoil = Polygon(aerofoil)
        domain = Rectangle(p0, p1) - aerofoil
        self.mesh = generate_mesh(domain, mesh_resolution)
        return (200)

    def create_formulation(self):
        #Creating Mesh according to geometry
        self.V = FunctionSpace(self.mesh, 'P', 1)
        # Define boundary condition expression
        expr_str = str(self.bc_upperwall_value)+'*x[1]'
        # print(expr_str)
        self.u_D = Expression(expr_str, degree=1)
        # Making a mark dictionary for bcs
        mark = {"generic": 0,"lower_wall": 1,"upper_wall": 2,"left": 3,"right": 4, "airfoil": 5 }
        subdomains = MeshFunction("size_t", self.mesh, 1)
        subdomains.set_all(mark["generic"])
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
        # Marking the subdomains
        left = Left()
        left.mark(subdomains, mark["left"])
        right = Right()
        right.mark(subdomains, mark["right"])
        upper_wall = UpperWall()
        upper_wall.mark(subdomains, mark["upper_wall"])
        lower_wall = LowerWall()
        lower_wall.mark(subdomains, mark["lower_wall"])
        airfoil = Airfoil()
        airfoil.mark(subdomains, mark["airfoil"])
        bc_lower_wall = DirichletBC(self.V, 0, subdomains, mark["lower_wall"])
        bc_upper_wall = DirichletBC(self.V, self.u_D((0,diameter)), subdomains, mark["upper_wall"])
        bc_left = DirichletBC(self.V, self.u_D, subdomains, mark["left"])
        bc_right = DirichletBC(self.V, self.u_D, subdomains, mark["right"])
        bc_airfoil = DirichletBC(self.V, self.u_D((self.baricenter[0],self.baricenter[1])), subdomains, mark["airfoil"])
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
# def calculate_lift():

    def visualize_solutions(self, plot_type):
        if not self.u:
            return('Please compute_solution() first!')
        else:
            aspect = 20
            pad_fraction = 0.5
            rho = 1.292
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
            if plot_type == 'matplotlib':
                plt.show()
            elif plot_type == 'web':
                fig1 = plt.gcf()
            #fig2
            plt.figure(figsize=(12,12))
            plt.title('Velocidade vetorial ($\\vec{v}$) $[m/s]$ \n')
            plt.xlabel('x $[m]$')
            plt.ylabel('y $[m]$')
            plt.xlim(self.baricenter[0] - 0.6, self.baricenter[0] + 0.6)
            plt.ylim(self.baricenter[1] - 0.6, self.baricenter[1] + 0.6)
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
            if plot_type == 'matplotlib':
                plt.show()
            elif plot_type == 'web':
                fig2 = plt.gcf()
                return fig1, fig2



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
