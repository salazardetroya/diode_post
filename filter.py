import firedrake as fd
from firedrake import sqrt, jump, dS, dx


def pde_filter(RHO, rho, filter_radius=fd.Constant(0.2), solver_parameters=None):
    mesh = RHO.ufl_domain()
    x, y = fd.SpatialCoordinate(mesh)
    x_ = fd.interpolate(x, RHO)
    y_ = fd.interpolate(y, RHO)
    Delta_h = sqrt(jump(x_) ** 2 + jump(y_) ** 2)
    af, b = fd.TrialFunction(RHO), fd.TestFunction(RHO)
    aH = filter_radius * jump(af) / Delta_h * jump(b) * dS + af * b * dx
    LH = rho * b * dx

    rhof = fd.Function(RHO)
    fd.solve(aH == LH, rhof, solver_parameters=solver_parameters)

    return rhof
