import argparse
import firedrake as fd
from firedrake import inner, dot, grad, div, dx, ds, pi
import firedrake_adjoint as fda
from pyadjoint.placeholder import Placeholder
from pyMMAopt import MMASolver, ReducedInequality
from filter import pde_filter
from penalization import ramp

DOMAIN = 5
INMOUTH = 6
OUTMOUTH = 7
INLET_WIDTH = 1.0
Y_INLET = 0.5


def navier_stokes_flow():

    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "--output-dir",
        dest="output_dir",
        type=str,
        action="store",
        default="./",
        help="Directory where to save the result",
    )

    opts = parser.parse_args()
    output_dir = opts.output_dir

    solver_parameters_direct = {
        "snes_atol": 1e-7,
        "mat_type": "aij",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": 200,
        "mat_mumps_icntl_24": 1,
    }

    boundary_conditions_forward = {
        "INLET": 1,
        "OUTLET": 2,
        "WALLS": 3,
        "ZERO_NORMAL": 4,
    }
    boundary_conditions_reverse = {
        "OUTLET": 1,
        "INLET": 2,
        "WALLS": 3,
        "ZERO_NORMAL": 4,
    }
    mesh = fd.Mesh("./diode.msh")

    RHO = fd.FunctionSpace(mesh, "DG", 0)
    V = fd.VectorFunctionSpace(mesh, "CG", 2)
    Q = fd.FunctionSpace(mesh, "CG", 1)
    W = V * Q
    ramp_p = fd.Constant(10.0)
    Re_val = 300.0
    Re = fd.Constant(Re_val)
    Da = fd.Constant(1e-4)
    Placeholder(Da)

    with fda.stop_annotating():
        rho = fd.interpolate(fd.Constant(0.5), RHO)

    rhof = pde_filter(
        RHO,
        rho,
        filter_radius=fd.Constant(1e-4),
        solver_parameters=solver_parameters_direct,
    )
    rhof_control = fda.Control(rhof)
    _, y = fd.SpatialCoordinate(W.ufl_domain())
    u_inflow = 1.0

    def inflow(u_inflow):
        return fd.as_vector(
            [
                -u_inflow * fd.sin(((y - (Y_INLET)) * pi) / INLET_WIDTH),
                0.0,
            ]
        )

    inflow_forward = inflow(u_inflow)
    up_forward = flow_problem(
        W,
        rhof,
        Re,
        Da,
        ramp_p,
        boundary_conditions_forward,
        inflow_forward,
        ramp,
        solver_parameters=solver_parameters_direct,
    )
    u_f, _ = fd.split(up_forward)

    inflow_reverse = inflow(-u_inflow)
    up_reverse = flow_problem(
        W,
        rhof,
        Re,
        Da,
        ramp_p,
        boundary_conditions_reverse,
        inflow_reverse,
        ramp,
        solver_parameters=solver_parameters_direct,
    )
    u_r, _ = fd.split(up_reverse)

    c = fda.Control(rho)

    plot_file = f"{output_dir}/design.pvd"
    controls_f = fd.File(plot_file)
    vel_pvd = fd.File("velocity.pvd")
    up_control_f = fda.Control(up_forward)
    up_control_r = fda.Control(up_reverse)
    rho_viz_f = fd.Function(RHO, name="rho")

    def deriv_cb(j, dj, rho):
        with fda.stop_annotating():
            rho_viz_f.assign(rhof_control.tape_value())
            u_plot, p_plot = up_control_f.tape_value().split()
            u_plot.rename("Velocity")
            p_plot.rename("Pressure")
            u_plot_r, p_plot_r = up_control_r.tape_value().split()
            u_plot_r.rename("Velocity reverse")
            p_plot_r.rename("Pressure reverse")
            vel_pvd.write(u_plot, p_plot, u_plot_r, p_plot_r)
            controls_f.write(rho_viz_f)

    G = power_dissipation(u_f, rhof, Re, Da, penalization=ramp, ramp_p=ramp_p)
    G_reverse = power_dissipation(u_r, rhof, Re, Da, penalization=ramp, ramp_p=ramp_p)
    Wf = 1.0
    Diodicity = fda.AdjFloat(1.0) / (G_reverse / G) + fd.assemble(
        Wf
        * alpha(rhof, Da, penalization=ramp, ramp_p=ramp_p)
        * inner(u_r, u_r)
        * dx(DOMAIN)
    )
    Ghat = fda.ReducedFunctional(Diodicity, c, derivative_cb_post=deriv_cb)

    Vol = fd.assemble((rhof) * dx(DOMAIN))
    Vhat = fda.ReducedFunctional(Vol, c)
    Vcontrol = fda.Control(Vol)
    with fda.stop_annotating():
        Vlimit = 0.7 * fd.assemble(fd.Constant(1.0) * dx(DOMAIN, domain=mesh))

    Da_arr = [1e-3, 5e-4, 2e-4]
    max_iter_arr = [50, 50, 1000]
    movlim_arr = [0.01, 0.01, 0.01]
    for Da_val, max_iter, movlim in zip(Da_arr, max_iter_arr, movlim_arr):
        Da.assign(Da_val)

        lb = 0.0
        ub = 1.0

        problem = fda.MinimizationProblem(
            Ghat,
            bounds=(lb, ub),
            constraints=[
                ReducedInequality(Vhat, Vlimit, Vcontrol),
            ],
        )

        parameters_mma = {
            "move": movlim,
            "maximum_iterations": max_iter,
            "m": 1,
            "norm": "L2",
            "gcmma": False,
        }
        solver = MMASolver(problem, parameters=parameters_mma)

        results = solver.solve()
        rho_opt = results["control"]
        with fda.stop_annotating():
            rho.assign(rho_opt)


def alpha(rho, Da, penalization=ramp, ramp_p=10.0):
    return (
        fd.Constant(1.0) / Da * penalization(rho, ramp_p=ramp_p, val_0=1.0, val_1=0.0)
    )


def GLS(u, v, p, q, rhof, Da, Re, penalization=ramp, ramp_p=10.0):
    """
    Outside of design domain DOMAIN:
    GLS = tau_gls * inner(R_U, theta_u) * dx(OUTSIDE)
    Inside of design domain
    GLS = tau_gls_alpha * inner(R_U + R_U_alpha, theta_U + theta_U_alpha) * dx(DOMAIN)
    """
    R_U = dot(u, grad(u)) - 1.0 / Re * div(grad(u)) + grad(p)
    R_U_alpha = alpha(rhof, Da, penalization=penalization, ramp_p=ramp_p) * u
    theta_U = dot(u, grad(v)) - 1.0 / Re * div(grad(v)) + grad(q)
    theta_U_alpha = alpha(rhof, Da, penalization=penalization, ramp_p=ramp_p) * v
    mesh = rhof.function_space().ufl_domain()
    h = fd.CellDiameter(mesh)

    beta_gls = 0.5
    # beta_gls = 30.0
    tau_gls = fd.Constant(beta_gls) * (
        (4.0 * dot(u, u) / h ** 2) + 9.0 * (4.0 / (Re * h ** 2)) ** 2
    ) ** (-0.5)
    tau_gls_alpha = fd.Constant(beta_gls) * (
        (4.0 * dot(u, u) / h ** 2)
        + 9.0 * (4.0 / (Re * h ** 2)) ** 2
        + (alpha(rhof, Da, penalization=penalization, ramp_p=ramp_p)) ** 2
    ) ** (-0.5)

    return tau_gls * inner(R_U, theta_U) * (
        dx(INMOUTH) + dx(OUTMOUTH)
    ) + tau_gls_alpha * inner(R_U + R_U_alpha, theta_U + theta_U_alpha) * dx(DOMAIN)


def flow_problem(
    W,
    rhof,
    Re,
    Da,
    ramp_p,
    boundary_conditions,
    inflow1,
    penalization,
    solver_parameters=None,
):

    v, q = fd.TestFunctions(W)

    up1 = fd.Function(W)
    u1, p1 = fd.split(up1)
    F1 = (
        1.0 / Re * inner(grad(u1), grad(v)) * dx
        + inner(dot(grad(u1), u1), v) * dx
        - p1 * div(v) * dx
        + div(u1) * q * dx
        + alpha(rhof, Da, penalization=penalization, ramp_p=ramp_p)
        * inner(u1, v)
        * dx(DOMAIN)
    )
    F1 = F1 + GLS(u1, v, p1, q, rhof, Da, Re, penalization=penalization, ramp_p=ramp_p)

    # Dirichelt boundary conditions
    noslip = fd.Constant((0.0, 0.0))
    bcs1_1 = fd.DirichletBC(W.sub(0), noslip, boundary_conditions["WALLS"])
    bcs1_2 = fd.DirichletBC(W.sub(0), inflow1, boundary_conditions["INLET"])
    bcs1_3 = fd.DirichletBC(
        W.sub(0).sub(1), fd.Constant(0.0), boundary_conditions["ZERO_NORMAL"]
    )
    bcs1 = [bcs1_1, bcs1_2, bcs1_3]

    problem = fd.NonlinearVariationalProblem(F1, up1, bcs=bcs1)
    solver = fd.NonlinearVariationalSolver(problem, solver_parameters=solver_parameters)
    solver.solve()

    return up1


def power_dissipation(u, rhof, Re, Da, penalization, ramp_p):
    return fd.assemble(
        1.0 / Re * inner(grad(u), grad(u)) * dx
        + alpha(rhof, Da, penalization=penalization, ramp_p=ramp_p)
        * inner(u, u)
        * dx(DOMAIN)
    )


if __name__ == "__main__":
    navier_stokes_flow()
