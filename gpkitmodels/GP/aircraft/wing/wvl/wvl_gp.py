from numpy import dot, pi, linspace
from gpkit import Model, parse_variables, VectorVariable

class WVL(Model):
    """ Weissenger Vortex Lattice Model

    Variables
    ---------
    CL       1.1         [-]         lift coefficient
    CDi                  [-]         induced drag coefficient
    xCL                  [m^2]       segment placeholder lift coefficient
    S        8           [m^2]       surface area

    Variables of length N
    ---------------------
    G                    [m]         normalized vortex filament strength

    """

    def setup(self, N):
        exec parse_variables(WVL.__doc__)

        B = VectorVariable([N, N], "B", "-", "spanwise distribution matrix")
        dy = VectorVariable(N, "dy", linspace(0.5, 8, N), "m", "segment length")

        constraints = [CL <= 2*N*xCL/S,
                       xCL <= G*dy,
                       CDi >= 2*dot(G, dot(B, G))/S]

        for i in range(N):
            for j in range(N):
                constraints.append(B[i][j] >= dy[j]/pi/dy[i])

        return constraints

if __name__ == "__main__":
    m = WVL(2)
    m.cost = m["CDi"]
    sol = m.solve("mosek")
    print sol.table()
