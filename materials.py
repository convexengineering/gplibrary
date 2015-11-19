"""Material library"""
from gpkit import Model, Variable


class Material(Model):
    # todo: what should be common?
    pass


class Aluminum6061(Material):
    def setup(self):
        rho = Variable(r"\rho", 2.70, "g/cm^3", "density")
        E = Variable("E", 69, "GPa", "Young's modulus")
        sigma_y = Variable(r"\sigma_{y}", 55, "MPa", "yield stress")
        # the objective and constraints are completely arbitrary.
        # use case is accessing Variables (and their associated values)
        return rho, [sigma_y <= E]


if __name__ == "__main__":
    pass
