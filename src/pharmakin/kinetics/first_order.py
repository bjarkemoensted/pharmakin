import sympy

from pharmakin import parameters


@parameters.half_life.formula
def t_half(volume_of_distribution, clearance):
    res = sympy.log(2)*volume_of_distribution/clearance
    return res


if __name__ == '__main__':
    pass
