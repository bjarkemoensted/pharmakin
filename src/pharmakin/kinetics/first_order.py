import sympy

from pharmakin import parameters


@parameters.half_life.formula
def t_half(volume_of_distribution, clearance):
    res = sympy.log(2)*volume_of_distribution/clearance
    return res


def k_el(half_life):
    res = sympy.log(2)/half_life
    return res


def t_half_from_k(k: float):
    res = sympy.log(2)/k
    return res


if __name__ == '__main__':
    pass
