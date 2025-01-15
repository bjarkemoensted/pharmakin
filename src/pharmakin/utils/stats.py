from functools import partial
import numpy as np


def _compute_rayleigh_scale(mean, scale):
    """Computes the sigma (scale) parameter for a Rayleigh distribution, given the scale or mean."""
    
    if sum(arg is None for arg in (mean, scale)) != 1:
        raise ValueError
    
    if scale is None:
        scale = np.sqrt(2/np.pi)*mean
    
    return dict(scale=scale)


def _make_renamer(mean: str=None, scale: str=None):
    """Makes a function which takes mean and scale and returns a dictionary with the provided names for mean and scale.
    None can be used to exclude one of the parameters.
    For example, _make_renamer(mean="loc") returns a function which takes mean and scale as arguments,
    and returns a dictionary like {'loc': value_of_mean}."""

    mean_rename = mean
    scale_rename = scale
    def inner(mean, scale):
        res = dict()
        for new, old in zip((mean_rename, scale_rename), (mean, scale), strict=True):
            if new is not None:
                res[new] = old
            #
        return res
    return inner


class Simulator:
    """Callable for simulating values drawn from various distributions, using a specified mean and possibly scale,
    depending on the distribution."""

    distributions = dict(
        rayleigh=_compute_rayleigh_scale,
        normal=_make_renamer(mean="loc", scale="scale")
    )
    
    def __init__(self, distribution: str, mean: float=None, scale: float=None, seed: int=None):
        """distribution is a string matching a distribution from numpy.random.Generator.<distribution>.
        mean and scale are the desired mean and scale parameters for the distribution.
        seed is an optional interger seed for constructing the numpy PRNG."""
        
        self.generator = np.random.default_rng(seed=seed)
        self.distribution = distribution
        self._func = self._make_func(mean=mean, scale=scale)
    
    def _make_func(self, mean, scale):
        """Makes a function which can be called to simulate values from the specified distribution"""
        
        # Determine the correct way to parse the input mean and scale parameters
        try:
            parser = self.distributions[self.distribution]
        except KeyError:
            dists = ", ".join(self.distributions)
            raise ValueError(f"Invalid distribution: {self.distribution}. Supported: {dists}.")
    
        # Parse mean/scale into the required parameters for the distribution
        d = parser(mean=mean, scale=scale)
        # Create a callable which uses the PNRG to simulate values from the distribution
        f = partial(getattr(self.generator, self.distribution), **d)
        return f
    
    def __call__(self, size=None):
        """Draws a number/numbers from the distribution.
        size denotes the number of values to simulate (defaults to a single value)."""
        return self._func(size=size)


if __name__ == '__main__':
    simulate = Simulator("rayleigh", mean=2.0)
    print(np.mean(simulate()))
    