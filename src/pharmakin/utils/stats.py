from functools import partial
import numpy as np
from typing import Callable


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


# Map the names of some distributions (must match a numpy.Generator distribution name) to a method for
# converting mean + scale into the required kwargs for the distribution
_distribution_parsers = dict(
    rayleigh=_compute_rayleigh_scale,
    normal=_make_renamer(mean="loc", scale="scale")
)


def parse_distribution_mean_scale(distribution: str, mean: float=None, scale: float=None) -> dict:
    """Takes a distribution name, and values for mean/scale. Computes the parameters needed to make a distrution
    with the specified mean/scale.
    The resulting dict d can be passed to numpy.Generator.<distribution>(**d) to sample."""

    # Determine the correct way to parse the input mean and scale parameters
    try:
        parser = _distribution_parsers[distribution]
    except KeyError:
        dists = ", ".join(sorted(_distribution_parsers.keys()))
        raise ValueError(f"Invalid distribution: {distribution}. Supported: {dists}.")

    d = parser(mean=mean, scale=scale)
    return d


class Simulator:
    """Callable for simulating values drawn from various distributions, using a specified mean and possibly scale,
    depending on the distribution."""
    
    def __init__(
            self,
            distribution: str,
            lower_bound=None,
            upper_bound=None,
            seed: int=None,
            **kwargs):
        """distribution is a string matching a distribution from numpy.random.Generator.<distribution>.
        seed is an optional interger seed for constructing the numpy PRNG.
        lower_bound and upper_bound denote the lower and upper bounds on sampled values.
        Values drawn outside provided bounds are truncated to the provided limits.
        Note that this might affect the mean and scale of the distribution.
        kwargs are keyword arguments to be passed to numpy.Generator.<distribution>"""
        
        assert any(arg is None for arg in (lower_bound, upper_bound)) or upper_bound >= lower_bound
        self.generator = np.random.default_rng(seed=seed)
        self.distribution = distribution
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self._func = self._make_func(**kwargs)
    
    def _make_func(self, **kwargs):
        """Makes a function which can be called to simulate values from the specified distribution"""
        
        # Create a callable which uses the PNRG to simulate values from the distribution
        f = partial(getattr(self.generator, self.distribution), **kwargs)

        return f
    
    def __call__(self, size=None):
        """Draws a number/numbers from the distribution.
        size denotes the number of values to simulate (defaults to a single value)."""
        
        sample = self._func(size=size)
        if size is None:
            sample = [sample]
        
        # If bounds are provied, force values to be within bounds
        res = np.clip(a=sample, a_min=self.lower_bound, a_max=self.upper_bound)
        
        # If sampling a single value, convert back to native data types (for consistency with numpy.Generator.<dist>)
        if size is None:
            res = res[0].item()
        
        return res


if __name__ == '__main__':
    d = parse_distribution_mean_scale("rayleigh", mean=2.0)
    simulate = Simulator("rayleigh", **d)
    print(simulate())
    