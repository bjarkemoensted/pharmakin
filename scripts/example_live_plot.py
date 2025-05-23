import numpy as np


def plot(ax):
    x = np.linspace(0, 12, 1000)
    y = np.sin(x)
    ax.plot(x, y)
    ax.set_title(f"Dynamic Reloaded Plott!")
