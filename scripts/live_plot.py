"""Script for displaying and continuously updating a plot in another file.
Using this for experimenting with viasualization.
Requires a target file (specified as a command line argument when calling this file) to have a plot(ax) method.
As an example, run
python scripts/live_plot.py scripts/example_live_plot.py"""

import argparse
import importlib.util
import os
import sys
import time
import threading

import matplotlib.pyplot as plt
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler


def load_plot_function(module_path):
    """Dynamically load a plot(ax) function from the given Python file path."""
    module_name = os.path.splitext(os.path.basename(module_path))[0]

    if module_name in sys.modules:
        del sys.modules[module_name]

    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not hasattr(module, "plot"):
        raise AttributeError(f"{module_path} must define a `plot(ax)` function")

    return module


class PlotReloader(FileSystemEventHandler):
    def __init__(self, path, update_event):
        self.module_path = os.path.abspath(path)
        self.update_event = update_event

    def on_modified(self, event):
        if os.path.abspath(event.src_path) == self.module_path:
            self.update_event.set()


def main():
    parser = argparse.ArgumentParser(description="Live-reload matplotlib plot from a module.")
    parser.add_argument("module", help="Path to the Python file that defines a `plot(ax)` function.")
    args = parser.parse_args()

    module_path = os.path.abspath(args.module)
    update_event = threading.Event()

    plt.ion()
    fig, ax = plt.subplots()

    def update_plot():
        try:
            mod = load_plot_function(module_path)
            ax.clear()
            mod.plot(ax)
            fig.canvas.draw()
            fig.canvas.flush_events()
        except Exception as e:
            print(f"[Error updating plot]: {e}")

    # Do initial plot
    update_plot()

    # Start the watchdog observer in a separate thread
    observer = Observer()
    observer.schedule(PlotReloader(module_path, update_event), path=os.path.dirname(module_path), recursive=False)
    observer.start()

    print(f"Watching {module_path}. Press Ctrl+C to stop.")
    try:
        while True:
            plt.pause(0.1)  # Required to keep GUI responsive
            if update_event.is_set():
                update_event.clear()
                update_plot()
    except KeyboardInterrupt:
        observer.stop()
    observer.join()


if __name__ == "__main__":
    main()
