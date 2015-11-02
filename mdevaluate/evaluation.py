import matplotlib.pyplot as plot
import numpy as np
import argparse
from .simulation import GromacsSimulationResult
import os
import sys

def cli(evaluation):
    parser = argparse.ArgumentParser(description=evaluation.__doc__)
    parser.add_argument('--evaluate', default=True, type=int)
    parser.add_argument('--plot', default=True, type=bool)
    parser.add_argument('--mdrun-dir', default="mdrun_dir")
    parser.add_argument('--xtc-file')
    parser.add_argument('--gro-file')
    parser.add_argument('--index', action='append')
    parser.add_argument('--result-name')

    args = parser.parse_args()
    simulation = GromacsSimulationResult(args.mdrun_dir, xtc_file=args.xtc_file, 
        gro_file=args.gro_file,
        index_files=args.index)
    
    
    
    evaluation.results_dir = args.mdrun_dir
    evaluation.simulation = simulation
    if args.result_name:
        evaluation.base_filename = args.result_name

    print(args.evaluate)
    if args.evaluate:
        evaluation.evaluate()
        evaluation.save_results()
    if args.plot:
        evaluation.load_results()
        evaluation.plot()
        evaluation.save_plot()



class Evaluation:
    results = None
    simulation = None
    results_dir = None
    _base_filename = None

    @property
    def base_filename(self):
        if self._base_filename is not None:
            return self._base_filename
        return self.__class__.__name__

    @base_filename.setter
    def base_filename(self, b):
        self._base_filename = b

    @property
    def result_filename(self):
        return self.base_filename + ".npz"

    @property
    def plot_filename(self):
        return self.base_filename + ".png"


    def publish_results(self, **kwargs):
        if self.results is None:
            self.results = kwargs
        else:
            self.results.update(kwargs)

    def load_results(self):
        self.results = np.load(os.path.join(self.results_dir,self.result_filename))


    def save_results(self):
        np.savez(os.path.join(self.results_dir,self.result_filename), **self.results)

    def save_plot(self):
        plot.savefig(os.path.join(self.results_dir,self.plot_filename))

    def run(self):
        self.evaluate()
        self.save_results()
        self.plot()

    def evaluate(self):
        pass

    def plot(self):
        pass
