__author__ = 'mbartelm'
import os
from glob import glob

from pygmx.gromacs import atoms_from_grofile, load_indices
from pygmx.gromacs import XTCReader, TRRReader


class GromacsSimulationResult:
    __slots__ = ['simulation_dir', 'trajectory', 'full_trajectory', 'atoms', 'gro_file', 'xtc_file', 'trr_file',
                'index_files', 'indices']
    def __init__(self, sim_dir,
                 xtc_file=None,
                 gro_file=None,
                 trr_file=None,
                 index_files=None):

        self.simulation_dir = sim_dir

        self.index_files = glob(os.path.join(sim_dir, "*.ndx"))


        if index_files is not None:
            self.index_files += index_files
        self._load_index_files()

        if xtc_file is None:
            xtc_files = glob(os.path.join(sim_dir, "*.xtc"))
            if len(xtc_files) == 1:
                xtc_file = xtc_files[0]

        if gro_file is None:
            gro_files = glob(os.path.join(sim_dir, "*.gro"))
            if len(gro_files) == 1:
                gro_file = gro_files[0]

        assert gro_file, "No GRO File selected"
        assert xtc_file, "No XTC File selected"

        self.gro_file = gro_file
        self.xtc_file = xtc_file



        self.full_trajectory = None
        #if self.trr_file:
        #    self.full_trajectory = TRRReader(self.trr_file)

        if self.xtc_file:
            self.trajectory = XTCReader(self.xtc_file)
        self.atoms = atoms_from_grofile(self.gro_file)

    def _load_index_files(self):
        self.indices = {}
        for index_file in self.index_files:
            indices = load_indices(index_file)
            self.indices.update(indices)
