import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np

class GROMOS(AnalysisBase):
    def __init__(self, atomgroup, reference, **kwargs):
        # sets up the self.atomgroup and self.reference variables
        self.atomgroup = atomgroup
        self.reference = reference
        super().__init__(atomgroup.universe.trajectory,
                         **kwargs)

    def _prepare(self):
        # called before iteration on the trajectory has begun.
        # initialize results
        # =====================================================
        # here, probably make an array that assigns frames to a 
        # cluster 
        # OR probably 
        matrix = np.zeros((self.atomgroup.universe.trajectory.n_frames, 
                           self.atomgroup.universe.trajectory.n_frames))

    def _single_frame(self):
        # called after the trajectory is moved onto each new frame.
        pass
    
    def _conclude(self):
        # put everything together
        pass
    
    def _rmsd_to_other_frames(self):
        pass
    
    def _neighbor_count(self):
        pass
    
    def _compute_rmsd(self):
        # actually compute rmsd 
        pass
        