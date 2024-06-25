import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.base import AnalysisBase

import numpy as np

class GROMOS(AnalysisBase):
    def __init__(self, atomgroup, **kwargs):
        self.atomgroup = atomgroup
        super().__init__(atomgroup.universe.trajectory,
                         **kwargs)
    def _prepare(self):
        # called before iteration on the trajectory has begun.
        # initialize results
        # =====================================================
        # here, probably make an array that assigns frames to a 
        # cluster 
        # OR probably 
        pass

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
    
    def _compute_rmsd(self, frame_atomgroup, ref_atomgroup):
        '''
        Helper function to compute RMSD between the current frame index and a reference.

        RMSD is calculated like this:

        Parameters
        ----------
        frame_atomgroup: AtomGroup
            Atom group for current frame
        ref_atomgroup: AtomGroup
            Atom group for reference frame

        Returns
        ----------
        rmsd : float

        '''
        # actually compute rmsd 
        return rms.rmsd(frame_atomgroup.positions,
                        ref_atomgroup.positions,
                        center=True,
                        superposition=True)
        