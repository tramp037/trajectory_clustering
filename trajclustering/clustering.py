import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.base import AnalysisBase

import numpy as np

class GROMOS(AnalysisBase):
    """
    Class for GROMOS clustering.

    Parameters
    ----------
    atomgroup : AtomGroup
        First atom group
    reference : AtomGroup
        Second atom group
    cutoff : float
        Cutoff for neighbors

    Returns
    -------
    matrix : np.ndarray
        Matrix of frames and clusters
    """
    def __init__(self, atomgroup, reference, cutoff, **kwargs):
        # sets up the self.atomgroup and self.reference variables
        self.atomgroup = atomgroup
        self.reference = reference
        self.cutoff = cutoff
        super().__init__(atomgroup.universe.trajectory,
                         **kwargs)

    def _prepare(self):
        # called before iteration on the trajectory has begun.
        # initialize results
        # =====================================================
        # here, probably make an array that assigns frames to a 
        # cluster 
        # OR probably 
        self.matrix = np.zeros((self.atomgroup.universe.trajectory.n_frames, 
                           self.atomgroup.universe.trajectory.n_frames))

    def _single_frame(self):
        for ts in self.reference.universe.trajectory:
            self.matrix[self._frame_index, ts.frame] = self._compute_rmsd(self.atomgroup, self.reference)
    
    def _conclude(self):
        # put everything together
        pass
    
    def _rmsd_to_other_frames(self):
        pass
    
    def _neighbor_count(self):
        pass
    
    def _compute_rmsd(self, atomgroup, ref_atomgroup):
        '''
        Helper function to compute RMSD between the current frame index and a reference.

        RMSD is calculated like this:

        Parameters
        ----------
        atomgroup: AtomGroup
            Atom group for current frame
        ref_atomgroup: AtomGroup
            Atom group for reference frame

        Returns
        ----------
        rmsd : float

        '''
        return rms.rmsd(atomgroup.positions,
                        ref_atomgroup.positions,
                        center=True,
                        superposition=True)
        