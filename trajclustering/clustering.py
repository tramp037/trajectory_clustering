import MDAnalysis as mda
from MDAnalysis.analysis import rms
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
        self.matrix = np.zeros((self.atomgroup.universe.trajectory.n_frames, 
                           self.atomgroup.universe.trajectory.n_frames))
        self.cluster_groups = {
            0:np.array([])
        } # not sure if this will be a dict or not but just for now 

    def _single_frame(self):
        pass
    
    def _conclude(self):
        # put everything together
        pass
    
    def _rmsd_to_other_frames(self):
        pass
    
    
    def _neighbor_count(self, structure_pool=None):
        '''
        Helper function that assigns frames to clusters. 

        Function finds all neighbors for each frame within a certain distance cutoff. The largest group is assigned to
        GROMOS.cluster_groups (dict) and removed from the pool of available structures. This function is ran recursively 
        until no more clusters are found.

        Parameters 
        -----------
        structure_pool : optional, None, np.ndarray 
            Indices of frames that are in the available structure pool
        
        Returns
        -----------
        structure_pool : np.ndarray
        '''
        if structure_pool is None:
            structure_pool = np.arange(0, self.atomgroup.n_atoms, 1) # create pool of initial frame indices
        max_neighbors = 0
        neighbors = np.array([])
        for i, row in enumerate(self.matrix):
            if i not in structure_pool:
                continue
            row_ = np.array([row[i] for i in range(len(row)) if i in structure_pool])
            n_neighbors = len(row_[row_ < self.cutoff])
            if n_neighbors > max_neighbors:
                max_neighbors = n_neighbors
                neighbors = np.where(row_[(row_ < self.cutoff)])
        cluster_members = neighbors[:]
        cluster_index = min(list(self.cluster_groups.keys())) + 1
        self.cluster_groups[cluster_index] = cluster_members
        return np.delete(structure_pool, cluster_members)
    
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
        