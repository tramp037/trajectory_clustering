"""
Unit and regression test for the trajclustering package.
"""

# Import package, test suite, and other packages as needed
import trajclustering
from trajclustering.clustering import GROMOS
import pytest
import sys

import MDAnalysis as mda
import MDAnalysisTests.datafiles as data
from MDAnalysis.analysis import rms
import numpy as np

@pytest.fixture
def u_gmx():
    return mda.Universe(data.TPR, data.XTC)

@pytest.fixture
def u_charmm():
    return mda.Universe(data.PSF, data.DCD)

def test_trajclustering_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "trajclustering" in sys.modules


def test_mdanalysis_logo_length(mdanalysis_logo_text):
    """Example test using a fixture defined in conftest.py"""
    logo_lines = mdanalysis_logo_text.split("\n")
    assert len(logo_lines) == 46, "Logo file does not have 46 lines!"

@pytest.mark.parametrize(
        "topology,trajectory",
        [
            (
                data.TPR,
                data.XTC
            ),
            (
                data.PSF,
                data.DCD
            )
        ]
)
def test_compute_rmsd(topology, trajectory):
    '''
    Test trajclustering.clustering.GROMOS._compute_rsmd()
    '''
    u = mda.Universe(topology, trajectory)
    ref = mda.Universe(topology, trajectory)
    u_atomgroup = u.select_atoms('backbone')
    ref_atomgroup = ref.select_atoms('backbone')
    cutoff = 2
    cluster = GROMOS(u_atomgroup, ref_atomgroup, cutoff)
    gromos_rmsd = cluster._compute_rmsd(u_atomgroup, ref_atomgroup)
    mda_rmsd = rms.rmsd(u_atomgroup.positions,
                        ref_atomgroup.positions,
                        center=True,
                        superposition=True)
    assert gromos_rmsd == pytest.approx(mda_rmsd)

@pytest.mark.parametrize(
        "topology,trajectory",
        [
            (
                data.TPR,
                data.XTC
            ),
            (
                data.PSF,
                data.DCD
            )
        ]
)
def test_gromos_run(topology, trajectory):
    '''
    Test trajclustering.clustering.GROMOS.run()
    '''
    u = mda.Universe(topology, trajectory)
    ref = mda.Universe(topology, trajectory)
    u_atomgroup = u.select_atoms('backbone')
    ref_atomgroup = ref.select_atoms('backbone')
    cutoff = 2
    cluster = GROMOS(u_atomgroup, ref_atomgroup, cutoff)
    cluster.run()
    assert cluster.matrix.shape == (u.trajectory.n_frames, u.trajectory.n_frames)
    # check that the diagonals are near zero
    assert cluster.matrix[0, 0] < 0.0001
    assert cluster.matrix[-1, -1] < 0.0001
    # check symmetry of the matrix
    assert np.allclose(cluster.matrix, cluster.matrix.T, rtol=0.001, atol=0.001)

from trajclustering import clustering

def test_structure_pool(u_gmx):
    u = u_gmx
    ref = u_gmx.copy()
    u_atomgroup = u.select_atoms('backbone')
    ref_atomgroup = ref.select_atoms('backbone')
    cluster = clustering.GROMOS(u.select_atoms("backbone"), ref.select_atoms("backbone"), cutoff=2)
    cluster = clustering.GROMOS(u_atomgroup, ref_atomgroup, cutoff=2)
    cluster.run()
    cluster._neighbor_count()
    assert len(cluster.cluster_groups.keys()) == 1

    cluster = GROMOS(u_atomgroup, ref_atomgroup, cutoff=1)
    cluster.run()
    assert len(cluster.cluster_groups.keys()) == 7
