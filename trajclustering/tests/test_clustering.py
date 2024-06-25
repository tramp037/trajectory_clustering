from trajclustering import clustering
import MDAnalysis as mda
from MDAnalysisTests import datafiles as data

import pytest

@pytest.fixture
def u():
    return mda.Universe(data.GRO, data.XTC)

def test_gromos_init(u):
    u.select_atoms("backbone")
    ref = u.select_atoms("backbone")

    gromos = clustering.GROMOS(u.atoms, ref)
    assert gromos.atomgroup == u.atoms
    assert gromos.reference == ref

def test_gromos_prepare(u):
    u.select_atoms("backbone")
    ref = u.select_atoms("backbone")

    gromos = clustering.GROMOS(u.atoms, ref)
    gromos._prepare()
    assert gromos.matrix.shape == (u.trajectory.n_frames, u.trajectory.n_frames)
    assert gromos.matrix.sum() == 0