from trajclustering import clustering
import MDAnalysis as mda
from MDAnalysisTests import datafiles as data

import pytest

@pytest.fixture
def u():
    return mda.Universe(data.GRO, data.XTC)

def test_gromos_init(u):
    atomgroup = u.select_atoms("backbone")
    ref = u.select_atoms("backbone")
    cutoff = 2

    gromos = clustering.GROMOS(atomgroup, ref, cutoff)
    assert gromos.atomgroup == atomgroup
    assert gromos.reference == ref
    assert gromos.cutoff == cutoff

def test_gromos_prepare(u):
    u.select_atoms("backbone")
    atomgroup = u.select_atoms("backbone")
    ref = u.select_atoms("backbone")
    cutoff = 2

    gromos = clustering.GROMOS(atomgroup, ref, cutoff)
    gromos._prepare()
    assert gromos.matrix.shape == (u.trajectory.n_frames, u.trajectory.n_frames)
    assert gromos.matrix.sum() == 0