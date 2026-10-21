"""
This file provides a parameterized implementation of the DT muon-system cluster
reconstruction efficiency for the CMS DT LLP-cluster search, as implemented in the
DtClusterEfficiency/DTClusterEff module of
cards/delphes_card_CMS_MDSCluster.tcl.

The efficiency maps are binned in (LLP EM_E, LLP Had_E, LLP decay-product multiplicity).
The efficiency maps depend also on the DT stations (MB2 vs MB34).
The LLP decay radius is used to select the DT stations of interest.
"""

import numpy as np

# multiplicity bin edges: [0,2), [2,7), [7,inf)
MULT_EDGES = [0.0, 2.0, 7.0]

MB2_MAPS = {
    0: {
        "Eem_edges": [0.0, 9.706285],
        "Ehad_edges": [0.0, 5.496247, 26.994982, 51.251625, 63.99268, 79.248497, 437.200378],
        "values": [
            [0.002, 0.073, 0.106, 0.129, 0.168, 0.166],
        ],
    },
    1: {
        "Eem_edges": [0.0, 5.674536, 18.46857, 36.249339, 52.044257, 64.906179, 81.518173, 797.095337],
        "Ehad_edges": [0.0, 9.723825, 25.879377, 485.057312],
        "values": [
            [0.02, 0.09, 0.23],
            [0.065, 0.13, 0.251],
            [0.125, 0.171, 0.255],
            [0.167, 0.179, 0.251],
            [0.2, 0.173, 0.231],
            [0.201, 0.223, 0.348],
            [0.183, 0.167, 0.281],
        ],
    },
    2: {
        "Eem_edges": [0.0, 5.449181, 9.906474, 14.366728, 19.392473, 25.706384, 36.042243, 604.979736],
        "Ehad_edges": [0.0, 26.697799, 36.106002, 43.363067, 50.42286, 58.515944, 72.289762, 746.49176],
        "values": [
            [0.132, 0.197, 0.256, 0.307, 0.327, 0.358, 0.334],
            [0.153, 0.254, 0.314, 0.3, 0.322, 0.327, 0.315],
            [0.174, 0.274, 0.289, 0.298, 0.341, 0.322, 0.327],
            [0.199, 0.271, 0.304, 0.32, 0.324, 0.341, 0.327],
            [0.219, 0.302, 0.313, 0.328, 0.325, 0.294, 0.341],
            [0.228, 0.288, 0.287, 0.326, 0.29, 0.319, 0.314],
            [0.23, 0.253, 0.274, 0.287, 0.249, 0.279, 0.346],
        ],
    },
}

MB34_MAPS = {
    0: {
        "Eem_edges": [0.0, 12.123251],
        "Ehad_edges": [0.0, 4.863147, 25.787655, 50.804605, 62.54412, 75.728984, 455.701416],
        "values": [
            [0.0, 0.03, 0.085, 0.111, 0.11, 0.144],
        ],
    },
    1: {
        "Eem_edges": [0.0, 5.456717, 17.911263, 36.554535, 51.7448, 63.398891, 77.948358, 719.496582],
        "Ehad_edges": [0.0, 9.284322, 25.028338, 426.601898],
        "values": [
            [0.008, 0.048, 0.171],
            [0.041, 0.086, 0.182],
            [0.094, 0.131, 0.196],
            [0.13, 0.143, 0.204],
            [0.15, 0.155, 0.218],
            [0.146, 0.142, 0.237],
            [0.149, 0.149, 0.235],
        ],
    },
    2: {
        "Eem_edges": [0.0, 5.335904, 9.624089, 13.935826, 18.765238, 24.840552, 34.457932, 379.359222],
        "Ehad_edges": [0.0, 25.828524, 34.966445, 41.976092, 48.682428, 55.938417, 67.775543, 820.406128],
        "values": [
            [0.064, 0.11, 0.17, 0.222, 0.268, 0.302, 0.337],
            [0.098, 0.179, 0.235, 0.279, 0.312, 0.31, 0.343],
            [0.109, 0.182, 0.238, 0.279, 0.337, 0.315, 0.36],
            [0.131, 0.212, 0.267, 0.291, 0.308, 0.298, 0.348],
            [0.16, 0.215, 0.276, 0.299, 0.296, 0.311, 0.357],
            [0.183, 0.243, 0.264, 0.27, 0.27, 0.291, 0.381],
            [0.19, 0.217, 0.251, 0.245, 0.253, 0.297, 0.39],
        ],
    },
}

EFF_MAPS = {"MB2": MB2_MAPS, "MB34": MB34_MAPS}

"""

The region function classifies the LLP decay radius (decayR, in mm, matching Delphes'
internal units and the tcl EfficiencyFormula) into the MB2 ([4500,5330) mm) or
MB34 ([5330,7380) mm) station region, or 'NA' if outside both (no DT cluster is possible there).

"""


def region(decayR):
    decayR = np.asarray(decayR, dtype=float)
    reg = np.full(decayR.shape, "NA", dtype=object)
    reg[(decayR >= 4500.0) & (decayR < 5330.0)] = "MB2"
    reg[(decayR >= 5330.0) & (decayR < 7380.0)] = "MB34"
    return reg


"""

The mult_bin function classifies the LLP decay-product multiplicity into bin index
0 ([0,2)), 1 ([2,7)), or 2 ([7,inf)), or -1 if negative (not physical).

"""


def mult_bin(multiplicity):
    multiplicity = np.asarray(multiplicity, dtype=float)
    mbin = np.full(multiplicity.shape, -1, dtype=int)
    mbin[(multiplicity >= MULT_EDGES[0]) & (multiplicity < MULT_EDGES[1])] = 0
    mbin[(multiplicity >= MULT_EDGES[1]) & (multiplicity < MULT_EDGES[2])] = 1
    mbin[multiplicity >= MULT_EDGES[2]] = 2
    return mbin


def _lookup(grid, eem, ehad):
    xedges = grid["Eem_edges"]
    yedges = grid["Ehad_edges"]
    xbin = min(max(np.searchsorted(xedges, eem, side="right") - 1, 0), len(xedges) - 2)
    ybin = min(max(np.searchsorted(yedges, ehad, side="right") - 1, 0), len(yedges) - 2)
    return grid["values"][xbin][ybin]


"""

input: list of LLP decay radius (decayR, mm), decay-product multiplicity, LLP EM energy (Eem,
GeV) and LLP hadronic energy (Ehad, GeV)
input type: decayR, multiplicity, Eem, Ehad are numpy arrays
output: array of DT cluster-reconstruction efficiency per event

"""


def dt_cluster_efficiency(decayR, multiplicity, Eem, Ehad):
    decayR = np.asarray(decayR, dtype=float)
    multiplicity = np.asarray(multiplicity, dtype=float)
    Eem = np.asarray(Eem, dtype=float)
    Ehad = np.asarray(Ehad, dtype=float)
    assert decayR.shape == multiplicity.shape == Eem.shape == Ehad.shape

    regions = region(decayR)
    mbins = mult_bin(multiplicity)

    eff = np.zeros(decayR.shape)
    for i in range(len(decayR)):
        reg = regions[i]
        mbin = mbins[i]
        if reg == "NA" or mbin < 0:
            eff[i] = 0.0
            continue
        eff[i] = _lookup(EFF_MAPS[reg][mbin], Eem[i], Ehad[i])
    return eff


if __name__ == "__main__":
    """
    dummy IO of input lists; 
    replace decayR, multiplicity, Eem, Ehad with real values from your ntuple.
    """

    decayR = np.array([4800.0, 5600.0, 6900.0, 3000.0])
    multiplicity = np.array([1, 3, 8, 5])
    Eem = np.array([12.4, 45.2, 6.1, 30.0])
    Ehad = np.array([15.0, 60.0, 40.0, 20.0])

    assert len(decayR) == len(multiplicity) == len(Eem) == len(Ehad)

    efficiencies = dt_cluster_efficiency(decayR, multiplicity, Eem, Ehad)
    for i in range(len(decayR)):
        print(
            f"decayR={decayR[i]:.1f} mm, mult={multiplicity[i]}, Eem={Eem[i]:.1f} GeV, "
            f"Ehad={Ehad[i]:.1f} GeV -> efficiency={efficiencies[i]:.3f}"
        )
