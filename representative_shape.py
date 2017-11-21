'''
Heuristic for extracting a representative shape from a cloud of 2d points
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.signal import savgol_filter

def representative_shape(pts):

    # first we select points that are near the edge and distributed radially
    com = np.mean(pts,axis=0)
    pts_dist = pts - com
    angles = [x+360 if x < 0 else x for x in np.degrees(np.arctan2(pts_dist[:,0], pts_dist[:,1]))]
    distance = np.linalg.norm(pts_dist, axis=1)
    small_bound = np.percentile(distance,0)
    upper_bound = np.percentile(distance,100)

    delta_angle = 5
    angle_range = np.arange(0,360,delta_angle)
    select_pts = [None for x in angle_range]

    for index, angle in enumerate(angle_range):
        for pt_index, pt_angle in enumerate(angles):
            if angle < pt_angle < angle+5:
                if select_pts[index] is None:
                    if small_bound < distance[pt_index] < upper_bound:
                        select_pts[index] = pt_index
                elif upper_bound > distance[pt_index] > distance[select_pts[index]]:
                        select_pts[index] = pt_index

    shape = list(filter(lambda x: x is not None, select_pts))
    shape.append(shape[0])
    (x,y)=(pts[shape,0],pts[shape,1])

    # apply Savitsky-Golay filter to smooth points
    fx = savgol_filter(x,5,1)
    fy = savgol_filter(y,5,1)
    fx[-1] = fx[0]
    fy[-1] = fy[0]

    # now do optional interpolation for final coordinates
    tck, u = interpolate.splprep([fx, fy],s=3)
    unew = np.arange(0, 1.01, 0.01)
    out = interpolate.splev(unew, tck)

    return out

if __name__ == '__main__':
    print("Running unit test...")

    from Bio.PDB import *
    parser = PDBParser()
    structure = parser.get_structure('4NOB', '4NOB.pdb')
    chain = structure[0]['A']
    coords = np.array([list(atom.get_vector()) for atom in chain.get_atoms()])
    pts = coords[:,:2]

    # center coordinates
    pts = pts - np.mean(pts, axis=0)

    out = representative_shape(pts)

    # output points
    np.savetxt('data/test_shape.csv', np.array(out).T, delimiter=',')

    # superimpose original points and representative shape
    plt.scatter(pts[:,0],pts[:,1])
    plt.plot(out[0],out[1],color='k')
    plt.show()
