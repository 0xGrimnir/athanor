"""Base-pair and base-pair-step geometry (the Cambridge convention, as computed by 3DNA).

Every base carries a reference frame (Olson et al., J. Mol. Biol. 313, 229, 2001): an origin and three
orthonormal axes, stored as a 3x3 matrix whose columns are x, y, z. Two frames are related by six
parameters - three translations and three rotations - defined by the El Hassan & Calladine (1995)
construction that 3DNA uses (Lu & Olson, Nucleic Acids Res. 31, 5108, 2003):

* the hinge is the axis about which the first z axis turns into the second (z1 x z2); the bending angle
  Gamma is the angle between them;
* each frame is rotated by half of Gamma about the hinge, in opposite senses, so the two z axes coincide;
  the middle frame takes that common z, and its y bisects the two rotated y axes;
* the angle between the rotated y axes about the common z is the twist (or, within a pair, the opening);
  Gamma splits into roll = Gamma cos(phi) and tilt = Gamma sin(phi), phi being the angle from the hinge
  to the middle y axis;
* the translations are the origin-to-origin vector in the middle frame.

For a step between base pairs the six numbers are shift, slide, rise (A) and tilt, roll, twist (degrees);
for the two bases of a pair, the same construction (with the second base's y and z reversed) gives
shear, stretch, stagger and buckle, propeller, opening.
"""
from __future__ import annotations

import numpy as np

FLIP = np.diag([1.0, -1.0, -1.0])       # strand II base frame -> its "reversed" frame


def rot_axis(axis: np.ndarray, deg: float) -> np.ndarray:
    """Right-handed rotation by ``deg`` about ``axis`` (Rodrigues)."""
    a = np.asarray(axis, dtype=float)
    a = a / np.linalg.norm(a)
    t = np.radians(deg)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * (K @ K)


def _signed_angle(a: np.ndarray, b: np.ndarray, axis: np.ndarray) -> float:
    """Angle from a to b about axis, in degrees (a and b projected onto the plane normal to axis)."""
    n = axis / np.linalg.norm(axis)
    a = a - np.dot(a, n) * n
    b = b - np.dot(b, n) * n
    return float(np.degrees(np.arctan2(np.dot(np.cross(a, b), n), np.dot(a, b))))


def step_parameters(R1: np.ndarray, o1: np.ndarray, R2: np.ndarray, o2: np.ndarray):
    """Six parameters relating frame 2 to frame 1, and the middle frame (R, origin).

    Returns ([shift, slide, rise, tilt, roll, twist], Rm, om)."""
    z1, z2 = R1[:, 2], R2[:, 2]
    hinge = np.cross(z1, z2)
    gamma = float(np.degrees(np.arccos(np.clip(np.dot(z1, z2), -1.0, 1.0))))
    if np.linalg.norm(hinge) < 1e-10:
        # z axes parallel: no bending; any hinge in the plane works - take the bisector of the y axes
        hinge = R1[:, 1] + R2[:, 1]
        if np.linalg.norm(hinge) < 1e-10:
            hinge = R1[:, 0]
        gamma = 0.0
    P2 = rot_axis(hinge, -gamma / 2) @ R2
    P1 = rot_axis(hinge, gamma / 2) @ R1
    zm = P1[:, 2] + P2[:, 2]
    zm /= np.linalg.norm(zm)
    ym = P1[:, 1] + P2[:, 1]
    ym /= np.linalg.norm(ym)
    xm = np.cross(ym, zm)
    Rm = np.column_stack([xm, ym, zm])
    om = (np.asarray(o1, float) + np.asarray(o2, float)) / 2
    twist = _signed_angle(P1[:, 1], P2[:, 1], zm)
    phi = np.radians(_signed_angle(hinge, ym, zm))
    roll = gamma * np.cos(phi)
    tilt = gamma * np.sin(phi)
    d = np.asarray(o2, float) - np.asarray(o1, float)
    shift, slide, rise = (float(np.dot(d, Rm[:, k])) for k in range(3))
    return np.array([shift, slide, rise, tilt, roll, twist]), Rm, om


def bz_check(R1, o1, R2, o2):
    """3DNA's treatment of left-handed and B-Z junction steps between two Watson-Crick pairs (analyze's
    default -bz): in Z-DNA the base-pair z axes point against the direction from one pair to the next,
    so both frames are turned 180 degrees about y (x and z reversed); at a B-Z junction, where the two z
    axes oppose each other, the frame pointing backwards is turned. Returns the (possibly) adjusted
    frames."""
    x1, y1, z1 = R1[:, 0], R1[:, 1], R1[:, 2]
    x2, y2, z2 = R2[:, 0], R2[:, 1], R2[:, 2]
    d = np.asarray(o2, float) - np.asarray(o1, float)
    turn = np.diag([-1.0, 1.0, -1.0])
    if np.dot(x1, x2) < 0 and np.dot(z1, z2) < 0 and np.dot(y1, y2) > 0:
        return (R1, R2 @ turn) if np.dot(d, z1) > 0 else (R1 @ turn, R2)
    if (np.dot(x1, x2) > 0 and np.dot(z1, z2) > 0 and np.dot(y1, y2) > 0 and np.dot(d, z1) < 0
            and np.dot(d, z2) < 0):
        return R1 @ turn, R2 @ turn
    return R1, R2


def step_transform(p) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """For parameters p = [shift, slide, rise, tilt, roll, twist], frame 2 in the coordinates of
    frame 1: (R, middle-frame rotation, displacement d) with R2 = R1 @ R, Rm = R1 @ Mrel and
    o2 = o1 + R1 @ Mrel @ [shift, slide, rise] (the inverse of :func:`step_parameters`)."""
    shift, slide, rise, tilt, roll, twist = (float(x) for x in p)
    gamma = float(np.hypot(roll, tilt))
    phi = float(np.degrees(np.arctan2(tilt, roll))) if gamma > 1e-12 else 0.0
    z = np.array([0.0, 0.0, 1.0])
    y = np.array([0.0, 1.0, 0.0])
    Rz = lambda a: rot_axis(z, a)  # noqa: E731
    R = Rz(twist / 2 - phi) @ rot_axis(y, gamma) @ Rz(twist / 2 + phi)
    M = Rz(twist / 2 - phi) @ rot_axis(y, gamma / 2) @ Rz(phi)
    return R, M, M @ np.array([shift, slide, rise])


def next_frame(R1, o1, p):
    """Frame 2 from frame 1 and the six step parameters."""
    R, _, d = step_transform(p)
    return R1 @ R, np.asarray(o1, float) + R1 @ d


def pair_parameters(Ri: np.ndarray, oi: np.ndarray, Rii: np.ndarray, oii: np.ndarray):
    """Base-pair parameters [shear, stretch, stagger, buckle, propeller, opening] for the strand I base
    frame (Ri, oi) and the strand II base frame (Rii, oii), and the base-pair frame (R, origin).
    As in 3DNA, the strand II frame is reversed (y and z negated) and the parameters relate it to the
    strand I frame."""
    p, Rm, om = step_parameters(Rii @ FLIP, oii, Ri, oi)
    return p, Rm, om


def pair_frames(Rm: np.ndarray, om: np.ndarray, p):
    """Base frames (strand I, strand II) for a base pair with frame (Rm, om) and parameters p:
    the inverse of :func:`pair_parameters`."""
    R, M, d = step_transform(p)
    Rx = Rm @ M.T                         # the reversed strand II frame
    ox = np.asarray(om, float) - 0.5 * (Rx @ d)
    Ri = Rx @ R
    oi = ox + Rx @ d
    return (Ri, oi), (Rx @ FLIP, ox)


def kabsch(P: np.ndarray, Q: np.ndarray):
    """Rotation R and translation t minimising |R @ P_k + t - Q_k| (P, Q: n x 3); and the RMSD."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    t = qc - R @ pc
    rmsd = float(np.sqrt(np.mean(np.sum((P @ R.T + t - Q) ** 2, axis=1))))
    return R, t, rmsd
