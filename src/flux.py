"""Sand flux roses and predicted dune crest orientations.

Two orientation modes (Courrech du Pont, Narteau and Gao 2014, Geology):

* bed instability (BI): crests take the orientation that maximises the gross
  sand flux normal to the crest, with the flux on the crest boosted by a
  speed-up factor that grows with the angle between wind and crest;
* fingering / elongation (F): dunes grow along the direction of the resultant
  of the boosted crest fluxes, so the crest is parallel to that resultant.

All directions are compass azimuths. Wind and flux directions point the way
the sand moves (downwind). Crest orientations are axial, in [0, 180).
"""
import numpy as np

KAPPA = 0.4
Z0 = 1e-3            # aerodynamic roughness of a sand bed (m)
G = 9.81
RHO_F, RHO_S = 1.225, 2650.0
D = 300e-6           # grain size (m)
GAMMA = 1.6          # crest speed-up factor (Courrech du Pont et al. 2014)
NBINS = 180          # 2 deg flux-rose bins
ALPHA = np.arange(0.0, 360.0, 0.5)

U_STAR_T = 0.1 * np.sqrt((RHO_S / RHO_F) * G * D)  # ~0.25 m/s, as in Gunn et al. 2022


def sand_flux(u, v, u_star_t=U_STAR_T):
    """Saturated flux magnitude (m^2/s) and downwind azimuth from 10 m wind components."""
    speed = np.hypot(u, v)
    u_star = speed * KAPPA / np.log(10.0 / Z0)
    q = 25.0 * (RHO_F / RHO_S) * np.sqrt(D / G) * np.clip(u_star**2 - u_star_t**2, 0.0, None)
    az = np.degrees(np.arctan2(u, v)) % 360.0
    return q, az


def flux_rose(u, v, axis=-1, **kw):
    """Mean flux per direction bin, shape (..., NBINS). Bin centres from rose_centres()."""
    q, az = sand_flux(u, v, **kw)
    idx = np.minimum((az / (360.0 / NBINS)).astype(int), NBINS - 1)
    q = np.moveaxis(q, axis, -1)
    idx = np.moveaxis(idx, axis, -1)
    flat_q = q.reshape(-1, q.shape[-1])
    flat_i = idx.reshape(-1, idx.shape[-1])
    out = np.zeros((flat_q.shape[0], NBINS))
    for r in range(flat_q.shape[0]):
        out[r] = np.bincount(flat_i[r], weights=flat_q[r], minlength=NBINS)
    return (out / q.shape[-1]).reshape(q.shape[:-1] + (NBINS,))


def rose_centres():
    return (np.arange(NBINS) + 0.5) * (360.0 / NBINS)


def resultant(rose):
    th = np.radians(rose_centres())
    e = (rose * np.sin(th)).sum(-1)
    n = (rose * np.cos(th)).sum(-1)
    rdp = np.hypot(e, n)
    dp = rose.sum(-1)
    return np.degrees(np.arctan2(e, n)) % 360.0, rdp, dp


def _crest_terms(rose, gamma):
    """sin and boosted flux of every rose bin against every candidate crest azimuth."""
    th = np.radians(rose_centres())[None, :]
    a = np.radians(ALPHA)[:, None]
    s = np.sin(th - a)                                   # (n_alpha, NBINS)
    boost = 1.0 + gamma * np.abs(s)
    return s, np.cos(th - a), boost


def orientation_bi(rose, gamma=GAMMA):
    """Crest azimuth (axial) maximising boosted gross bedform-normal transport."""
    s, _, boost = _crest_terms(rose, gamma)
    gbnt = rose @ (boost * np.abs(s)).T                     # (..., n_alpha)
    return ALPHA[np.argmax(gbnt, axis=-1)] % 180.0


def orientation_f(rose, gamma=GAMMA):
    """Elongation (fingering) crest azimuth, axial.

    The elongation direction alpha is where the boosted crest fluxes have no net
    component across the crest and a positive component along it. Among the
    roots we take the one with the largest along-crest flux.
    """
    s, c, boost = _crest_terms(rose, gamma)
    normal = rose @ (boost * s).T
    along = rose @ (boost * c).T
    sign = np.sign(normal)
    cross = (sign != np.roll(sign, -1, axis=-1)) & (along > 0)
    score = np.where(cross, along, -np.inf)
    j = np.argmax(score, axis=-1)
    has = np.isfinite(np.take_along_axis(score, j[..., None], -1)[..., 0])
    out = ALPHA[j] % 180.0
    return np.where(has, out, np.nan)


def axial_diff(a, b):
    d = np.abs(np.asarray(a) - np.asarray(b)) % 180.0
    return np.minimum(d, 180.0 - d)


def signed_axial(a, b):
    """Signed axial difference a - b in (-90, 90]."""
    return (np.asarray(a) - np.asarray(b) + 90.0) % 180.0 - 90.0
