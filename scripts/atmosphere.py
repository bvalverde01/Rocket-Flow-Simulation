import numpy as np
from scipy.integrate import solve_ivp

_g0 = 9.80665
_RE = 6.371e6
_M  = 0.0289644
_R0 = 8.314462618
_P0 = 101325

_layers = [
    (    0, 288.15, -0.0065),
    (11000, 216.65,  0.0),
    (20000, 216.65, +0.0010),
    (32000, 228.65, +0.0028),
    (47000, 270.65,  0.0),
    (51000, 270.65, -0.0028),
    (71000, 214.65, -0.0020)
]
_max_P_h = 100000


def g(h):

    h_arr = np.asarray(h, dtype=float)
    return _g0 * (_RE/( _RE + h_arr))**2


def T(h):

    def _T_scalar(h0):
        for hi, Ti, Li in reversed(_layers):
            if h0 >= hi:
                return Ti + Li*(h0 - hi)
        # below lowest layer
        return _layers[0][1] + _layers[0][2]*h0

    h_arr = np.asarray(h, dtype=float)
    vec = np.vectorize(_T_scalar)
    return vec(h_arr)


def P(h):
    
    h_arr = np.atleast_1d(h).astype(float)
    # limit to [0, _max_P_h]
    h_arr = np.clip(h_arr, 0, _max_P_h)
    # solve ODE on the full range needed
    sol = solve_ivp(
        lambda x, y: -y*( _M * g(x)) / (_R0 * T(x)),
        (0, max(h_arr)),
        [_P0],
        t_eval=np.sort(h_arr),
        rtol=1e-6,
        atol=1e-6
    )
    sorted_h = sol.t
    sorted_P = sol.y[0]
    P_interp = np.interp(h_arr, sorted_h, sorted_P)
    if np.isscalar(h):
        return float(P_interp)
    return P_interp

__all__ = ['g', 'T', 'P']
