import numpy as np
from scipy.integrate import solve_ivp

# ————————————————————————————————
#  Constants
# ————————————————————————————————
_g0    = 9.80665         # sea-level gravity [m/s²]
_RE    = 6.371e6         # Earth radius [m]
_M     = 0.0289644       # molar mass of air [kg/mol]
_R0    = 8.314462618     # universal gas constant [J/(mol·K)]
_P0    = 101325          # sea-level pressure [Pa]

_layers = [
    (    0, 288.15, -0.0065),
    (11000, 216.65,  0.0),
    (20000, 216.65, +0.0010),
    (32000, 228.65, +0.0028),
    (47000, 270.65,  0.0),
    (51000, 270.65, -0.0028),
    (71000, 214.65, -0.0020)
]
_max_P_h = 100000  # altitude cap for pressure integration [m]


def gravity(h):
    """Return local gravity [m/s²] at altitude h [m]."""
    h_arr = np.asarray(h, float)
    return _g0 * (_RE / (_RE + h_arr))**2


def temperature(h):
    """Return standard-atmosphere temperature [K] at altitude h [m]."""
    def _T_scalar(h0):
        for hi, Ti, Li in reversed(_layers):
            if h0 >= hi:
                return Ti + Li*(h0 - hi)
        return _layers[0][1] + _layers[0][2]*h0

    h_arr = np.asarray(h, float)
    return np.vectorize(_T_scalar)(h_arr)


def pressure(h):
    """
    Return standard-atmosphere pressure [Pa] at altitude h [m].
    Solves dP/dh = –[M·g(h)/(R0·T(h))]·P via integrate.
    """
    h_arr = np.atleast_1d(h).astype(float)
    h_arr = np.clip(h_arr, 0, _max_P_h)

    sol = solve_ivp(
        lambda x, y: -y * (_M * gravity(x)) / (_R0 * temperature(x)),
        (0, float(h_arr.max())),
        [ _P0 ],
        t_eval=np.sort(h_arr),
        rtol=1e-6, atol=1e-6
    )
    P_full = sol.y[0]
    P_at_h = np.interp(h_arr, sol.t, P_full)

    return float(P_at_h) if np.isscalar(h) else P_at_h


def profile(h_array):
    """
    Given altitudes h_array [m], returns three arrays of the same shape:
      (gravity, temperature, pressure)
    """
    h_arr = np.asarray(h_array, float)
    return gravity(h_arr), temperature(h_arr), pressure(h_arr)


def grid_profile(n, h_min=0.0, h_max=84000.0):
    """
    Generate an altitude grid of n points from h_min to h_max,
    and return a dict with arrays:
      'h': altitudes,
      'g': gravity,
      'T': temperature,
      'P': pressure
    """
    h_arr = np.linspace(h_min, h_max, n)
    g_arr = gravity(h_arr)
    T_arr = temperature(h_arr)
    P_arr = pressure(h_arr)
    return {'h': h_arr, 'g': g_arr, 'T': T_arr, 'P': P_arr }

__all__ = ['gravity', 'temperature', 'pressure', 'profile', 'grid_profile']
