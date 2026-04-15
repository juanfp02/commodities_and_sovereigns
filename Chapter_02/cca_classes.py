from scipy.optimize import least_squares, fsolve
from scipy.special import gammaln
from scipy.stats import norm
import numpy as np

###############################################################################################
# Common functions
###############################################################################################

def compute_lcl_usd(M_bn, Bd_bn, fx_rate, r_d, r_f, T=1.0):
    """
    Local-currency liabilities in USD.
    
    LCL$ = (M·e^{r_d·T} + B_d) · e^{-r_f·T} / X_F
    
    All inputs in billions of local currency. fx_rate = LC per USD.
    r_d, r_f in decimal. Returns USD billions.
    """
    lcl_lc = M_bn * np.exp(r_d * T) + Bd_bn
    return lcl_lc * np.exp(-r_f * T) / fx_rate


def compute_barrier_kvm(debt, r_f, T=1.0):
    """
    KMV-style distress barrier: B_f = ST + 0.5·LT + interest.
    All in USD millions.
    """
    interest = (debt) * r_f * T
    return (debt + interest)

###############################################################################################
# M0/M1/M2 Solvers 
###############################################################################################

class BaselineCCAPricer:
    """
    M0: Standard sovereign CCA baseline.
    Geometric Brownian motion, no jumps, no convenience yield.
    """

    def _equations(self, V, sigma_V, LCL_usd, sigma_lcl, B_f, r_f, T):
        sqt = np.sqrt(T)
        d1  = (np.log(V / B_f) + (r_f + 0.5 * sigma_V**2) * T) / (sigma_V * sqt)
        d2  = d1 - sigma_V * sqt

        eq1 = V * norm.cdf(d1) - B_f * np.exp(-r_f * T) * norm.cdf(d2) - LCL_usd
        eq2 = V * sigma_V * norm.cdf(d1) - LCL_usd * sigma_lcl

        return np.array([eq1, eq2])

    def CCA_system_M0(self, log_unknowns, LCL_usd, sigma_lcl, B_f, r_f, T):
        V       = np.exp(log_unknowns[0])
        sigma_V = np.exp(log_unknowns[1])
        return self._equations(V, sigma_V, LCL_usd, sigma_lcl, B_f, r_f, T)

    def _guesses(self, LCL_usd, sigma_lcl, B_f):
        V_base   = LCL_usd + B_f
        sig_base = sigma_lcl * LCL_usd / V_base

        return [
            (V_base,        sig_base),
            (V_base * 1.5,  sig_base),
            (V_base * 0.5,  sig_base),
            (V_base,        sig_base * 2.0),
            (V_base,        sig_base * 0.5),
            (LCL_usd * 2,   sig_base),
            (B_f * 1.5,     sig_base),
        ]

    def solve_CCA_M0(self, LCL_usd, sigma_lcl, B_f, r_f, T):
        if any(np.isnan(x) or x <= 0 for x in [LCL_usd, sigma_lcl, B_f]):
            return {'implied_V': np.nan, 'implied_sigma_V': np.nan, 'converged': False}

        best       = {'implied_V': np.nan, 'implied_sigma_V': np.nan, 'converged': False}
        best_resid = np.inf

        for V0, sig0 in self._guesses(LCL_usd, sigma_lcl, B_f):
            if V0 <= 0 or sig0 <= 0:
                continue

            # --- Method 1: fsolve on log-transformed unknowns ---
            try:
                sol, info, ier, _ = fsolve(
                    self.CCA_system_M0,
                    x0=[np.log(V0), np.log(sig0)],
                    args=(LCL_usd, sigma_lcl, B_f, r_f, T),
                    full_output=True
                )
                V, sig = np.exp(sol[0]), np.exp(sol[1])
                resid  = np.sum(info['fvec']**2)
                if ier == 1 and V > 0 and sig > 0 and resid < best_resid:
                    best       = {'implied_V': V, 'implied_sigma_V': sig, 'converged': True}
                    best_resid = resid
                    if resid < 1e-12:
                        return best
            except Exception:
                pass

            # --- Method 2: least_squares with bounds ---
            try:
                sol = least_squares(
                    lambda x: self._equations(x[0], x[1], LCL_usd, sigma_lcl, B_f, r_f, T),
                    x0=[V0, sig0],
                    bounds=([LCL_usd * 0.1, 1e-4], [V0 * 20, 5.0]),
                    method='trf',
                    xtol=1e-10, ftol=1e-10,
                    max_nfev=5000
                )
                V, sig = sol.x
                resid  = np.sum(sol.fun**2)
                if sol.success and V > 0 and sig > 0 and resid < best_resid:
                    best       = {'implied_V': V, 'implied_sigma_V': sig, 'converged': True}
                    best_resid = resid
                    if resid < 1e-12:
                        return best
            except Exception:
                pass

        if best_resid < 1e-6 and not best['converged']:
            best['converged'] = True

        return best

class ConvenienceYieldCCAPricer:
    """
    M1: Convenience yield adjustment to sovereign asset drift.
    Solves for implied V and sigma_V given the oil futures convenience yield.
    """

    def __init__(self, gamma):
        self.gamma = gamma

    def _sigma_total(self, sigma_V, sigma_y):
        return np.sqrt(sigma_V**2 + (self.gamma * sigma_y)**2)

    def _equations(self, V, sigma_V, LCL_usd, sigma_lcl, B_f, r_f, y, sigma_y, T):
        gy          = self.gamma * y
        gs          = self.gamma * sigma_y
        sigma_total = np.sqrt(sigma_V**2 + gs**2)
        sqt         = np.sqrt(T)

        # Geske (1978) stochastic dividend discount
        phi_inv = np.exp(-gy * T + 0.5 * gs**2 * T)
        Veff    = V * phi_inv

        d1 = (np.log(Veff / B_f) + (r_f + 0.5 * sigma_total**2) * T) / (sigma_total * sqt)
        d2 = d1 - sigma_total * sqt

        eq1 = Veff * norm.cdf(d1) - B_f * np.exp(-r_f * T) * norm.cdf(d2) - LCL_usd
        eq2 = Veff * sigma_total * norm.cdf(d1) - LCL_usd * sigma_lcl

        return np.array([eq1, eq2])

    def CCA_system_M1(self, log_unknowns, LCL_usd, sigma_lcl, B_f, r_f, y, sigma_y, T):
        V       = np.exp(log_unknowns[0])
        sigma_V = np.exp(log_unknowns[1])
        return self._equations(V, sigma_V, LCL_usd, sigma_lcl, B_f, r_f, y, sigma_y, T)

    def _guesses(self, LCL_usd, sigma_lcl, B_f, y, T):
        gy      = self.gamma * y
        V_base  = LCL_usd + B_f
        sig_base = sigma_lcl * LCL_usd / V_base
        V_cy    = V_base * np.exp(gy * T)
        sig_cy  = sig_base * np.exp(abs(gy) * T)   # fixed: no longer redundant

        return [
            (V_cy,        sig_cy),
            (V_base,      sig_base),
            (V_cy  * 1.5, sig_cy),
            (V_cy  * 0.5, sig_cy),
            (V_base,      sig_base * 2.0),
            (V_base,      sig_base * 0.5),
            (LCL_usd * 2, sig_base),
            (B_f   * 1.5, sig_base),
        ]

    def solve_CCA_M1(self, LCL_usd, sigma_lcl, B_f, r_f, y, sigma_y, T):
            if any(np.isnan(x) or x <= 0 for x in [LCL_usd, sigma_lcl, B_f]):
                return {'implied_V': np.nan, 'implied_sigma_V': np.nan, 'converged': False, 'convenience_yield': np.nan}
            if np.isnan(y):
                return {'implied_V': np.nan, 'implied_sigma_V': np.nan, 'converged': False, 'convenience_yield': np.nan}

            best        = {'implied_V': np.nan, 'implied_sigma_V': np.nan, 'converged': False, 'convenience_yield': y}
            best_resid  = np.inf

            for V0, sig0 in self._guesses(LCL_usd, sigma_lcl, B_f, y, T):
                if V0 <= 0 or sig0 <= 0:
                    continue

                try:
                    sol, info, ier, _ = fsolve(
                        self.CCA_system_M1,
                        x0=[np.log(V0), np.log(sig0)],
                        args=(LCL_usd, sigma_lcl, B_f, r_f, y, sigma_y, T),
                        full_output=True
                    )
                    V, sig = np.exp(sol[0]), np.exp(sol[1])
                    resid  = np.sum(info['fvec']**2)
                    if ier == 1 and V > 0 and sig > 0 and resid < best_resid:
                        best       = {'implied_V': V,
                                    'implied_sigma_V': self._sigma_total(sig, sigma_y),
                                    'converged': True,
                                    'convenience_yield': y}
                        best_resid = resid
                        if resid < 1e-12:
                            return best
                except Exception:
                    pass

                try:
                    def resid_func(x):
                        V, sig = x
                        if V <= 0 or sig <= 0:
                            return np.array([1e10, 1e10])
                        return self._equations(V, sig, LCL_usd, sigma_lcl,
                                            B_f, r_f, y, sigma_y, T)

                    sol   = least_squares(
                        resid_func,
                        x0=[V0, sig0],
                        bounds=([LCL_usd * 0.1, 1e-4], [V0 * 20, 5.0]),
                        method='trf',
                        xtol=1e-10, ftol=1e-10,
                        max_nfev=5000
                    )
                    V, sig = sol.x
                    resid  = np.sum(sol.fun**2)
                    if sol.success and V > 0 and sig > 0 and resid < best_resid:
                        best       = {'implied_V': V,
                                    'implied_sigma_V': self._sigma_total(sig, sigma_y),
                                    'converged': True,
                                    'convenience_yield': y}
                        best_resid = resid
                        if resid < 1e-12:
                            return best
                except Exception:
                    pass

            if best_resid < 1e-6 and not best['converged']:
                best['converged'] = True

            return best

class FixedJumpCCAPricer:
    """
    M2: Fixed jump size J, OVX-driven intensity λ.
    Merton (1976) series with corrected volatility transfer equation.
    """

    def __init__(self, J, max_terms=150):
        self.J         = J
        self.max_terms = max_terms
        self.log1J     = np.log(1 + J)
        self.k         = J
        self._cache    = {}   # missing

    def sigma_total(self, sigma_diff, lam):
        return np.sqrt(sigma_diff**2 + lam * self.log1J**2)

    def _series(self, V, B, r, T, sigma_diff, lam):
        lam_prime = lam * (1 + self.k)
        lam_T     = lam_prime * T
        sqt       = np.sqrt(T)

        sigma_pois = np.sqrt(lam_T)
        n_low  = max(0, int(lam_T - 6 * sigma_pois))
        n_high = int(lam_T + 6 * sigma_pois) + 1
        ns     = np.arange(n_low, n_high + 1)

        log_w   = -lam_T + ns * np.log(lam_T + 1e-300) - gammaln(ns + 1)
        weights = np.exp(log_w)
        weights[weights < 1e-15] = 0.0

        r_ns = r - lam * self.k + ns * self.log1J / T
        d1s  = (np.log(V / B) + (r_ns + 0.5 * sigma_diff**2) * T) / (sigma_diff * sqt)
        d2s  = d1s - sigma_diff * sqt

        calls  = V * norm.cdf(d1s) - B * np.exp(-r_ns * T) * norm.cdf(d2s)
        deltas = norm.cdf(d1s)

        call  = float(weights @ calls)
        delta = float(weights @ deltas)

        return call, delta
    
    def _series_cached(self, V, B, r, T, sigma_diff, lam):
        key = (round(V, 8), round(sigma_diff, 8), round(lam, 6))
        if key not in self._cache:
            self._cache[key] = self._series(V, B, r, T, sigma_diff, lam)
        return self._cache[key]

    def CCA_system_M2(self, unknowns, LCL_usd, sigma_lcl, B_f, r_f, T, lam):
        V, sigma_diff = unknowns
        if V <= 0 or sigma_diff <= 0:
            return np.array([1e10, 1e10])

        call,   delta = self._series_cached(V,                B_f, r_f, T, sigma_diff, lam)
        call_J, _     = self._series_cached(V * (1 + self.k), B_f, r_f, T, sigma_diff, lam)

        if len(self._cache) > 500:
            self._cache.clear()

        diffusion_term = sigma_diff * V * delta
        jump_term      = np.sqrt(lam) * (call_J - call)

        eq1 = call - LCL_usd
        eq2 = sigma_diff * V * delta + lam * (call_J - call) - LCL_usd * sigma_lcl
        
        return np.array([eq1, eq2])
    

    def solve_CCA_M2(self, LCL_usd, sigma_lcl, B_f, r_f, T, lam,
                    v_guess=None, sig_guess=None, prev_solution=None):

        if any(np.isnan(x) or x <= 0 for x in [LCL_usd, sigma_lcl, B_f]):
            return {'implied_V': np.nan, 'implied_sigma_diff': np.nan,
                    'implied_sigma_V': np.nan, 'converged': False}

        V0   = v_guess   if v_guess   else LCL_usd + B_f
        sig0 = sig_guess if sig_guess else sigma_lcl * LCL_usd / V0

        guesses = [
            (V0,         sig0),
            (V0 * 2.0,   sig0),
            (V0 * 0.5,   sig0),
            (V0,         sig0 * 2.0),
            (V0,         sig0 * 0.5),
            (B_f * 2.0,  sig0),
            (B_f * 3.0,  sig0),
            (B_f * 5.0,  sig0 * 0.5),
            (B_f * 10.0, sig0 * 0.3),
        ]

        # previous solution goes first — best possible warm start
        if prev_solution is not None:
            prev_V   = prev_solution.get('implied_V')
            prev_sig = prev_solution.get('implied_sigma_diff')
            if prev_V and prev_sig and prev_V > 0 and prev_sig > 0:
                guesses.insert(0, (prev_V, prev_sig))

        best_resid = np.inf
        best       = None

        for V_init, sig_init in guesses:
            if V_init <= 0 or sig_init <= 0:
                continue
            try:
                sol = least_squares(
                    lambda x: self.CCA_system_M2(x, LCL_usd, sigma_lcl, B_f, r_f, T, lam),
                    x0=[V_init, sig_init],
                    bounds=([1e-6, 1e-6], [np.inf, np.inf]),
                    method='trf',
                    xtol=1e-10, ftol=1e-10, gtol=1e-10,
                    max_nfev=10000,
                    x_scale=[V_init, sig_init]
                )
                resid = np.sum(sol.fun**2)
                if resid < best_resid and sol.x[0] > 0 and sol.x[1] > 0:
                    best_resid = resid
                    best       = sol.x
                if best_resid < 1e-12:
                    break
            except Exception:
                continue

        if best is not None and best_resid < 1.0:
            V, sigma_diff = best
            return {
                'implied_V'        : V,
                'implied_sigma_diff': sigma_diff,
                'implied_sigma_V'  : self.sigma_total(sigma_diff, lam),
                'converged'        : True
            }

        return {'implied_V': np.nan, 'implied_sigma_diff': np.nan,
                'implied_sigma_V': np.nan, 'converged': False}
    

class LogNormalJumpCCAPricer:
    """
    M2: Log-normal jump size, OVX-driven intensity λ.
    Merton (1976) closed-form series with ln(1+J) ~ N(μ_J, σ_J²).
    """

    def __init__(self, mu_J, sigma_J, max_terms=150):
        self.mu_J      = mu_J
        self.sigma_J   = sigma_J
        self.k         = np.exp(mu_J + 0.5 * sigma_J**2) - 1  # E[J] risk-neutral
        self._cache    = {}

    def sigma_total(self, sigma_diff, lam):
        return np.sqrt(sigma_diff**2 + lam * (self.mu_J**2 + self.sigma_J**2))

    def _series(self, V, B, r, T, sigma_diff, lam):
        lam_prime = lam * (1 + self.k)
        lam_T     = lam_prime * T
        sqt       = np.sqrt(T)

        sigma_pois = np.sqrt(lam_T) if lam_T > 0 else 1.0
        n_low  = max(0, int(lam_T - 6 * sigma_pois))
        n_high = int(lam_T + 6 * sigma_pois) + 1
        ns     = np.arange(n_low, n_high + 1)

        log_w   = -lam_T + ns * np.log(lam_T + 1e-300) - gammaln(ns + 1)
        weights = np.exp(log_w)
        weights[weights < 1e-15] = 0.0

        # per-term parameters — key difference from fixed-J
        r_ns     = r - lam * self.k + ns * (self.mu_J + 0.5 * self.sigma_J**2) / T
        sigma_ns = np.sqrt(sigma_diff**2 + ns * self.sigma_J**2 / T)

        d1s = (np.log(V / B) + (r_ns + 0.5 * sigma_ns**2) * T) / (sigma_ns * sqt)
        d2s = d1s - sigma_ns * sqt

        calls  = V * norm.cdf(d1s) - B * np.exp(-r_ns * T) * norm.cdf(d2s)
        deltas = norm.cdf(d1s)

        call  = float(weights @ calls)
        delta = float(weights @ deltas)

        return call, delta

    def _series_cached(self, V, B, r, T, sigma_diff, lam):
        key = (round(V, 8), round(sigma_diff, 8), round(lam, 6))
        if key not in self._cache:
            self._cache[key] = self._series(V, B, r, T, sigma_diff, lam)
        return self._cache[key]

    def CCA_system_M2(self, unknowns, LCL_usd, sigma_lcl, B_f, r_f, T, lam):
        V, sigma_diff = unknowns
        if V <= 0 or sigma_diff <= 0:
            return np.array([1e10, 1e10])

        call,  delta = self._series_cached(V, B_f, r_f, T, sigma_diff, lam)

        if len(self._cache) > 500:
            self._cache.clear()

        diffusion_term = sigma_diff * V * delta
        jump_term      = lam * self.k * V * delta   # linear approximation

        eq1 = call - LCL_usd
        eq2 = diffusion_term + jump_term - LCL_usd * sigma_lcl

        return np.array([eq1, eq2])

    def solve_CCA_M2(self, LCL_usd, sigma_lcl, B_f, r_f, T, lam,
                     v_guess=None, sig_guess=None, prev_solution=None):

        if any(np.isnan(x) or x <= 0 for x in [LCL_usd, sigma_lcl, B_f]):
            return {'implied_V': np.nan, 'implied_sigma_diff': np.nan,
                    'implied_sigma_V': np.nan, 'converged': False}

        V0   = v_guess   if v_guess   else LCL_usd + B_f
        sig0 = sig_guess if sig_guess else sigma_lcl * LCL_usd / V0

        guesses = [
            (V0,         sig0),
            (V0 * 2.0,   sig0),
            (V0 * 0.5,   sig0),
            (V0,         sig0 * 2.0),
            (V0,         sig0 * 0.5),
            (B_f * 1.1,  sig0),
            (B_f * 1.2,  sig0),
            (B_f * 2.0,  sig0),
            (B_f * 3.0,  sig0),
            (B_f * 5.0,  sig0 * 0.5),
            (B_f * 10.0, sig0 * 0.3),
        ]

        if prev_solution is not None:
            prev_V   = prev_solution.get('implied_V')
            prev_sig = prev_solution.get('implied_sigma_diff')
            if prev_V and prev_sig and prev_V > 0 and prev_sig > 0:
                guesses.insert(0, (prev_V, prev_sig))

        best_resid = np.inf
        best       = None

        for V_init, sig_init in guesses:
            if V_init <= 0 or sig_init <= 0:
                continue
            try:
                sol = least_squares(
                    lambda x: self.CCA_system_M2(
                        x, LCL_usd, sigma_lcl, B_f, r_f, T, lam),
                    x0=[V_init, sig_init],
                    bounds=([1e-6, 1e-6], [np.inf, np.inf]),
                    method='trf',
                    xtol=1e-10, ftol=1e-10, gtol=1e-10,
                    max_nfev=10000,
                    x_scale=[V_init, sig_init]
                )
                resid = np.sum(sol.fun**2)
                if resid < best_resid and sol.x[0] > 0 and sol.x[1] > 0:
                    best_resid = resid
                    best       = sol.x
                if best_resid < 1e-12:
                    break
            except Exception:
                continue

        if best is not None and best_resid < 1.0:
            V, sigma_diff = best
            return {
                'implied_V'        : V,
                'implied_sigma_diff': sigma_diff,
                'implied_sigma_V'  : self.sigma_total(sigma_diff, lam),
                'converged'        : True
            }

        return {'implied_V': np.nan, 'implied_sigma_diff': np.nan,
                'implied_sigma_V': np.nan, 'converged': False}