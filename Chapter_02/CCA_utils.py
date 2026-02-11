from pyexpat import model
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import fsolve, brentq
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.optimize import curve_fit
import statsmodels.api as sm
from scipy.special import gammaln
from numba import jit, prange, float64, int64
import numba as nb


#Core functions
def CCA_system(unknowns, LCL_usd, sigma_lcl, B_f, r_f,T):

    """
    Returns the values of the two CCA equations:

    Eq1: Call pricing: LCL_usd = V N(d1) - B_f e^{-rT} N(d2) 
    Eq2: Volatility transfer: LCL_usd·sigma_LCL = V sigma_V N(d1)

    """

    V, sigma_V = unknowns

    if V<=0 or sigma_V<=0:
        return (1e10, 1e10)
    
    d1 = (np.log(V / B_f) + (r_f + 0.5 * sigma_V**2) * T) / (sigma_V * np.sqrt(T))
    d2 = d1 - sigma_V * np.sqrt(T)
    
    eq1 = V * norm.cdf(d1) - B_f * np.exp(-r_f * T) * norm.cdf(d2) - LCL_usd
    eq2 = V * sigma_V * norm.cdf(d1) - LCL_usd * sigma_lcl
    
    return (eq1, eq2)

def solve_CCA(LCL_usd, sigma_lcl, B_f, r_f, T):
    """
    Solve for implied sovereign asset value (V) and volatility (sigma_V).
    
    Returns dict: {'V', 'sigma_V', 'converged'}
    """
    if any(np.isnan(x) or x <= 0 for x in [LCL_usd, sigma_lcl, B_f]):
        return {'V': np.nan, 'sigma_V': np.nan, 'converged': False}
    
    # Initial guess: V ≈ LCL + B_f, σ_V ≈ de-levered σ_LCL
    V0 = LCL_usd + B_f
    sig0 = sigma_lcl * LCL_usd / V0
    
    try:
        sol, info, ier, msg = fsolve(
            CCA_system, x0=(V0, sig0),
            args=(LCL_usd, sigma_lcl, B_f, r_f, T),
            full_output=True
        )
        V, sig = sol
        ok = (ier == 1) and (V > 0) and (sig > 0)
        return {'V': V if ok else np.nan, 'sigma_V': sig if ok else np.nan, 'converged': ok}
    except Exception:
        return {'V': np.nan, 'sigma_V': np.nan, 'converged': False}

def compute_risk(V, sigma_V, B_f, r_f, T=1.0):
    """
    Compute risk indicators from solved CCA.
    
    Returns dict: {'d2', 'default_prob', 'credit_spread_bps', 
                   'put_value', 'risky_debt', 'leverage'}
    """
    nans = {'d2': np.nan, 'default_prob': np.nan, 'credit_spread_bps': np.nan,
            'put_value': np.nan, 'risky_debt': np.nan, 'leverage': np.nan}
    
    if any(np.isnan(x) for x in [V, sigma_V, B_f]) or V <= 0 or sigma_V <= 0:
        return nans
    
    d1 = (np.log(V / B_f) + (r_f + 0.5 * sigma_V**2) * T) / (sigma_V * np.sqrt(T))
    d2 = d1 - sigma_V * np.sqrt(T)
    
    default_prob = norm.cdf(-d2)
    
    # Put = implicit guarantee against default
    put = B_f * np.exp(-r_f * T) * norm.cdf(-d2) - V * norm.cdf(-d1)
    put = max(put, 0.0)
    
    # Risky debt = default-free debt - put
    df_debt = B_f * np.exp(-r_f * T)
    risky_debt = df_debt - put
    
    # Credit spread: D = B·e^{-yT} → y = ln(B/D)/T → spread = y - r
    if risky_debt > 0:
        y = np.log(B_f / risky_debt) / T
        spread_bps = (y - r_f) * 10000
    else:
        spread_bps = np.nan
    
    return {
        'd2': d2,
        'default_prob': default_prob,
        'credit_spread_bps': spread_bps,
        'put_value': put,
        'risky_debt': risky_debt,
        'leverage': B_f / V
    }


# ============================================================
# OBSERVABLE CONSTRUCTION
# ============================================================

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


def compute_fx_volatility(fx_series, window=52):
    """
    Rolling annualized FX volatility from weekly log returns.
    Input: pd.Series of FX rates (LC per USD).
    Returns: pd.Series of annualized vol.
    """
    log_ret = np.log(fx_series / fx_series.shift(1))
    return log_ret.rolling(window=window).std() * np.sqrt(52)

########################################################
# Model 1: LCL augmentation and eta estimation
########################################################
def augment_lcl(LCL_standard, P_oil, P_oil_ref, theta, eta):

    if P_oil <= 0 or P_oil_ref <= 0:
        return LCL_standard, 1.0

    multiplier = theta
    oil_return = (P_oil / P_oil_ref - 1)

    return LCL_standard * (1 + oil_return)**multiplier, multiplier


def estimate_eta(r_lcl, r_oil):

    reg_df = pd.DataFrame({
        'r_lcl': r_lcl,
        'r_oil': r_oil
    }).dropna()

    X = sm.add_constant(reg_df['r_oil'])   # adds α (intercept)
    y = reg_df['r_lcl']

    model = sm.OLS(y, X).fit()

    return model.params['r_oil']

def pv_of_future_oil_revenue(government_take, prod_qty, reference_price, current_price, discount_rate, mean_reversion_spread):

    delta = current_price - reference_price
    discounting_factor = 1/(discount_rate+mean_reversion_spread)

    return government_take * prod_qty * delta * discounting_factor



########################################################
# Model 2: Adding jumps
########################################################

@jit(nopython=True, cache=True)
def montecarlo_jump_paths_numba(S, T, r, sigma, lam, eta, steps, n_paths):
    
    dt = T / steps
    drift_term = (r - 0.5 * sigma**2 + lam / (eta + 1.0)) * dt
    vol_term = sigma * np.sqrt(dt)

    paths = np.empty((steps, n_paths))

    for i in range(steps):
        normal_diffusion = np.random.normal(0.0, 1.0, n_paths)
        poisson_counts = np.random.poisson(lam * dt, n_paths)

        step_jump = np.zeros(n_paths)
        for j in range(n_paths):
            k = poisson_counts[j]
            if k > 0:
                step_jump[j] = -np.random.gamma(k, 1.0 / eta)  # sum of k exponentials, negated

        if i == 0:
            paths[i] = drift_term + vol_term * normal_diffusion + step_jump
        else:
            paths[i] = paths[i - 1] + drift_term + vol_term * normal_diffusion + step_jump

    return S * np.exp(paths)

@jit(nopython=True, cache=True)
def montecarlo_option_pricer_numba(paths, T, K, r_f, option_type):
    """Numba-optimized option pricing"""
    final_prices = paths[-1]
    
    if option_type == 'Call':
        payoffs = np.maximum(final_prices - K, 0.0)
    else:  # Put
        payoffs = np.maximum(K - final_prices, 0.0)
    
    price = np.mean(payoffs) * np.exp(-r_f * T)
    return price

@jit(nopython=True, cache=True)
def compute_black_scholes_jump(V, B_f, r_f, T, sigma_total):
    """Black-Scholes formulas for jump diffusion"""
    d1 = (np.log(V / B_f) + (r_f + 0.5 * sigma_total**2) * T) / (sigma_total * np.sqrt(T))
    d2 = d1 - sigma_total * np.sqrt(T)
    return d1, d2

from scipy.optimize import fsolve
from scipy.stats import norm
import numpy as np

def sigma_total_from_exp_jumps(sigma_diff, lam, eta):
    # instantaneous variance add-on from compound Poisson with Y=-Exp(eta)
    return np.sqrt(sigma_diff**2 + lam * (2.0 / (eta**2)))

class JumpDiffusionPricer:
    def __init__(self, steps=252, n_paths=50000):
        self.steps = steps
        self.n_paths = n_paths

    def montecarlo_jump_paths(self, S, T, r, sigma, lam, eta, seed=None):
        if seed is not None:
            np.random.seed(seed)
        return montecarlo_jump_paths_numba(S, T, r, sigma, lam, eta, self.steps, self.n_paths)

    def CCA_system_jd(self, unknowns, LCL_usd, sigma_lcl, B_f, r_f, T, lam, eta, seed):
        V, sigma_diff = unknowns

        if V <= 0 or sigma_diff <= 0:
            return np.array([1e10, 1e10])

        # Eq1: equity = call(V, B)
        paths = self.montecarlo_jump_paths(V, T, r_f, sigma_diff, lam, eta, seed)
        call_jd = montecarlo_option_pricer_numba(paths, T, B_f, r_f, 'Call')
        eq1 = call_jd - LCL_usd

        # Eq2: volatility matching (CCA-style)
        sigma_total = sigma_total_from_exp_jumps(sigma_diff, lam, eta)
        d1 = (np.log(V / B_f) + (r_f + 0.5 * sigma_total**2) * T) / (sigma_total * np.sqrt(T))
        eq2 = V * sigma_total * norm.cdf(d1) - LCL_usd * sigma_lcl

        return np.array([eq1, eq2])

    def solve_CCA_jd(self, LCL_usd, sigma_lcl, B_f, r_f, T, lam, eta, seed=42):
        if any(np.isnan(x) or x <= 0 for x in [LCL_usd, sigma_lcl, B_f]):
            return self._nan_result()

        V0 = LCL_usd + B_f
        sig0 = sigma_lcl * LCL_usd / V0

        try:
            sol = fsolve(
                lambda x: self.CCA_system_jd(x, LCL_usd, sigma_lcl, B_f, r_f, T, lam, eta, seed),
                x0=[V0, sig0],
                full_output=True
            )
            V, sig_diff = sol[0]

            if (sol[2] == 1) and (V > 0) and (sig_diff > 0):
                sigma_total = sigma_total_from_exp_jumps(sig_diff, lam, eta)
                return {'V': V, 'sigma_diff': sig_diff, 'sigma_total': sigma_total, 'converged': True}

            return self._nan_result()
        except:
            return self._nan_result()

    def compute_risk_jd(self, V, sigma_diff, B_f, r_f, T, lam, eta, seed=42):
        if any(np.isnan(x) for x in [V, sigma_diff, B_f]) or V <= 0 or sigma_diff <= 0:
            return self._nan_risk_result()

        try:
            paths = self.montecarlo_jump_paths(V, T, r_f, sigma_diff, lam, eta, seed)
            put_jd  = montecarlo_option_pricer_numba(paths, T, B_f, r_f, 'Put')

            sigma_total = sigma_total_from_exp_jumps(sigma_diff, lam, eta)
            d1 = (np.log(V / B_f) + (r_f + 0.5 * sigma_total**2) * T) / (sigma_total * np.sqrt(T))
            d2 = d1 - sigma_total * np.sqrt(T)

            df_debt = B_f * np.exp(-r_f * T)
            risky_debt = df_debt - put_jd

            spread = np.nan
            if risky_debt > 0:
                spread = (np.log(B_f / risky_debt) / T - r_f) * 10000

            return {
                'sigma_V': sigma_total, 
                'd2': d2,
                'default_prob': norm.cdf(-d2),
                'credit_spread_bps': spread,
                'put_value': put_jd,
                'risky_debt': risky_debt,
                'leverage': B_f / V
            }
        except:
            return self._nan_risk_result()

    @staticmethod
    def _nan_result():
        return {'V': np.nan, 'sigma_diff': np.nan, 'sigma_total': np.nan, 'converged': False}

    @staticmethod
    def _nan_risk_result():
        return {k: np.nan for k in ['d2', 'default_prob', 'credit_spread_bps',
                                    'put_value', 'risky_debt', 'leverage']}



@jit(nopython=True, parallel=True, cache=True)
def batch_montecarlo_jump_paths(S0, T, r, sigma, lam, eta, steps, n_paths, n_batches, K):
    results = np.empty((n_batches, 2))

    for i in prange(n_batches):
        paths = montecarlo_jump_paths_numba(S0, T, r, sigma, lam, eta, steps, n_paths)
        ST = paths[-1]

        call_price = np.mean(np.maximum(ST - K, 0.0)) * np.exp(-r * T)
        put_price  = np.mean(np.maximum(K - ST, 0.0)) * np.exp(-r * T)

        results[i, 0] = call_price
        results[i, 1] = put_price

    return results

########################################################
# Graphs
########################################################
def plot_d2_vs_cds(df, countries=None, ncols=3, save_path=None):
    if countries is None:
        countries = sorted(df['country'].unique())
    
    
    nrows = int(np.ceil(len(countries) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 3.5*nrows))
    axes = np.atleast_2d(axes).flatten()
    
    for i, country in enumerate(countries):
        ax1 = axes[i]
        ax2 = ax1.twinx()
        d = df[df['country'] == country].sort_values('date')
        d['date'] = pd.to_datetime(d['date'])
        
        l1, = ax1.plot(d['date'], d['distance_to_distress'], color='#2C3E50', lw=1, label='d2')
        l2, = ax2.plot(d['date'], d['cds_spread'], color='#E74C3C', lw=1, alpha=0.8, label='CDS')
        
        ax1.set_ylabel('d2', fontsize=8, color='#2C3E50')
        ax2.set_ylabel('CDS (bps)', fontsize=8, color='#E74C3C')
        ax1.set_title(country, fontsize=10, fontweight='bold')
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%y'))
        ax1.tick_params(labelsize=7)
        ax2.tick_params(labelsize=7)
        ax1.grid(True, alpha=0.2)
        
        # Invert d2 axis so lower d2 (more risk) aligns with higher CDS
        ax1.invert_yaxis()
        
        # Correlation
        valid = d.dropna(subset=['distance_to_distress', 'cds_spread'])
        if len(valid) > 10:
            corr = valid['distance_to_distress'].corr(valid['cds_spread'])
            ax1.text(0.03, 0.08, f'ρ={corr:.2f}', transform=ax1.transAxes, fontsize=7,
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    for j in range(len(countries), len(axes)):
        axes[j].set_visible(False)
    
    axes[0].legend([l1, l2], ['Distance-to-Distress (d2)', 'Market CDS'], fontsize=7, loc='upper left')
    fig.suptitle('Distance-to-Distress vs Market CDS Spreads', fontsize=13, fontweight='bold')
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def _power_func(x, a, b):
    return a * np.power(x, b)


def plot_d2_scatter(df, countries=None, ncols=3, save_path=None):
    if countries is None:
        countries = sorted(df['country'].unique())
    

    nrows = int(np.ceil(len(countries) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.atleast_2d(axes).flatten()

    for i, country in enumerate(countries):
        ax = axes[i]
        d = df[df['country'] == country].dropna(subset=['distance_to_distress', 'cds_spread'])
        d = d[(d['distance_to_distress'] > 0) & (d['cds_spread'] > 0)]

        if len(d) < 20:
            ax.set_title(f'{country} (insufficient data)', fontsize=10)
            ax.set_visible(False)
            continue

        x = d['distance_to_distress'].values
        y = d['cds_spread'].values

        ax.scatter(x, y, s=10, alpha=0.5, color='#2C5F8A', edgecolors='none')

        # Fit power curve: spread = a * d2^b
        try:
            popt, _ = curve_fit(_power_func, x, y, p0=[500, -2], maxfev=5000)
            x_fit = np.linspace(max(x.min(), 0.01), x.max(), 200)
            y_fit = _power_func(x_fit, *popt)
            ax.plot(x_fit, y_fit, 'k-', lw=2)
        except Exception:
            pass

        # R² from log-log OLS
        lx, ly = np.log(x), np.log(y)
        mask = np.isfinite(lx) & np.isfinite(ly)
        if mask.sum() > 10:
            ss_res = np.sum((ly[mask] - np.polyval(np.polyfit(lx[mask], ly[mask], 1), lx[mask])) ** 2)
            ss_tot = np.sum((ly[mask] - ly[mask].mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
            ax.text(0.95, 0.92, f'R² = {r2:.2f}', transform=ax.transAxes, fontsize=9,
                    ha='right', bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8))

        ax.set_title(country, fontsize=10, fontweight='bold')
        ax.set_xlabel('Distance to Distress', fontsize=8)
        ax.set_ylabel('CDS Spread (bps)', fontsize=8)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.2)

    for j in range(len(countries), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle('Distance-to-Distress vs CDS Spread', fontsize=14, fontweight='bold')
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def plot_changes(df,
                 countries=None, ncols=3, save_path=None):
    if countries is None:
        countries = sorted(df['country'].unique())

    nrows = int(np.ceil(len(countries) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.atleast_2d(axes).flatten()

    for i, country in enumerate(countries):
        ax = axes[i]
        d = df[df['country'] == country].sort_values('date').dropna(subset=['distance_to_distress', 'cds_spread'])

        if len(d) < 20:
            ax.set_visible(False)
            continue

        dx = d['distance_to_distress'].diff()
        dy = d['cds_spread'].diff()
        mask = np.isfinite(dx) & np.isfinite(dy)
        dx, dy = dx[mask].values, dy[mask].values

        ax.scatter(dx, dy, s=8, alpha=0.4, color='#2C5F8A', edgecolors='none')

        # OLS fit
        slope, intercept = np.polyfit(dx, dy, 1)
        x_line = np.array([dx.min(), dx.max()])
        ax.plot(x_line, slope * x_line + intercept, 'k-', lw=1.5)

        # R²
        ss_res = np.sum((dy - (slope * dx + intercept)) ** 2)
        ss_tot = np.sum((dy - dy.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

        # Directional accuracy: d2 up → CDS should go down (and vice versa)
        nonzero = (dx != 0) & (dy != 0)
        if nonzero.sum() > 0:
            dir_acc = np.mean((dx[nonzero] > 0) == (dy[nonzero] < 0))
        else:
            dir_acc = np.nan

        ax.text(0.95, 0.95,
                f'R² = {r2:.2f}\nβ = {slope:.1f}\nDir. acc = {dir_acc:.0%}',
                transform=ax.transAxes, fontsize=8, ha='right', va='top',
                bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8))

        ax.axhline(0, color='gray', lw=0.5, ls='--')
        ax.axvline(0, color='gray', lw=0.5, ls='--')
        ax.set_title(country, fontsize=10, fontweight='bold')
        ax.set_xlabel('Δ d2', fontsize=8)
        ax.set_ylabel('Δ CDS (bps)', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.15)

    for j in range(len(countries), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle('Change in Distance-to-Distress vs Change in CDS Spread',
                 fontsize=13, fontweight='bold')
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def correlation_table(df, dd_col='distance_to_distress', cds_col='cds_spread',
                      horizons=[1, 4, 13], horizon_labels=None,
                      group_map=None, save_path=None):
    """
    Duyvesteyn & Martens (2015) Table 4 style correlation table.
    
    Correlations between Δd2 and ΔCDS at multiple horizons,
    per country, with significance stars and directional accuracy.
    
    Parameters
    ----------
    df : DataFrame with 'country', 'date', cds_col, and dd_col
    dd_col : d2 column name
    cds_col : CDS spread column name
    horizons : list of ints, periods for differencing (1=1-week, 4≈1-month, 13≈3-month)
    horizon_labels : list of str labels for horizons
    group_map : dict mapping country -> group label (e.g. {'Saudi Arabia': 'Oil Exporter'})
    save_path : if provided, saves table as CSV
    
    Returns
    -------
    display_df : formatted DataFrame
    raw_df : DataFrame with raw numeric values for further analysis
    """
    from scipy import stats

    if horizon_labels is None:
        horizon_labels = [f'{h}w' for h in horizons]

    countries = sorted(df['country'].unique())

    rows = []
    for country in countries:
        d = df[df['country'] == country].sort_values('date').copy()
        d['date'] = pd.to_datetime(d['date'])

        group = group_map.get(country, 'Control') if group_map else ''

        row = {'Country': country, 'Group': group, 'N': d[dd_col].notna().sum()}

        for h, h_label in zip(horizons, horizon_labels):
            d_cds = d[cds_col].diff(h)
            d_dd = d[dd_col].diff(h)

            valid = d_cds.notna() & d_dd.notna()
            n = valid.sum()

            if n < 20:
                row[f'ρ ({h_label})'] = np.nan
                row[f'p ({h_label})'] = np.nan
                row[f'Dir ({h_label})'] = np.nan
                continue

            x = d_dd[valid].values
            y = d_cds[valid].values

            rho, pval = stats.pearsonr(x, y)

            # Directional accuracy: d2 up → CDS down
            nonzero = (x != 0) & (y != 0)
            dir_acc = np.mean((x[nonzero] > 0) == (y[nonzero] < 0)) if nonzero.sum() > 0 else np.nan

            row[f'ρ ({h_label})'] = rho
            row[f'p ({h_label})'] = pval
            row[f'Dir ({h_label})'] = dir_acc

        rows.append(row)

    raw = pd.DataFrame(rows)

    # ── Group averages ──
    if group_map:
        for group_name in sorted(set(group_map.values())):
            group_rows = raw[raw['Group'] == group_name]
            if len(group_rows) == 0:
                continue
            avg_row = {'Country': f'── {group_name} avg ──', 'Group': group_name, 'N': ''}
            for h_label in horizon_labels:
                rho_vals = group_rows[f'ρ ({h_label})'].dropna()
                dir_vals = group_rows[f'Dir ({h_label})'].dropna()
                avg_row[f'ρ ({h_label})'] = rho_vals.mean() if len(rho_vals) > 0 else np.nan
                avg_row[f'p ({h_label})'] = np.nan
                avg_row[f'Dir ({h_label})'] = dir_vals.mean() if len(dir_vals) > 0 else np.nan
            rows.append(avg_row)

        # Overall average
        avg_row = {'Country': '── All avg ──', 'Group': '', 'N': ''}
        for h_label in horizon_labels:
            rho_vals = raw[f'ρ ({h_label})'].dropna()
            dir_vals = raw[f'Dir ({h_label})'].dropna()
            avg_row[f'ρ ({h_label})'] = rho_vals.mean() if len(rho_vals) > 0 else np.nan
            avg_row[f'p ({h_label})'] = np.nan
            avg_row[f'Dir ({h_label})'] = dir_vals.mean() if len(dir_vals) > 0 else np.nan
        rows.append(avg_row)

    full = pd.DataFrame(rows)

    # ── Format display version ──
    def _fmt_rho(rho, pval):
        if np.isnan(rho):
            return ''
        stars = ''
        if pval < 0.01:
            stars = '***'
        elif pval < 0.05:
            stars = '**'
        elif pval < 0.10:
            stars = '*'
        return f'{rho:.2f}{stars}'

    def _fmt_dir(d):
        if np.isnan(d):
            return ''
        return f'{d:.0%}'

    display = full[['Country', 'Group', 'N']].copy()
    for h_label in horizon_labels:
        display[f'ρ ({h_label})'] = [
            _fmt_rho(r, p) for r, p in zip(full[f'ρ ({h_label})'], full[f'p ({h_label})'])
        ]
        display[f'Dir ({h_label})'] = full[f'Dir ({h_label})'].apply(_fmt_dir)

    if not group_map:
        display = display.drop(columns=['Group'])

    if save_path:
        display.to_csv(save_path, index=False)

    return display, raw