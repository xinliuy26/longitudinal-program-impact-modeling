# ==========================================
# 8. 森林图 Forest plot + 功效曲线 Power curves
# Paste as a new section at the end of longitudinal_impact_modeling.ipynb
# (after section 7, PSM). Requires df_drp, df, att_pairs from sections 1-7.
# Recomputes every number from the fitted models - nothing is hard-coded.
# Saves forest_plot.png and power_curve.png next to the notebook.
# ==========================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# --- 8.1 Forest plot: SPIES coefficient across M1-M5 (DRP) + PSM ATT ---
formulas = {
    'M1\n(pretest only)': 'drp_post_z ~ spies + drp_pre_z',
    'M2\n(+ ethnicity)': 'drp_post_z ~ spies + drp_pre_z + white',
    'M3\n(+ gender)': 'drp_post_z ~ spies + drp_pre_z + white + female',
    'M4\n(+ IEP)': 'drp_post_z ~ spies + drp_pre_z + white + female + iep',
    'M5\n(full model)': 'drp_post_z ~ spies + drp_pre_z + white + female + iep + ell',
}
rows = []
for label, f in formulas.items():
    m = smf.ols(f, data=df_drp).fit()
    ci = m.conf_int().loc['spies']
    rows.append((label, m.params['spies'], ci[0], ci[1], m.pvalues['spies']))
# PSM ATT with 95% CI from the paired differences (section 7)
att_se = att_pairs.std(ddof=1) / np.sqrt(len(att_pairs))
att_t = stats.ttest_1samp(att_pairs, 0)
rows.append(('PSM ATT\n(55 matched pairs)', att, att - 1.96 * att_se,
             att + 1.96 * att_se, att_t.pvalue))
for r in rows:
    print(f"{r[0].replace(chr(10), ' '):<28} coef={r[1]:.3f}  95% CI=[{r[2]:.3f}, {r[3]:.3f}]  p={r[4]:.3f}")

fig, ax = plt.subplots(figsize=(7.2, 4.2))
y = np.arange(len(rows))[::-1]
coef = np.array([r[1] for r in rows])
lo = np.array([r[2] for r in rows])
hi = np.array([r[3] for r in rows])
ax.errorbar(coef, y, xerr=[coef - lo, hi - coef], fmt='o', color='#2471a3',
            ecolor='#2471a3', elinewidth=1.6, capsize=4, markersize=6)
ax.scatter([coef[-1]], [y[-1]], color='#1e8449', s=60, zorder=3)
ax.axvline(0, color='grey', lw=1.0)
ax.axvline(0.25, color='#c0392b', lw=1.0, ls='--')
ax.text(0.25, -0.75, 'paper estimate 0.25 SD', color='#c0392b', fontsize=8, ha='center')
ax.set_yticks(y)
ax.set_yticklabels([r[0] for r in rows], fontsize=9)
for yi, r in zip(y, rows):
    ax.text(r[1] + 0.03, yi + 0.28,
            f'{r[1]:.2f}  [{r[2]:.2f}, {r[3]:.2f}]  p = {r[4]:.3f}', fontsize=8, va='center')
ax.set_xlabel('SPIES participation effect on DRP post-test z-score (SD)')
ax.set_title('Forest plot — SPIES effect across progressive OLS models (M1–M5) and PSM\n'
             'Simulated data, DRP sample N = 131 (23 treated / 108 control); 95% CIs', fontsize=10)
ax.set_xlim(-0.35, 1.05)
ax.set_ylim(-1.1, 5.6)
fig.tight_layout()
fig.savefig('forest_plot.png', dpi=150)
print('Saved: forest_plot.png')

# --- 8.2 Power curves: DRP vs PSSA, same MDE formula as section 6 ---
def power(delta, n_treat, n_control, r2, alpha=0.05):
    se = np.sqrt((1 - r2) * (1 / n_treat + 1 / n_control))
    z = stats.norm.ppf(1 - alpha / 2)
    return stats.norm.cdf(delta / se - z) + stats.norm.cdf(-delta / se - z)

d = np.linspace(0, 1.2, 241)
p_drp = power(d, 23, 108, 0.273)
p_pssa = power(d, 11, 41, 0.308)
print(f"Power at the paper effect (0.25 SD): DRP = {power(0.25, 23, 108, 0.273):.3f}, "
      f"PSSA at 0.06 SD = {power(0.06, 11, 41, 0.308):.3f}")

fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.plot(d, p_drp, color='#2471a3', lw=2, label='DRP  (N = 131: 23 / 108, R² = 0.273)')
ax.plot(d, p_pssa, color='#c0392b', lw=2, label='PSSA  (N = 52: 11 / 41, R² = 0.308)')
ax.axhline(0.80, color='grey', lw=1.0, ls='--')
ax.axvline(0.25, color='grey', lw=1.0, ls=':')
ax.plot([0.549], [0.80], 'o', color='#2471a3')
ax.plot([0.791], [0.80], 'o', color='#c0392b')
ax.annotate('DRP MDE 0.55 SD', (0.549, 0.80), textcoords='offset points', xytext=(8, 6),
            fontsize=9, color='#2471a3')
ax.annotate('PSSA MDE 0.79 SD', (0.791, 0.80), textcoords='offset points', xytext=(8, -16),
            fontsize=9, color='#c0392b')
ax.annotate('paper effect 0.25 SD\n→ power ≈ 25% (DRP)', (0.25, 0.248),
            textcoords='offset points', xytext=(10, 8), fontsize=9)
ax.plot([0.25], [0.248], 'o', color='grey')
ax.set_xlabel('True effect size (SD)')
ax.set_ylabel('Statistical power (two-sided α = 0.05)')
ax.set_title('Power curves — DRP vs PSSA samples\nNormal-approximation MDE, same formula as notebook §6', fontsize=10)
ax.set_xlim(0, 1.2)
ax.set_ylim(0, 1.02)
ax.legend(loc='lower right', fontsize=9)
fig.tight_layout()
fig.savefig('power_curve.png', dpi=150)
print('Saved: power_curve.png')
