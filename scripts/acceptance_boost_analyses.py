#!/usr/bin/env python3
"""
Acceptance-boost analyses using existing dataset only (no external validation, no wet-lab).

Outputs:
  - results/acceptance_covariate_associations.csv
  - results/acceptance_multivariable_cox.csv
  - results/acceptance_ml_metrics.csv
  - results/acceptance_ml_small_panel_metrics.csv
  - results/acceptance_network_subtype_specific_summary.csv
  - results/acceptance_network_permutation_summary.csv

Also prints a concise textual summary to stdout.
"""

from __future__ import annotations

import os
import math
import numpy as np
import pandas as pd

from scipy import stats
from statsmodels.stats.multitest import multipletests

from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    balanced_accuracy_score,
    matthews_corrcoef,
    brier_score_loss,
)
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from lifelines import CoxPHFitter


ROOT = "/Users/senol/Desktop/GBM"
PROJ = f"{ROOT}/Lipidomics_Metabolomics_Integration"
RESULTS_DIR = f"{PROJ}/results"
os.makedirs(RESULTS_DIR, exist_ok=True)


def read_quant(name: str) -> pd.DataFrame:
    return pd.read_csv(f"{ROOT}/{name}", sep="\t")


def extract_sample_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if isinstance(c, str) and c.startswith("C3")]


def prepare_matrix(df: pd.DataFrame, sample_cols: list[str]) -> pd.DataFrame:
    # index by feature id/name (first column)
    idx_col = "ENTITY_STABLE_ID" if "ENTITY_STABLE_ID" in df.columns else df.columns[0]
    df = df.set_index(idx_col)
    m = df[sample_cols].T.apply(pd.to_numeric, errors="coerce")
    m = m.loc[:, m.isna().mean() < 0.5]
    m = m.fillna(m.median())
    return m


def rank_df(df: pd.DataFrame) -> pd.DataFrame:
    # Spearman via rank transform then Pearson
    return df.rank(axis=0, method="average")


def corr_ranked(x_rank: np.ndarray, y_rank: np.ndarray) -> np.ndarray:
    """
    Pearson correlation between columns of x_rank (n x p) and y_rank (n x q),
    returning p x q correlation matrix.
    """
    # center
    x = x_rank - x_rank.mean(axis=0, keepdims=True)
    y = y_rank - y_rank.mean(axis=0, keepdims=True)
    # normalize
    x_norm = np.sqrt((x * x).sum(axis=0, keepdims=True))
    y_norm = np.sqrt((y * y).sum(axis=0, keepdims=True))
    x = x / np.where(x_norm == 0, 1.0, x_norm)
    y = y / np.where(y_norm == 0, 1.0, y_norm)
    return x.T @ y


def normal_approx_pvals(r: np.ndarray, n: int) -> np.ndarray:
    """
    Fast p-value approximation for correlation coefficients using t distribution.
    """
    df = max(n - 2, 1)
    r = np.clip(r, -0.999999, 0.999999)
    t = r * np.sqrt(df / (1.0 - r * r))
    p = 2.0 * stats.t.sf(np.abs(t), df=df)
    return p


def cohen_d(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    nx = np.sum(np.isfinite(x))
    ny = np.sum(np.isfinite(y))
    if nx < 2 or ny < 2:
        return np.nan
    vx = np.nanvar(x, ddof=1)
    vy = np.nanvar(y, ddof=1)
    s = np.sqrt(((nx - 1) * vx + (ny - 1) * vy) / max(nx + ny - 2, 1))
    if s == 0:
        return 0.0
    return (np.nanmean(x) - np.nanmean(y)) / s


def main() -> None:
    # ---------------------------
    # Load matrices + clinical
    # ---------------------------
    metabolome = read_quant("metabolome_quantification.txt")
    lipid_pos = read_quant("lipidome_positive_quantification.txt")
    lipid_neg = read_quant("lipidome_negative_quantification.txt")

    clinical = pd.read_csv(f"{ROOT}/clinical_patient.txt", sep="\t", skiprows=4)
    clinical.columns = clinical.columns.str.strip()
    clinical = clinical[clinical["PATIENT_ID"].notna()].copy()

    # subtype assignments (use NMF as primary)
    subtype_df = pd.read_csv(f"{PROJ}/results/cluster_assignments.csv")
    subtype_df = subtype_df.set_index("sample_id")
    subtype = subtype_df["nmf_subtype"].astype(int)

    # common samples across omics + clinical + subtype
    m_s = set(extract_sample_cols(metabolome))
    lp_s = set(extract_sample_cols(lipid_pos))
    ln_s = set(extract_sample_cols(lipid_neg))
    c_s = set(clinical["PATIENT_ID"].astype(str).tolist())
    s_s = set(subtype.index.astype(str).tolist())
    common = sorted(list(m_s & lp_s & ln_s & c_s & s_s))

    m_mat = prepare_matrix(metabolome, common)
    lp_mat = prepare_matrix(lipid_pos, common)
    ln_mat = prepare_matrix(lipid_neg, common)
    l_mat = pd.concat([lp_mat, ln_mat], axis=1)

    # build combined feature matrix X and label y
    X = pd.concat([m_mat, l_mat], axis=1).loc[common]
    y = subtype.loc[common].values

    # clinical covariates aligned
    c = clinical.set_index("PATIENT_ID").loc[common].copy()
    c["subtype"] = y

    # harmonize numeric columns
    for col in ["AGE", "BMI", "HEIGHT", "WEIGHT", "PATH_DIAG_TO_LAST_CONTACT_DAYS", "PATH_DIAG_TO_DEATH_DAYS", "PATH_DIAG_TO_TUMOR_REOCCUR_DAYS"]:
        if col in c.columns:
            c[col] = pd.to_numeric(c[col], errors="coerce")

    # encode sex
    c["SEX_BIN"] = c["SEX"].astype(str).str.lower().map({"male": 1, "female": 0})

    # recurrence
    def _to_bool(v: object) -> float:
        if pd.isna(v):
            return np.nan
        s = str(v).strip().upper()
        if s == "TRUE":
            return 1.0
        if s == "FALSE":
            return 0.0
        return np.nan

    c["RECURRENCE_BIN"] = c["TUMOR_REOCCUR_AFTER_TREATMENT"].map(_to_bool) if "TUMOR_REOCCUR_AFTER_TREATMENT" in c.columns else np.nan

    # survival time/event
    time = np.where(c["PATH_DIAG_TO_DEATH_DAYS"].notna(), c["PATH_DIAG_TO_DEATH_DAYS"], c["PATH_DIAG_TO_LAST_CONTACT_DAYS"])
    event = np.where(c["PATH_DIAG_TO_DEATH_DAYS"].notna(), 1, 0)
    surv = pd.DataFrame({"time": time, "event": event}, index=c.index)
    surv = surv.dropna(subset=["time"])

    # ---------------------------
    # (A) Covariate associations
    # ---------------------------
    rows = []
    s0 = c[c["subtype"] == 0]
    s1 = c[c["subtype"] == 1]

    # numeric: AGE, BMI
    for col in ["AGE", "BMI"]:
        if col not in c.columns:
            continue
        x0 = s0[col].dropna().values
        x1 = s1[col].dropna().values
        # Welch t-test + MWU (report both; choose t-test p as primary)
        t_p = stats.ttest_ind(x0, x1, equal_var=False, nan_policy="omit").pvalue if len(x0) and len(x1) else np.nan
        try:
            mwu_p = stats.mannwhitneyu(x0, x1, alternative="two-sided").pvalue if len(x0) and len(x1) else np.nan
        except Exception:
            mwu_p = np.nan
        rows.append(
            {
                "variable": col,
                "type": "numeric",
                "n_subtype0": int(np.sum(np.isfinite(x0))),
                "n_subtype1": int(np.sum(np.isfinite(x1))),
                "mean_subtype0": float(np.nanmean(x0)) if len(x0) else np.nan,
                "mean_subtype1": float(np.nanmean(x1)) if len(x1) else np.nan,
                "cohens_d": cohen_d(x0, x1),
                "p_welch_t": t_p,
                "p_mannwhitney": mwu_p,
            }
        )

    # categorical: SEX, RECURRENCE
    def cat_test(col: str, a: pd.Series, b: pd.Series) -> dict:
        # 2x2 Fisher exact
        a = a.dropna().astype(int)
        b = b.dropna().astype(int)
        if a.empty or b.empty:
            return {"p_fisher": np.nan, "or_fisher": np.nan, "n0": int(a.shape[0]), "n1": int(b.shape[0])}
        n00 = int((a == 0).sum())
        n01 = int((a == 1).sum())
        n10 = int((b == 0).sum())
        n11 = int((b == 1).sum())
        table = np.array([[n01, n00], [n11, n10]])  # [yes,no] by subtype; orient for OR
        try:
            orv, pv = stats.fisher_exact(table)
        except Exception:
            orv, pv = (np.nan, np.nan)
        return {
            "p_fisher": pv,
            "or_fisher": orv,
            "n_subtype0": int(a.shape[0]),
            "n_subtype1": int(b.shape[0]),
            "prop_subtype0": float((a == 1).mean()) if a.shape[0] else np.nan,
            "prop_subtype1": float((b == 1).mean()) if b.shape[0] else np.nan,
        }

    for col, name in [("SEX_BIN", "SEX"), ("RECURRENCE_BIN", "RECURRENCE")]:
        if col not in c.columns:
            continue
        out = cat_test(col, s0[col], s1[col])
        rows.append(
            {
                "variable": name,
                "type": "binary",
                "n_subtype0": out.get("n_subtype0"),
                "n_subtype1": out.get("n_subtype1"),
                "mean_subtype0": out.get("prop_subtype0"),
                "mean_subtype1": out.get("prop_subtype1"),
                "cohens_d": np.nan,
                "p_welch_t": np.nan,
                "p_mannwhitney": np.nan,
                "p_fisher": out.get("p_fisher"),
                "or_fisher": out.get("or_fisher"),
            }
        )

    cov_df = pd.DataFrame(rows)
    # multiple testing correction (use primary p per row)
    primary_p = cov_df["p_welch_t"].copy()
    primary_p = primary_p.where(primary_p.notna(), cov_df.get("p_fisher"))
    cov_df["p_primary"] = primary_p
    cov_df["fdr_primary"] = np.nan
    mask = cov_df["p_primary"].notna()
    if mask.any():
        cov_df.loc[mask, "fdr_primary"] = multipletests(cov_df.loc[mask, "p_primary"], method="fdr_bh")[1]
    cov_df.to_csv(f"{RESULTS_DIR}/acceptance_covariate_associations.csv", index=False)

    # ---------------------------
    # (B) Multivariable Cox: subtype + age + sex + BMI
    # ---------------------------
    cox_df = pd.DataFrame(index=surv.index)
    cox_df["time"] = surv.loc[cox_df.index, "time"].astype(float)
    cox_df["event"] = surv.loc[cox_df.index, "event"].astype(int)
    cox_df["subtype"] = c.loc[cox_df.index, "subtype"].astype(int)
    # add covariates if sufficiently complete
    for col, out_col in [("AGE", "age"), ("SEX_BIN", "sex_male"), ("BMI", "bmi")]:
        if col in c.columns:
            cox_df[out_col] = c.loc[cox_df.index, col].astype(float)
    # drop rows with missing covariates
    cox_df_clean = cox_df.dropna()
    cox_out = None
    if cox_df_clean.shape[0] >= 30 and cox_df_clean["event"].sum() >= 10:
        cph = CoxPHFitter()
        cph.fit(cox_df_clean, duration_col="time", event_col="event")
        cox_out = cph.summary.reset_index().rename(columns={"index": "covariate"})
        cox_out.to_csv(f"{RESULTS_DIR}/acceptance_multivariable_cox.csv", index=False)
    else:
        pd.DataFrame({"note": ["Insufficient complete cases for multivariable Cox"]}).to_csv(
            f"{RESULTS_DIR}/acceptance_multivariable_cox.csv", index=False
        )

    # ---------------------------
    # (C) ML metrics: repeated stratified CV
    # ---------------------------
    rkf = RepeatedStratifiedKFold(n_splits=5, n_repeats=50, random_state=42)
    model = Pipeline(
        steps=[
            ("scaler", StandardScaler(with_mean=True, with_std=True)),
            ("clf", LogisticRegression(max_iter=5000, class_weight="balanced", solver="liblinear")),
        ]
    )

    ml_rows = []
    y_true_all = []
    y_prob_all = []
    y_pred_all = []

    Xv = X.values.astype(float)
    for split_i, (tr, te) in enumerate(rkf.split(Xv, y)):
        Xtr, Xte = Xv[tr], Xv[te]
        ytr, yte = y[tr], y[te]
        model.fit(Xtr, ytr)
        prob = model.predict_proba(Xte)[:, 1]
        pred = (prob >= 0.5).astype(int)
        # metrics
        ml_rows.append(
            {
                "split": split_i,
                "roc_auc": roc_auc_score(yte, prob) if len(np.unique(yte)) == 2 else np.nan,
                "pr_auc": average_precision_score(yte, prob) if len(np.unique(yte)) == 2 else np.nan,
                "balanced_accuracy": balanced_accuracy_score(yte, pred),
                "mcc": matthews_corrcoef(yte, pred),
                "brier": brier_score_loss(yte, prob),
            }
        )
        y_true_all.append(yte)
        y_prob_all.append(prob)
        y_pred_all.append(pred)

    ml_df = pd.DataFrame(ml_rows)
    ml_df.to_csv(f"{RESULTS_DIR}/acceptance_ml_metrics.csv", index=False)

    # ---------------------------
    # (D) Small panel (top-k by effect size within train fold)
    # ---------------------------
    k_panel = 10
    panel_rows = []
    for split_i, (tr, te) in enumerate(rkf.split(Xv, y)):
        Xtr, Xte = X.iloc[tr], X.iloc[te]
        ytr, yte = y[tr], y[te]

        # effect size by Cohen's d on training set
        d = []
        for j, col in enumerate(X.columns):
            x0 = Xtr.loc[ytr == 0, col].values
            x1 = Xtr.loc[ytr == 1, col].values
            d.append(abs(cohen_d(x0, x1)))
        d = np.asarray(d)
        top_idx = np.argsort(-np.nan_to_num(d, nan=-1.0))[:k_panel]
        cols = X.columns[top_idx].tolist()

        small_model = Pipeline(
            steps=[
                ("scaler", StandardScaler(with_mean=True, with_std=True)),
                ("clf", LogisticRegression(max_iter=5000, class_weight="balanced", solver="liblinear")),
            ]
        )
        small_model.fit(Xtr[cols].values, ytr)
        prob = small_model.predict_proba(Xte[cols].values)[:, 1]
        pred = (prob >= 0.5).astype(int)
        panel_rows.append(
            {
                "split": split_i,
                "k": k_panel,
                "roc_auc": roc_auc_score(yte, prob) if len(np.unique(yte)) == 2 else np.nan,
                "pr_auc": average_precision_score(yte, prob) if len(np.unique(yte)) == 2 else np.nan,
                "balanced_accuracy": balanced_accuracy_score(yte, pred),
                "mcc": matthews_corrcoef(yte, pred),
                "brier": brier_score_loss(yte, prob),
            }
        )

    panel_df = pd.DataFrame(panel_rows)
    panel_df.to_csv(f"{RESULTS_DIR}/acceptance_ml_small_panel_metrics.csv", index=False)

    # ---------------------------
    # (E) Network: subtype-specific summary + permutation test for edge count
    # ---------------------------
    # Use only metabolite vs lipid edges, Spearman threshold |r|>=0.3, FDR<0.05
    m_rank = rank_df(m_mat.loc[common])
    l_rank = rank_df(l_mat.loc[common])

    def edge_count_for_samples(sample_ids: list[str]) -> dict:
        mr = m_rank.loc[sample_ids].values
        lr = l_rank.loc[sample_ids].values
        n = mr.shape[0]
        R = corr_ranked(lr, mr)  # lipids x metabolites
        P = normal_approx_pvals(R, n=n)
        p_flat = P.reshape(-1)
        fdr = multipletests(p_flat, method="fdr_bh")[1].reshape(P.shape)
        sig = (np.abs(R) >= 0.3) & (fdr < 0.05)
        return {
            "n_samples": int(n),
            "n_lipids": int(lr.shape[1]),
            "n_metabolites": int(mr.shape[1]),
            "n_edges_sig": int(sig.sum()),
            "r_median_abs_sig": float(np.median(np.abs(R[sig]))) if sig.any() else np.nan,
            "r_max_abs": float(np.max(np.abs(R))),
        }

    s0_ids = [sid for sid in common if subtype.loc[sid] == 0]
    s1_ids = [sid for sid in common if subtype.loc[sid] == 1]

    net0 = edge_count_for_samples(s0_ids)
    net1 = edge_count_for_samples(s1_ids)

    net_summary = pd.DataFrame(
        [
            {"group": "Subtype0", **net0},
            {"group": "Subtype1", **net1},
            {"group": "All", **edge_count_for_samples(common)},
        ]
    )
    net_summary.to_csv(f"{RESULTS_DIR}/acceptance_network_subtype_specific_summary.csv", index=False)

    # permutation: shuffle metabolome sample labels relative to lipids (break cross-omics pairing)
    rng = np.random.default_rng(123)
    n_perm = 200
    obs_edges = int(net_summary.loc[net_summary["group"] == "All", "n_edges_sig"].values[0])
    perm_edges = []

    lr_all = l_rank.loc[common].values
    mr_all = m_rank.loc[common].values
    n_all = lr_all.shape[0]

    for _ in range(n_perm):
        perm = rng.permutation(n_all)
        mr_perm = mr_all[perm, :]
        R = corr_ranked(lr_all, mr_perm)
        P = normal_approx_pvals(R, n=n_all)
        fdr = multipletests(P.reshape(-1), method="fdr_bh")[1].reshape(P.shape)
        sig = (np.abs(R) >= 0.3) & (fdr < 0.05)
        perm_edges.append(int(sig.sum()))

    perm_edges = np.asarray(perm_edges, dtype=int)
    p_empirical = float((np.sum(perm_edges >= obs_edges) + 1) / (n_perm + 1))
    perm_summary = pd.DataFrame(
        [
            {
                "n_perm": n_perm,
                "observed_edges": obs_edges,
                "perm_edges_mean": float(np.mean(perm_edges)),
                "perm_edges_sd": float(np.std(perm_edges, ddof=1)),
                "perm_edges_max": int(np.max(perm_edges)),
                "empirical_p_perm_ge_obs": p_empirical,
            }
        ]
    )
    perm_summary.to_csv(f"{RESULTS_DIR}/acceptance_network_permutation_summary.csv", index=False)

    # ---------------------------
    # Print summary for the user
    # ---------------------------
    def _fmt(x: float) -> str:
        if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
            return "NA"
        if abs(x) >= 1000:
            return f"{x:.0f}"
        if abs(x) >= 10:
            return f"{x:.2f}"
        return f"{x:.3f}"

    print("\n" + "=" * 80)
    print("ACCEPTANCE-BOOST ANALYSES SUMMARY (dataset-only)")
    print("=" * 80)
    print(f"Samples with subtype + clinical + omics: n={len(common)} (Subtype0={len(s0_ids)}, Subtype1={len(s1_ids)})")
    print("\n[Covariates by subtype] (primary p with BH-FDR in results/acceptance_covariate_associations.csv)")
    for v in ["AGE", "BMI", "SEX", "RECURRENCE"]:
        r = cov_df[cov_df["variable"] == v]
        if r.empty:
            continue
        r = r.iloc[0].to_dict()
        print(
            f"  - {v}: mean0={_fmt(r.get('mean_subtype0'))}, mean1={_fmt(r.get('mean_subtype1'))}, "
            f"p={_fmt(r.get('p_primary'))}, FDR={_fmt(r.get('fdr_primary'))}"
        )

    if cox_out is not None and not cox_out.empty:
        st = cox_out[cox_out["covariate"] == "subtype"]
        if not st.empty:
            st = st.iloc[0].to_dict()
            print("\n[Multivariable Cox OS] subtype adjusted for age/sex/BMI (complete-case)")
            print(f"  - HR(subtype1 vs 0)={_fmt(st.get('exp(coef)'))}, p={_fmt(st.get('p'))}")
        print(f"  - Complete cases used: n={cox_df_clean.shape[0]}, events={int(cox_df_clean['event'].sum())}")
    else:
        print("\n[Multivariable Cox OS] Not run (insufficient complete cases/events).")

    def summarize_metrics(df: pd.DataFrame, label: str) -> None:
        print(f"\n[{label}] 5x50 repeated stratified CV (logistic regression, class_weight=balanced)")
        for m in ["roc_auc", "pr_auc", "balanced_accuracy", "mcc", "brier"]:
            print(f"  - {m}: mean={_fmt(float(df[m].mean()))}, sd={_fmt(float(df[m].std(ddof=1)))}")

    summarize_metrics(ml_df, "ML (all features)")
    summarize_metrics(panel_df, f"ML (small panel k={k_panel})")

    print("\n[Network edges |r|>=0.3 & FDR<0.05] (metabolite-lipid)")
    for _, row in net_summary.iterrows():
        print(
            f"  - {row['group']}: n={int(row['n_samples'])}, edges={int(row['n_edges_sig'])}, "
            f"median|r|_sig={_fmt(row['r_median_abs_sig'])}"
        )
    print("\n[Network permutation] shuffle cross-omics pairing (n=200)")
    print(
        f"  - observed_edges={obs_edges}, perm_mean={_fmt(perm_summary.loc[0,'perm_edges_mean'])}, "
        f"perm_max={int(perm_summary.loc[0,'perm_edges_max'])}, empirical_p={_fmt(p_empirical)}"
    )
    print("\nSaved outputs under results/: acceptance_*.csv")


if __name__ == "__main__":
    main()

