import os
import numpy as np
import pandas as pd
from collections import Counter
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve
from lightgbm import LGBMClassifier
import matplotlib.pyplot as plt
import matplotlib as mpl


# =========================
# 0) Helpers
# =========================
def read_table_auto(path: str) -> pd.DataFrame:
    """Read csv/tsv with auto separator (',' or '\\t')."""
    df = pd.read_csv(path, sep=None, engine="python")
    return df

def normalize_dict(d: dict) -> dict:
    """Normalize values to sum=1 (safe)."""
    s = float(sum(d.values()))
    return d if s == 0 else {k: v / s for k, v in d.items()}

def infer_feature_cols(df: pd.DataFrame, y_col: str, drop_cols=("eid",)) -> list:
    """Infer feature columns: numeric cols excluding y and id-like cols."""
    drop_set = set([y_col] + list(drop_cols))
    cols = [c for c in df.columns if c not in drop_set]
    # keep numeric only (robust across datasets)
    num_cols = df[cols].select_dtypes(include=[np.number]).columns.tolist()
    return num_cols

def fit_lgbm(params: dict) -> LGBMClassifier:
    """Build a LightGBM binary classifier with sane defaults."""
    base = dict(
        objective="binary",
        metric="auc",
        verbosity=-1,
        n_jobs=4,
        random_state=2022,
    )
    base.update(params)
    return LGBMClassifier(**base)


# =========================
# 1) CV feature importance
# =========================
def cv_feature_importance(
    X: pd.DataFrame,
    y: pd.Series,
    params: dict,
    n_splits: int = 5,
    seed: int = 2022,
):
    """Compute CV-summed normalized importance (gain + split)."""
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    gain_sum, split_sum = Counter(), Counter()

    for tr_idx, _ in skf.split(X, y):
        model = fit_lgbm(params)
        model.fit(X.iloc[tr_idx], y.iloc[tr_idx])

        booster = model.booster_
        names = booster.feature_name()
        gain = booster.feature_importance(importance_type="gain")
        split = booster.feature_importance(importance_type="split")

        gain_sum += Counter(normalize_dict(dict(zip(names, gain))))
        split_sum += Counter(normalize_dict(dict(zip(names, split))))

    imp = pd.DataFrame({"Feature": list(gain_sum.keys())})
    imp["Gain_cv"] = imp["Feature"].map(gain_sum).fillna(0.0)
    imp["Split_cv"] = imp["Feature"].map(split_sum).fillna(0.0)
    imp = imp.sort_values("Gain_cv", ascending=False).reset_index(drop=True)
    return imp


# =========================
# 2) Stepwise AUC
# =========================
def stepwise_auc(
    X: pd.DataFrame,
    y: pd.Series,
    feature_order: list,
    params: dict,
    top_k: int = 50,
    n_splits: int = 5,
    seed: int = 2022,
):
    """
    Add features one-by-one (in given order) and compute CV AUC.
    Returns a table: Feature, AUC_mean, AUC_std, AUC0..AUC(n-1)
    """
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    chosen, rows = [], []

    for f in feature_order[:top_k]:
        chosen.append(f)
        Xk = X[chosen]

        aucs = []
        for tr_idx, te_idx in skf.split(Xk, y):
            model = fit_lgbm(params)
            model.fit(Xk.iloc[tr_idx], y.iloc[tr_idx])
            prob = model.predict_proba(Xk.iloc[te_idx])[:, 1]
            aucs.append(roc_auc_score(y.iloc[te_idx], prob))

        rows.append([f, np.mean(aucs), np.std(aucs)] + aucs)

    cols = ["Feature", "AUC_mean", "AUC_std"] + [f"AUC{i}" for i in range(n_splits)]
    out = pd.DataFrame(rows, columns=cols)
    return out


# =========================
# 3) Mean ROC with band
# =========================
def cv_mean_roc(
    X: pd.DataFrame,
    y: pd.Series,
    params: dict,
    n_splits: int = 5,
    seed: int = 2022,
    n_grid: int = 200,
    out_pdf: str = None,
):
    """Plot mean ROC and ±2*std band using CV."""
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    mean_fpr = np.linspace(0, 1, n_grid)
    tprs, aucs = [], []

    for tr_idx, te_idx in skf.split(X, y):
        model = fit_lgbm(params)
        model.fit(X.iloc[tr_idx], y.iloc[tr_idx])
        prob = model.predict_proba(X.iloc[te_idx])[:, 1]

        fpr, tpr, _ = roc_curve(y.iloc[te_idx], prob)
        tpr_i = np.interp(mean_fpr, fpr, tpr)
        tpr_i[0] = 0.0
        tprs.append(tpr_i)
        aucs.append(roc_auc_score(y.iloc[te_idx], prob))

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0

    tpr_std = np.std(tprs, axis=0)
    upper = np.minimum(mean_tpr + 2 * tpr_std, 1)
    lower = np.maximum(mean_tpr - 2 * tpr_std, 0)

    plt.figure(figsize=(6, 6))
    plt.plot(mean_fpr, mean_tpr, linewidth=2,
             label=f"Mean ROC (AUC={np.mean(aucs):.2f} ± {np.std(aucs):.2f})")
    plt.plot([0, 1], [0, 1], "k--", linewidth=1)
    plt.fill_between(mean_fpr, lower, upper, alpha=0.2)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.legend(loc="lower right")
    plt.tight_layout()

    if out_pdf:
        plt.savefig(out_pdf)
    plt.show()

    return dict(mean_auc=float(np.mean(aucs)), std_auc=float(np.std(aucs)))


# =========================
# 4) Plot AUC-step + importance aligned
# =========================
def plot_auc_with_importance(
    imp_df: pd.DataFrame,
    auc_df: pd.DataFrame,
    imp_col: str = "Gain_cv",   # or "Split_cv"
    nb_core: int = 33,
    out_pdf: str = None,
    out_png: str = None,
):
    """Bar = importance, line = stepwise AUC (aligned by feature order)."""
    # Align importance to auc order
    imp_map = imp_df.set_index("Feature")[imp_col]
    plot_df = auc_df.copy()
    plot_df["Imp"] = plot_df["Feature"].map(imp_map)

    # CI using std (simple; matches your logic)
    plot_df["AUC_lower"] = plot_df["AUC_mean"] - 1.96 * plot_df["AUC_std"]
    plot_df["AUC_upper"] = plot_df["AUC_mean"] + 1.96 * plot_df["AUC_std"]

    x = np.arange(len(plot_df))
    fig, ax1 = plt.subplots(figsize=(18, 7))

    # Gradient blue bars
    cmap = mpl.colormaps["Blues"]
    colors = [cmap(0.25 + 0.65 * (1 - i / max(1, len(plot_df) - 1))) for i in range(len(plot_df))]
    ax1.bar(x, plot_df["Imp"], color=colors, edgecolor="none", alpha=0.9)
    ax1.set_ylabel("Predictor Importance", fontsize=16, fontweight="bold")
    ax1.grid(axis="y", linestyle="--", alpha=0.4)
    ax1.set_axisbelow(True)

    ax1.set_xticks(x)
    ax1.set_xticklabels(plot_df["Feature"], rotation=30, ha="right", fontsize=10)
    for i, tick in enumerate(ax1.get_xticklabels()):
        tick.set_color("red" if i < nb_core else "black")

    # AUC curve (right axis)
    ax2 = ax1.twinx()
    ax2.plot(x, plot_df["AUC_mean"], marker="o", linewidth=2.5)
    ax2.fill_between(x, plot_df["AUC_lower"], plot_df["AUC_upper"], alpha=0.2)
    ax2.set_ylabel("Cumulative AUC", fontsize=16, fontweight="bold")
    ax2.set_ylim(
        min(plot_df["AUC_lower"].min(), plot_df["AUC_mean"].min()) - 0.02,
        max(plot_df["AUC_upper"].max(), plot_df["AUC_mean"].max()) + 0.02
    )

    plt.tight_layout()
    if out_pdf:
        plt.savefig(out_pdf)
    if out_png:
        plt.savefig(out_png, dpi=300)
    plt.show()

    return plot_df


# =========================================================
# MAIN (minimal edits needed)
# =========================================================
if __name__ == "__main__":

    os.chdir("/home/ug1268u4/5.UKB/Project/1.蛋白组_migraine/1.Prediction")

    in_path = "resid_ukb_protein_pre2.tsv"
    y_col = "migraine_status"
    id_cols = ("eid",)  # add more if needed

    df = read_table_auto(in_path)

    # Build X/y robustly
    y = df[y_col].astype(int)
    feat_cols = infer_feature_cols(df, y_col=y_col, drop_cols=id_cols)
    X = df[feat_cols]

    os.makedirs("result_migraine", exist_ok=True)

    # ---- A) CV importance ----
    imp_params = dict(n_estimators=500, max_depth=15, num_leaves=10,
                      subsample=0.7, learning_rate=0.01, colsample_bytree=0.7)
    imp = cv_feature_importance(X, y, params=imp_params)
    imp.to_csv("result_migraine/importance_cv.csv", index=False)

    # Choose ordering for stepwise
    order_by = "Gain_cv"  # "Split_cv" also ok
    feat_order = imp.sort_values(order_by, ascending=False)["Feature"].tolist()

    # ---- B) Stepwise AUC (topK) ----
    step_params = dict(n_estimators=800, max_depth=15, num_leaves=25,
                       subsample=0.7, learning_rate=0.01, colsample_bytree=0.7)
    topK = 50
    auc_step = stepwise_auc(X, y, feature_order=feat_order, params=step_params, top_k=topK)
    auc_step.to_csv("result_migraine/auc_stepwise.csv", index=False)

    # ---- C) Final ROC using first nb_f features ----
    nb_f = 33  # choose by elbow/turning point
    final_features = feat_order[:nb_f]
    final_params = dict(n_estimators=1000, max_depth=5, num_leaves=30,
                        subsample=0.7, learning_rate=0.01, colsample_bytree=0.7,
                        is_unbalance=True)

    cv_mean_roc(X[final_features], y, params=final_params,
                out_pdf="result_migraine/ROC_mean_cv.pdf")

    # ---- D) Plot AUC-step + importance ----
    plot_df = plot_auc_with_importance(
        imp_df=imp,
        auc_df=auc_step,
        imp_col=order_by,
        nb_core=nb_f,
        out_pdf="result_migraine/Fig_importance_cumAUC.pdf",
        out_png="result_migraine/Fig_importance_cumAUC.png"
    )
    plot_df.to_csv("result_migraine/plot_df_auc_imp_aligned.csv", index=False)