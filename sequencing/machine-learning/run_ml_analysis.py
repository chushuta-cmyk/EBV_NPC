"""Run the EBV immune-cell machine-learning analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.base import clone
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.naive_bayes import GaussianNB
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


HIGH_CP = {
    "GSM7024387", "GSM7024389", "GSM7024391", "GSM7024392", "GSM7024394",
    "GSM7024396", "GSM7024398", "GSM7024404", "GSM7024407", "GSM7024413",
    "GSM7024414", "GSM7024422", "GSM7024385", "GSM7024400", "GSM7024420",
    "GSM7024409", "GSM7024405",
}
LOW_CP = {"GSM7024418", "GSM7024402", "GSM7024424", "GSM7024411", "GSM7024416"}


def parse_args() -> argparse.Namespace:
    project_root = Path(__file__).resolve().parents[3]
    data_dir = project_root / "EBV_recoding"
    parser = argparse.ArgumentParser()
    parser.add_argument("--cibersort", type=Path, default=data_dir / "CIBERSORT_Results_Filtered.csv")
    parser.add_argument("--genes", type=Path, default=data_dir / "Hub_Genes_Expression.csv")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "output")
    return parser.parse_args()


def load_data(path: Path) -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame]:
    data = pd.read_csv(path)
    first_column = data.columns[0]
    if not first_column or first_column.startswith("Unnamed"):
        data = data.rename(columns={first_column: "Sample_ID"})
    if "Sample_ID" not in data.columns:
        raise ValueError("The CIBERSORT file needs a Sample_ID column.")

    labels = data["Sample_ID"].map(
        lambda sample: 1 if sample in HIGH_CP else (0 if sample in LOW_CP else np.nan)
    )
    clean = data.loc[labels.notna()].copy()
    y = labels.loc[labels.notna()].astype(int).to_numpy()
    excluded = {"Sample_ID", "P-value", "Correlation", "RMSE", "Group"}
    feature_names = [column for column in clean.columns if column not in excluded]
    features = clean[feature_names].apply(pd.to_numeric, errors="coerce")
    if features.isna().any().any():
        bad = features.columns[features.isna().any()].tolist()
        raise ValueError(f"Non-numeric immune features: {', '.join(bad)}")
    if len(np.unique(y)) != 2:
        raise ValueError("Both High_CP and Low_CP samples are required.")
    return features, y, clean


def model_set() -> dict[str, object]:
    def scaled(model: object) -> Pipeline:
        return Pipeline([("scale", StandardScaler()), ("model", model)])

    return {
        "LASSO": scaled(LogisticRegression(
            penalty="l1", solver="liblinear", C=0.5, class_weight="balanced", random_state=42
        )),
        "Elastic Net": scaled(LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=0.5, C=0.5,
            class_weight="balanced", max_iter=5000, random_state=42
        )),
        "Ridge": scaled(RidgeClassifier(alpha=1.0, class_weight="balanced")),
        "Linear SVM": scaled(SVC(
            kernel="linear", C=0.1, class_weight="balanced", probability=True, random_state=42
        )),
        "Extra Trees": ExtraTreesClassifier(
            n_estimators=500, max_depth=2, min_samples_leaf=2,
            class_weight="balanced", random_state=42, n_jobs=-1
        ),
        "Random Forest": RandomForestClassifier(
            n_estimators=500, max_depth=2, min_samples_leaf=2,
            class_weight="balanced", random_state=42, n_jobs=-1
        ),
        "Naive Bayes": scaled(GaussianNB()),
    }


def prediction_scores(model: object, features: pd.DataFrame, y: np.ndarray, cv: LeaveOneOut) -> np.ndarray:
    if hasattr(model, "predict_proba"):
        return cross_val_predict(model, features, y, cv=cv, method="predict_proba", n_jobs=-1)[:, 1]
    return cross_val_predict(model, features, y, cv=cv, method="decision_function", n_jobs=-1)


def benchmark(features: pd.DataFrame, y: np.ndarray) -> pd.DataFrame:
    cv = LeaveOneOut()
    rows = []
    for name, model in model_set().items():
        predictions = cross_val_predict(model, features, y, cv=cv, method="predict", n_jobs=-1)
        scores = prediction_scores(model, features, y, cv)
        matrix = confusion_matrix(y, predictions, labels=[0, 1])
        rows.append({
            "Model": name,
            "LOO_Accuracy": accuracy_score(y, predictions),
            "LOO_Balanced_Accuracy": balanced_accuracy_score(y, predictions),
            "LOO_ROC_AUC": roc_auc_score(y, scores),
            "True_Low": int(matrix[0, 0]),
            "False_High": int(matrix[0, 1]),
            "False_Low": int(matrix[1, 0]),
            "True_High": int(matrix[1, 1]),
        })
    result = pd.DataFrame(rows)
    return result.sort_values(["LOO_Balanced_Accuracy", "LOO_ROC_AUC"], ascending=False)


def plot_performance(performance: pd.DataFrame, output_dir: Path) -> None:
    sns.set_theme(style="ticks")
    plot_data = performance.sort_values("LOO_Balanced_Accuracy")
    figure, axis = plt.subplots(figsize=(7, 4.8), dpi=300)
    sns.barplot(
        data=plot_data, x="LOO_Balanced_Accuracy", y="Model",
        color="#4A90E2", edgecolor="black", ax=axis
    )
    axis.axvline(0.5, color="#B2182B", linestyle="--", linewidth=1, label="Chance level")
    axis.set(title="Model comparison with leave-one-out validation", xlabel="Balanced accuracy", ylabel="")
    axis.set_xlim(0, 1.02)
    axis.legend(loc="lower right", frameon=False)
    sns.despine(ax=axis)
    figure.tight_layout()
    figure.savefig(output_dir / "ML_Model_Comparison.pdf", bbox_inches="tight")
    figure.savefig(output_dir / "ML_Model_Comparison.png", dpi=300, bbox_inches="tight")
    plt.close(figure)


def fit_lasso(features: pd.DataFrame, y: np.ndarray) -> pd.DataFrame:
    model = clone(model_set()["LASSO"])
    model.fit(features, y)
    weights = model.named_steps["model"].coef_[0]
    ranking = pd.DataFrame({"Feature": features.columns, "Weight": weights})
    ranking["Importance"] = ranking["Weight"].abs()
    return ranking.sort_values("Importance", ascending=False).reset_index(drop=True)


def adjust_bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def gene_immune_analysis(
    immune_data: pd.DataFrame,
    features: pd.DataFrame,
    gene_path: Path,
    output_dir: Path,
) -> None:
    genes = pd.read_csv(gene_path)
    if "Sample_ID" not in genes.columns:
        raise ValueError("The gene file needs a Sample_ID column.")
    merged = immune_data[["Sample_ID"] + features.columns.tolist()].merge(genes, on="Sample_ID")
    gene_names = [column for column in genes.columns if column != "Sample_ID"]
    correlations = pd.DataFrame(index=features.columns, columns=gene_names, dtype=float)
    p_values = correlations.copy()
    for immune_feature in features.columns:
        for gene in gene_names:
            correlation, p_value = spearmanr(merged[immune_feature], merged[gene])
            correlations.loc[immune_feature, gene] = correlation
            p_values.loc[immune_feature, gene] = p_value
    q_values = pd.DataFrame(
        adjust_bh(p_values.to_numpy().ravel()).reshape(p_values.shape),
        index=p_values.index,
        columns=p_values.columns,
    )

    correlations.to_csv(output_dir / "Gene_Immune_Spearman_Correlation.csv")
    p_values.to_csv(output_dir / "Gene_Immune_P_Values.csv")
    q_values.to_csv(output_dir / "Gene_Immune_Adjusted_P_Values.csv")

    figure, axis = plt.subplots(figsize=(11, 8), dpi=300)
    sns.heatmap(
        correlations, cmap="RdBu_r", center=0, vmin=-0.8, vmax=0.8,
        linewidths=0.4, linecolor="white",
        cbar_kws={"label": "Spearman correlation"}, ax=axis
    )
    for row in range(q_values.shape[0]):
        for column in range(q_values.shape[1]):
            q_value = q_values.iat[row, column]
            label = "**" if q_value < 0.01 else ("*" if q_value < 0.05 else "")
            if label:
                axis.text(column + 0.5, row + 0.5, label, ha="center", va="center", fontsize=9)
    axis.set(title="Hub gene and immune-cell associations", xlabel="Hub gene", ylabel="Immune cell")
    figure.tight_layout()
    figure.savefig(output_dir / "Gene_Immune_Correlation_Heatmap.pdf", bbox_inches="tight")
    figure.savefig(output_dir / "Gene_Immune_Correlation_Heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    features, y, immune_data = load_data(args.cibersort)
    performance = benchmark(features, y)
    performance.to_csv(args.output_dir / "ML_Model_Performance.csv", index=False)
    plot_performance(performance, args.output_dir)

    lasso_ranking = fit_lasso(features, y)
    lasso_ranking.to_csv(args.output_dir / "LASSO_Feature_Ranking.csv", index=False)
    gene_immune_analysis(immune_data, features, args.genes, args.output_dir)

    class_counts = pd.Series(y).value_counts().sort_index()
    print(f"Samples: {len(y)} (Low_CP={class_counts.get(0, 0)}, High_CP={class_counts.get(1, 0)})")
    print(performance.to_string(index=False))
    print(f"Output: {args.output_dir}")


if __name__ == "__main__":
    main()
