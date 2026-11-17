# ============================================================
# run_dl_models_simulation_updated.py
#
# Simulation DL script aligned with the thesis main analysis:
#   - LSTM-NARX, TCN, ODE-LSTM
#   - train / validation / test patient-level split from R
#   - target = last evaluable observation per patient within each split
#   - metrics on original scale
#   - metrics by Overall / Control / Statin
# ============================================================

import os
import sys
import gc
import argparse
import copy
import random
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

os.environ["TF_DETERMINISTIC_OPS"] = "1"

import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, LSTM, Dense, Concatenate, Dropout, Conv1D, GlobalAveragePooling1D
from tensorflow.keras.callbacks import EarlyStopping

try:
    from tcn import TCN
except Exception as e:
    TCN = None
    TCN_IMPORT_ERROR = e

import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader


OUTCOMES = ["HbA1c", "Glu_AC", "LDL", "TG"]

OUTCOME_LABELS = {
    "HbA1c": "HbA1c",
    "Glu_AC": "Glu-AC",
    "LDL": "LDL",
    "TG": "TG"
}

GROUP_LABELS = {
    "Overall": "整體",
    "Control": "Control 組",
    "Statin": "Statin 組"
}

MODEL_ORDER = ["LSTM-NARX", "TCN", "ODE-LSTM"]
AGE_VAR = "Age_baseline"
BMI_VAR = "BMI"
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Path to one simulated dataset CSV")
    ap.add_argument("--output", required=True, help="Path to output metrics CSV")
    ap.add_argument("--seed", type=int, default=5014, help="Seed for model initialization/training")
    ap.add_argument("--max_len", type=int, default=7)
    ap.add_argument("--epochs", type=int, default=50)
    ap.add_argument("--batch_size", type=int, default=32)
    ap.add_argument("--patience", type=int, default=5)
    ap.add_argument("--scenario", required=True)
    ap.add_argument("--run", required=True, type=int)
    return ap.parse_args()


def set_seed(seed: int):
    np.random.seed(seed)
    random.seed(seed)
    tf.random.set_seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    try:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    except Exception:
        pass


def standardize_outcome_names(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    if "HbA1C" in df.columns and "HbA1c" not in df.columns:
        rename_map["HbA1C"] = "HbA1c"
    if "Glucose_AC" in df.columns and "Glu_AC" not in df.columns:
        rename_map["Glucose_AC"] = "Glu_AC"
    if "Triglycerides" in df.columns and "TG" not in df.columns:
        rename_map["Triglycerides"] = "TG"
    if rename_map:
        df = df.rename(columns=rename_map)
    return df


def make_log_col(outcome: str) -> str:
    return f"log_{outcome}"


def make_z_col(col: str) -> str:
    return f"{col}_z"


def regression_metrics(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)

    keep = np.isfinite(y_true) & np.isfinite(y_pred) & (y_true > 0)
    yt = y_true[keep]
    yp = y_pred[keep]

    n_test = len(yt)
    if n_test == 0:
        return np.nan, np.nan, np.nan, np.nan, 0

    err = yp - yt
    return (
        float(np.sqrt(np.mean(err ** 2))),
        float(np.mean(np.abs(err))),
        float(np.mean(np.abs(err / yt)) * 100),
        float(np.mean(err)),
        int(n_test)
    )


def calc_metrics_by_group(pred_df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    group_defs = [("Overall", pred_df)]
    for g in ["Control", "Statin"]:
        group_defs.append((g, pred_df[pred_df["Group"].astype(str) == g]))

    for group_eval, sub in group_defs:
        rmse, mae, mape, bias, n_test = regression_metrics(sub["y_true"], sub["y_pred"])
        rows.append({
            "Group_eval": group_eval,
            "Group_label": GROUP_LABELS[group_eval],
            "RMSE": rmse,
            "MAE": mae,
            "MAPE": mape,
            "Bias": bias,
            "n_test": n_test
        })

    return pd.DataFrame(rows)


def pad_left_1d(seq, max_len, pad_value=0.0):
    seq = np.asarray(seq, dtype=float)
    T = len(seq)
    if T >= max_len:
        return seq[-max_len:]
    pad = np.full((max_len - T,), pad_value, dtype=float)
    return np.concatenate([pad, seq])


def preprocess_for_model(df: pd.DataFrame) -> pd.DataFrame:
    df = standardize_outcome_names(df)
    df = df.copy().sort_values(["ID", "Time_in_Years"]).reset_index(drop=True)

    required_cols = ["ID", "Group", "split", "Time_in_Years", AGE_VAR, "Sex", BMI_VAR] + OUTCOMES
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    df["ID"] = df["ID"].astype(str)
    df["Group"] = df["Group"].astype(str)
    df["split"] = df["split"].astype(str)
    df["Time_in_Years"] = pd.to_numeric(df["Time_in_Years"], errors="coerce")

    # Sex coding
    if df["Sex"].dtype == object:
        s = df["Sex"].astype(str)
        unique_s = set(s.dropna().unique())
        if unique_s.issubset({"0", "1"}):
            df["Sex_num"] = s.astype(int)
        elif unique_s.issubset({"Female", "Male"}):
            df["Sex_num"] = (s == "Male").astype(int)
        elif unique_s.issubset({"F", "M"}):
            df["Sex_num"] = (s == "M").astype(int)
        else:
            df["Sex_num"] = pd.factorize(s)[0]
    else:
        df["Sex_num"] = df["Sex"].astype(int)

    df["Group_Statin"] = (df["Group"] == "Statin").astype(int)

    for col in [AGE_VAR, BMI_VAR]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Irregular-time feature
    df["Delta_Time"] = (
        df.groupby("ID")["Time_in_Years"]
        .diff()
        .fillna(0.0)
        .clip(lower=0.0)
    )

    for outcome in OUTCOMES:
        df[outcome] = pd.to_numeric(df[outcome], errors="coerce")
        if (df[outcome] <= 0).any():
            raise ValueError(f"{outcome} has non-positive values; log-transform is invalid.")
        df[make_log_col(outcome)] = np.log(df[outcome].astype(float))

    return df


def validate_patient_split(df: pd.DataFrame) -> pd.DataFrame:
    valid_values = {"train", "validation", "test"}
    observed_values = set(df["split"].astype(str).unique())

    if not observed_values.issubset(valid_values):
        raise ValueError(f"Invalid split values: {observed_values}. Expected only: train, validation, test.")

    split_per_id = df[["ID", "split"]].drop_duplicates().groupby("ID")["split"].nunique()
    if (split_per_id > 1).any():
        raise ValueError("At least one patient ID appears in more than one split.")

    return df


def fit_scaler_from_train(df_train: pd.DataFrame, outcome: str) -> dict:
    continuous_cols = [
        make_log_col(outcome),
        "Time_in_Years",
        "Delta_Time",
        AGE_VAR,
        BMI_VAR
    ]

    scaler = {}
    for col in continuous_cols:
        mu = float(df_train[col].mean())
        sd = float(df_train[col].std())
        if sd == 0 or np.isnan(sd):
            sd = 1.0
        scaler[col] = {"mean": mu, "sd": sd}
    return scaler


def apply_scaler(df: pd.DataFrame, scaler: dict, outcome: str) -> pd.DataFrame:
    df = df.copy()
    continuous_cols = [
        make_log_col(outcome),
        "Time_in_Years",
        "Delta_Time",
        AGE_VAR,
        BMI_VAR
    ]
    for col in continuous_cols:
        mu = scaler[col]["mean"]
        sd = scaler[col]["sd"]
        df[make_z_col(col)] = (df[col] - mu) / sd
    return df


def inverse_transform_log_outcome(y_std, scaler: dict, outcome: str):
    y_col = make_log_col(outcome)
    mu = scaler[y_col]["mean"]
    sd = scaler[y_col]["sd"]
    y_log = np.asarray(y_std) * sd + mu
    return np.exp(y_log)


def build_sequence_dataset(df_scaled: pd.DataFrame, split_name: str, outcome: str, max_len: int = 7, min_history: int = 1):
    y_z_col = make_z_col(make_log_col(outcome))
    dt_z_col = make_z_col("Delta_Time")
    time_z_col = make_z_col("Time_in_Years")
    age_z_col = make_z_col(AGE_VAR)
    bmi_z_col = make_z_col(BMI_VAR)

    exog_cols = [time_z_col, age_z_col, "Sex_num", bmi_z_col, "Group_Statin"]

    df_split = (
        df_scaled[df_scaled["split"].astype(str) == split_name]
        .copy()
        .sort_values(["ID", "Time_in_Years"])
    )

    X_seq, X_dt, X_exog, y, mask, meta = [], [], [], [], [], []

    for pid, g in df_split.groupby("ID", sort=False):
        g = g.sort_values("Time_in_Years").reset_index(drop=True)

        if len(g) <= min_history:
            continue

        needed_target = [y_z_col] + exog_cols + ["Time_in_Years"]
        valid_target_mask = g[needed_target].notna().all(axis=1).values
        valid_target_pos = np.where(valid_target_mask)[0]

        if len(valid_target_pos) == 0:
            continue

        target_pos = int(valid_target_pos[-1])
        if target_pos < min_history:
            continue

        hist = g.iloc[:target_pos].copy()
        target = g.iloc[target_pos].copy()

        hist = hist.dropna(subset=[y_z_col, dt_z_col, "Delta_Time"])
        if len(hist) < min_history:
            continue

        hist_y = hist[y_z_col].values.astype(float)
        hist_dt_z = hist[dt_z_col].values.astype(float)
        hist_dt_raw = hist["Delta_Time"].values.astype(float)

        valid_len = min(len(hist), max_len)
        mask_vec = np.concatenate([
            np.zeros(max_len - valid_len, dtype=float),
            np.ones(valid_len, dtype=float)
        ])

        seq_features = np.column_stack([
            pad_left_1d(hist_y, max_len, pad_value=0.0),
            pad_left_1d(hist_dt_z, max_len, pad_value=0.0),
            mask_vec
        ])

        dt_features = pad_left_1d(hist_dt_raw, max_len, pad_value=0.0).reshape(-1, 1)

        X_seq.append(seq_features)
        X_dt.append(dt_features)
        X_exog.append(target[exog_cols].values.astype(float))
        y.append(float(target[y_z_col]))
        mask.append(mask_vec)

        meta.append({
            "row_id": int(target["row_id"]) if "row_id" in target.index else -1,
            "ID": str(target["ID"]),
            "Group": str(target["Group"]),
            "Time_in_Years": float(target["Time_in_Years"]),
            "target_definition": "last_simulated_test_target_with_history"
        })

    if len(y) == 0:
        return (
            np.empty((0, max_len, 3), dtype=np.float32),
            np.empty((0, max_len, 1), dtype=np.float32),
            np.empty((0, len(exog_cols)), dtype=np.float32),
            np.empty((0,), dtype=np.float32),
            np.empty((0, max_len), dtype=np.float32),
            pd.DataFrame()
        )

    return (
        np.asarray(X_seq, dtype=np.float32),
        np.asarray(X_dt, dtype=np.float32),
        np.asarray(X_exog, dtype=np.float32),
        np.asarray(y, dtype=np.float32),
        np.asarray(mask, dtype=np.float32),
        pd.DataFrame(meta)
    )


def build_lstm_narx(seq_shape, exog_shape):
    """
    LSTM-NARX:
    - autoregressive branch: past standardized log outcome and delta-time sequence
    - exogenous branch: target-time Group, time, age, sex, BMI
    This structure follows the main thesis DL notebook.
    """
    seq_in = Input(shape=seq_shape, name="seq_in")
    exog_in = Input(shape=exog_shape, name="exog_in")

    x = LSTM(64, name="lstm_ar_branch")(seq_in)
    x = Concatenate(name="concat_ar_exog")([x, exog_in])
    x = Dense(32, activation="relu")(x)
    x = Dropout(0.20)(x)
    out = Dense(1, name="output")(x)

    model = Model([seq_in, exog_in], out, name="LSTM_NARX")
    model.compile(optimizer="adam", loss="mse")
    return model


def build_tcn_dual_input(seq_shape, static_shape):
    """
    TCN:
    - uses keras-tcn when available
    - falls back to causal Conv1D stack when keras-tcn is unavailable,
      following the main thesis DL notebook.
    """
    seq_in = Input(shape=seq_shape, name="seq_in")
    static_in = Input(shape=static_shape, name="static_in")

    if TCN is not None:
        x = TCN(
            nb_filters=32,
            kernel_size=3,
            dilations=[1, 2, 4],
            dropout_rate=0.10,
            return_sequences=False,
            name="tcn_branch"
        )(seq_in)
    else:
        x = Conv1D(32, kernel_size=3, padding="causal", activation="relu", dilation_rate=1)(seq_in)
        x = Dropout(0.10)(x)
        x = Conv1D(32, kernel_size=3, padding="causal", activation="relu", dilation_rate=2)(x)
        x = Dropout(0.10)(x)
        x = Conv1D(32, kernel_size=3, padding="causal", activation="relu", dilation_rate=4)(x)
        x = GlobalAveragePooling1D()(x)

    x = Concatenate(name="concat_tcn_exog")([x, static_in])
    x = Dense(32, activation="relu")(x)
    x = Dropout(0.20)(x)
    out = Dense(1, name="output")(x)

    model = Model([seq_in, static_in], out, name="TCN")
    model.compile(optimizer="adam", loss="mse")
    return model


class ODEFunc(nn.Module):
    def __init__(self, hidden_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, hidden_dim)
        )

    def forward(self, h):
        return self.net(h)


class ODELSTMModel(nn.Module):
    """
    ODE-LSTM for simulation:
    - sequence features: [standardized log outcome, standardized delta time, valid mask]
    - uses each patient's own dt through an Euler ODE hidden transition
    - avoids batch-average dt, which would blur individual irregular-time information
    """
    def __init__(self, input_dim=3, exog_dim=5, hidden_dim=64):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.odefunc = ODEFunc(hidden_dim)
        self.lstm_cell = nn.LSTMCell(input_dim, hidden_dim)
        self.fc = nn.Sequential(
            nn.Linear(hidden_dim + exog_dim, 64),
            nn.ReLU(),
            nn.Dropout(0.20),
            nn.Linear(64, 1)
        )

    def forward(self, x_seq, x_dt, x_exog):
        B, T, _ = x_seq.shape
        h = torch.zeros(B, self.hidden_dim, device=x_seq.device)
        c = torch.zeros(B, self.hidden_dim, device=x_seq.device)

        for t in range(T):
            x_t = x_seq[:, t, :]
            dt_t = x_dt[:, t, 0].clamp(min=0.0, max=5.0).unsqueeze(1)
            mask_t = x_seq[:, t, 2].unsqueeze(1)

            h_ode = h + dt_t * self.odefunc(h)
            h_new, c_new = self.lstm_cell(x_t, (h_ode, c))

            h = mask_t * h_new + (1 - mask_t) * h
            c = mask_t * c_new + (1 - mask_t) * c

        return self.fc(torch.cat([h, x_exog], dim=1)).squeeze(-1)


def train_one_odelstm(model, train_loader, val_loader, epochs=50, patience=5, lr=1e-3):
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss()

    best_model = None
    best_val = np.inf
    wait = 0

    for _ in range(epochs):
        model.train()
        for x_seq, x_dt, x_exog, y in train_loader:
            x_seq = x_seq.to(DEVICE)
            x_dt = x_dt.to(DEVICE)
            x_exog = x_exog.to(DEVICE)
            y = y.to(DEVICE)

            optimizer.zero_grad()
            pred = model(x_seq, x_dt, x_exog)
            loss = criterion(pred, y)
            loss.backward()
            optimizer.step()

        model.eval()
        val_losses = []
        with torch.no_grad():
            for x_seq, x_dt, x_exog, y in val_loader:
                x_seq = x_seq.to(DEVICE)
                x_dt = x_dt.to(DEVICE)
                x_exog = x_exog.to(DEVICE)
                y = y.to(DEVICE)

                pred = model(x_seq, x_dt, x_exog)
                val_losses.append(criterion(pred, y).item())

        val_loss = float(np.mean(val_losses)) if len(val_losses) > 0 else np.inf

        if val_loss < best_val:
            best_val = val_loss
            best_model = copy.deepcopy(model.state_dict())
            wait = 0
        else:
            wait += 1
            if wait >= patience:
                break

    if best_model is not None:
        model.load_state_dict(best_model)

    return model


def prepare_scaled_data(df_raw: pd.DataFrame, outcome: str):
    df_train_raw = df_raw[df_raw["split"] == "train"].copy()
    scaler = fit_scaler_from_train(df_train_raw, outcome)
    df_scaled = apply_scaler(df_raw, scaler, outcome)
    return df_scaled, scaler


def run_lstm_narx(df_raw, outcome, args):
    df_scaled, scaler = prepare_scaled_data(df_raw, outcome)

    Xs_tr, Xdt_tr, Xe_tr, y_tr, mask_tr, meta_tr = build_sequence_dataset(df_scaled, "train", outcome, max_len=args.max_len)
    Xs_va, Xdt_va, Xe_va, y_va, mask_va, meta_va = build_sequence_dataset(df_scaled, "validation", outcome, max_len=args.max_len)
    Xs_te, Xdt_te, Xe_te, y_te, mask_te, meta_te = build_sequence_dataset(df_scaled, "test", outcome, max_len=args.max_len)

    if len(y_tr) == 0 or len(y_va) == 0 or len(y_te) == 0:
        return pd.DataFrame()

    model = build_lstm_narx(seq_shape=Xs_tr.shape[1:], exog_shape=Xe_tr.shape[1:])
    es = EarlyStopping(monitor="val_loss", patience=args.patience, restore_best_weights=True)

    model.fit(
        [Xs_tr, Xe_tr],
        y_tr,
        validation_data=([Xs_va, Xe_va], y_va),
        epochs=args.epochs,
        batch_size=args.batch_size,
        callbacks=[es],
        verbose=0
    )

    y_pred_std = model([Xs_te, Xe_te], training=False).numpy().ravel()
    y_true = inverse_transform_log_outcome(y_te, scaler, outcome)
    y_pred = inverse_transform_log_outcome(y_pred_std, scaler, outcome)

    pred_df = meta_te.copy()
    pred_df["y_true"] = y_true
    pred_df["y_pred"] = y_pred

    tf.keras.backend.clear_session()
    return pred_df


def run_tcn(df_raw, outcome, args):
    df_scaled, scaler = prepare_scaled_data(df_raw, outcome)

    Xs_tr, Xdt_tr, Xe_tr, y_tr, mask_tr, meta_tr = build_sequence_dataset(df_scaled, "train", outcome, max_len=args.max_len)
    Xs_va, Xdt_va, Xe_va, y_va, mask_va, meta_va = build_sequence_dataset(df_scaled, "validation", outcome, max_len=args.max_len)
    Xs_te, Xdt_te, Xe_te, y_te, mask_te, meta_te = build_sequence_dataset(df_scaled, "test", outcome, max_len=args.max_len)

    if len(y_tr) == 0 or len(y_va) == 0 or len(y_te) == 0:
        return pd.DataFrame()

    model = build_tcn_dual_input(seq_shape=Xs_tr.shape[1:], static_shape=Xe_tr.shape[1:])
    es = EarlyStopping(monitor="val_loss", patience=args.patience, restore_best_weights=True)

    model.fit(
        [Xs_tr, Xe_tr],
        y_tr,
        validation_data=([Xs_va, Xe_va], y_va),
        epochs=args.epochs,
        batch_size=args.batch_size,
        callbacks=[es],
        verbose=0
    )

    y_pred_std = model([Xs_te, Xe_te], training=False).numpy().ravel()
    y_true = inverse_transform_log_outcome(y_te, scaler, outcome)
    y_pred = inverse_transform_log_outcome(y_pred_std, scaler, outcome)

    pred_df = meta_te.copy()
    pred_df["y_true"] = y_true
    pred_df["y_pred"] = y_pred

    tf.keras.backend.clear_session()
    return pred_df


def run_odelstm(df_raw, outcome, args):
    df_scaled, scaler = prepare_scaled_data(df_raw, outcome)

    Xs_tr, Xdt_tr, Xe_tr, y_tr, mask_tr, meta_tr = build_sequence_dataset(df_scaled, "train", outcome, max_len=args.max_len)
    Xs_va, Xdt_va, Xe_va, y_va, mask_va, meta_va = build_sequence_dataset(df_scaled, "validation", outcome, max_len=args.max_len)
    Xs_te, Xdt_te, Xe_te, y_te, mask_te, meta_te = build_sequence_dataset(df_scaled, "test", outcome, max_len=args.max_len)

    if len(y_tr) == 0 or len(y_va) == 0 or len(y_te) == 0:
        return pd.DataFrame()

    train_ds = TensorDataset(
        torch.from_numpy(Xs_tr).float(),
        torch.from_numpy(Xdt_tr).float(),
        torch.from_numpy(Xe_tr).float(),
        torch.from_numpy(y_tr).float()
    )
    val_ds = TensorDataset(
        torch.from_numpy(Xs_va).float(),
        torch.from_numpy(Xdt_va).float(),
        torch.from_numpy(Xe_va).float(),
        torch.from_numpy(y_va).float()
    )
    test_ds = TensorDataset(
        torch.from_numpy(Xs_te).float(),
        torch.from_numpy(Xdt_te).float(),
        torch.from_numpy(Xe_te).float(),
        torch.from_numpy(y_te).float()
    )

    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False)
    test_loader = DataLoader(test_ds, batch_size=args.batch_size, shuffle=False)

    model = ODELSTMModel(input_dim=Xs_tr.shape[2], exog_dim=Xe_tr.shape[1], hidden_dim=64).to(DEVICE)
    model = train_one_odelstm(model, train_loader, val_loader, epochs=args.epochs, patience=args.patience, lr=1e-3)

    model.eval()
    preds = []
    with torch.no_grad():
        for x_seq, x_dt, x_exog, y in test_loader:
            pred = model(x_seq.to(DEVICE), x_dt.to(DEVICE), x_exog.to(DEVICE))
            preds.append(pred.detach().cpu().numpy())

    y_pred_std = np.concatenate(preds)
    y_true = inverse_transform_log_outcome(y_te, scaler, outcome)
    y_pred = inverse_transform_log_outcome(y_pred_std, scaler, outcome)

    pred_df = meta_te.copy()
    pred_df["y_true"] = y_true
    pred_df["y_pred"] = y_pred

    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    return pred_df


def make_metric_rows(pred_df: pd.DataFrame, outcome: str, model_name: str, scenario: str, run: int) -> pd.DataFrame:
    if pred_df is None or len(pred_df) == 0:
        return pd.DataFrame([{
            "scenario": scenario,
            "run": run,
            "outcome": outcome,
            "outcome_label": OUTCOME_LABELS[outcome],
            "model": model_name,
            "Group_eval": g,
            "Group_label": GROUP_LABELS[g],
            "RMSE": np.nan,
            "MAE": np.nan,
            "MAPE": np.nan,
            "Bias": np.nan,
            "n_test": 0
        } for g in ["Overall", "Control", "Statin"]])

    metrics = calc_metrics_by_group(pred_df)
    metrics.insert(0, "model", model_name)
    metrics.insert(0, "outcome_label", OUTCOME_LABELS[outcome])
    metrics.insert(0, "outcome", outcome)
    metrics.insert(0, "run", run)
    metrics.insert(0, "scenario", scenario)
    return metrics



def validate_final_output(res_df: pd.DataFrame) -> None:
    """Require 4 outcomes × 3 models × 3 groups = 36 complete metric rows."""
    expected = pd.MultiIndex.from_product(
        [OUTCOMES, MODEL_ORDER, ["Overall", "Control", "Statin"]],
        names=["outcome", "model", "Group_eval"]
    ).to_frame(index=False)

    observed = res_df[["outcome", "model", "Group_eval"]].drop_duplicates()
    missing = expected.merge(
        observed,
        on=["outcome", "model", "Group_eval"],
        how="left",
        indicator=True
    )
    missing = missing[missing["_merge"] == "left_only"]

    metric_cols = ["RMSE", "MAE", "MAPE", "Bias"]
    invalid_rows = (
        res_df[metric_cols].isna().any(axis=1)
        | ~np.isfinite(res_df[metric_cols]).all(axis=1)
        | (res_df["n_test"] <= 0)
    )

    if len(res_df) != 36 or len(missing) > 0 or invalid_rows.any():
        raise RuntimeError(
            "Incomplete DL output: "
            f"rows={len(res_df)}, "
            f"missing_combinations={len(missing)}, "
            f"invalid_rows={int(invalid_rows.sum())}"
        )

def main():
    args = parse_args()
    set_seed(args.seed)

    outdir = os.path.dirname(args.output)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    print(f"Using device: {DEVICE}", flush=True)

    df = pd.read_csv(args.input)
    df = preprocess_for_model(df)
    df = validate_patient_split(df)

    model_runners = {
        "LSTM-NARX": run_lstm_narx,
        "TCN": run_tcn,
        "ODE-LSTM": run_odelstm
    }

    all_rows = []
    failures = []

    for outcome_index, outcome in enumerate(OUTCOMES):
        for model_index, model_name in enumerate(MODEL_ORDER):
            local_seed = args.seed + outcome_index * 100 + model_index * 10
            set_seed(local_seed)

            runner = model_runners[model_name]
            print(f"Running {model_name} - {outcome}", flush=True)

            try:
                pred_df = runner(df, outcome, args)
                if pred_df is None or len(pred_df) == 0:
                    raise ValueError("Empty prediction result.")

                rows = make_metric_rows(
                    pred_df=pred_df,
                    outcome=outcome,
                    model_name=model_name,
                    scenario=args.scenario,
                    run=args.run
                )
                all_rows.append(rows)

            except Exception as exc:
                failures.append(f"{model_name} - {outcome}: {exc}")
                print(
                    f"FAILED: {model_name} - {outcome}: {exc}",
                    file=sys.stderr,
                    flush=True
                )

            tf.keras.backend.clear_session()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            gc.collect()

    if all_rows:
        res_df = pd.concat(all_rows, ignore_index=True)
        res_df.to_csv(args.output, index=False)
    else:
        res_df = pd.DataFrame()

    if failures:
        raise RuntimeError(
            "One or more DL tasks failed:\n" + "\n".join(failures)
        )

    validate_final_output(res_df)
    print("Saved complete metrics to:", args.output, flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr, flush=True)
        sys.exit(1)
        