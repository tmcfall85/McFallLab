import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import (
    TensorDataset,
    DataLoader,
    random_split,
    Subset,
    SubsetRandomSampler,
)
from typing import Sequence, Tuple
import numpy as np

from sklearn.model_selection import StratifiedShuffleSplit

# Reproducibility
RNG_SEED = 42
np.random.seed(RNG_SEED)
torch.manual_seed(RNG_SEED)


class FCNet(nn.Module):
    def __init__(
        self,
        input_dim: int,
        hidden_sizes: Sequence[int],
        n_classes: int,
        dropout: float = 0.2,
    ):
        super().__init__()
        layers = []
        last_dim = input_dim
        for h in hidden_sizes:
            layers.append(nn.Linear(last_dim, h))
            layers.append(nn.ReLU(inplace=True))
            layers.append(nn.BatchNorm1d(h))
            layers.append(nn.Dropout(dropout))
            last_dim = h
        layers.append(nn.Linear(last_dim, n_classes))  # logits
        layers.append(nn.Softmax(1))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


class FCNetReg(nn.Module):
    def __init__(
        self,
        input_dim: int,
        hidden_sizes: Sequence[int],
        n_classes: int = 1,
        dropout: float = 0.2,
    ):
        super().__init__()
        layers = []
        last_dim = input_dim
        for h in hidden_sizes:
            layers.append(nn.Linear(last_dim, h))
            layers.append(nn.ReLU(inplace=True))
            layers.append(nn.BatchNorm1d(h))
            layers.append(nn.Dropout(dropout))
            last_dim = h
        # single output for regression
        layers.append(nn.Linear(last_dim, 1))
        self.net = nn.Sequential(*layers)
        # Softplus ensures positive outputs while remaining smooth
        self.out_act = nn.Softplus()

    def forward(self, x):
        out = self.net(x)
        out = self.out_act(out)
        return out.squeeze(-1)


def train_fc_network(
    X_np: np.ndarray,
    y_np: np.ndarray,
    hidden_sizes: Sequence[int] = (512, 128),
    n_classes: int = 4,
    epochs: int = 5,
    batch_size: int = 4,
    lr: float = 1e-3,
    weight_decay: float = 1e-5,
    val_fraction: float = 0.2,
    n_splits: int = 5,
) -> Tuple[nn.Module, dict]:
    """
    Train a fully connected network on numpy data.

    Args:
        X_np: shape (N, D) numpy array of features (float).
        y_np: shape (N,) numpy array of integer labels in [0, n_classes-1].
        hidden_sizes: sizes of hidden layers.
        n_classes: number of classes.
        epochs, batch_size, lr: training hyperparameters.
        val_fraction: fraction of data used for validation.
        device: "cpu" or "cuda" (auto-detected if None).

    Returns:
        model, history dict with train/val loss and accuracy lists.
    """
    device = "cuda" if torch.cuda.is_available() else "cpu"

    X = np.asarray(X_np, dtype=np.float32)
    if n_classes > 1:
        y = np.asarray(y_np, dtype=np.int64).ravel()
    else:
        y = np.asarray(y_np, dtype=np.float32).ravel()

    assert X.ndim == 3, "X must be 3D (N, D1, D2)"
    assert X.shape[0] == y.shape[0], "X and y must have same first dimension"

    # Flatten each sample's (D1, D2) feature map into a vector (D1*D2)
    N, D1, D2 = X.shape
    X = X.reshape(N, D1 * D2)
    assert y.ndim == 1, "y must be 1D (N,)"
    assert X.shape[0] == y.shape[0], "X and y must have same first dimension"
    if n_classes > 1:
        assert np.all(
            (y >= 0) & (y < n_classes)
        ), "labels must be between 0 and n_classes-1"
    print(np.shape(X), np.shape(y))
    # Simple standardization (per-feature)
    mean = X.mean(axis=0, keepdims=True)
    std = X.std(axis=0, keepdims=True)
    std[std == 0] = 1.0
    X = (X - mean) / std

    X_t = torch.from_numpy(X)
    y_t = torch.from_numpy(y)

    dataset = TensorDataset(X_t, y_t)

    n_val_max = np.ceil(len(dataset) * (1 / n_splits))
    n_val_min = np.floor(len(dataset) * (1 / n_splits))
    print(n_val_min, n_val_max)
    if (len(dataset) - n_val_max) % 2 == 0:
        n_val = int(n_val_max)
    else:
        n_val = int(n_val_min)
    print(len(dataset) - n_val, n_val)

    n_train = len(dataset) - n_val
    # train_ds, val_ds = random_split(dataset, [n_train, n_val])

    # Assuming you have your features X and target y as numpy arrays or pandas DataFrames/Series
    # X and y must have the same number of samples (e.g., n_samples = 200)
    # Replace n_samples with the actual number of samples in your dataset
    n_samples = len(X)
    # validation_size = 24

    # Calculate n_splits required if needed, or simply use StratifiedShuffleSplit
    # For KFold behavior (using all data points over splits), you'd need careful calculation

    # Using StratifiedShuffleSplit to get a validation set of exactly 31 samples for each split
    ss = StratifiedShuffleSplit(n_splits=n_splits, test_size=n_val)

    # Iterate through the splits
    all_preds = []
    all_preds_proba = []
    all_y_vals = []
    histories = []
    models = []

    for i, (train_index, val_index) in enumerate(ss.split(X, y)):
        print(f"Fold {i+1}:")
        print(f"  Train set size: {len(train_index)}")
        print(f"  Validation set size: {len(val_index)}")

        # Access the data for training and validation
        # X_train, X_val = X[train_index], X[val_index]
        # y_train, y_val = y[train_index], y[val_index]
        train_sampler = SubsetRandomSampler(train_index)
        validation_sampler = SubsetRandomSampler(val_index)

        train_loader = DataLoader(dataset, batch_size=batch_size, sampler=train_sampler)
        val_loader = DataLoader(
            dataset, batch_size=batch_size, sampler=validation_sampler
        )

        # train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
        # val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False)

        if n_classes == 1:

            model = FCNetReg(
                input_dim=X.shape[1], hidden_sizes=hidden_sizes, n_classes=n_classes
            ).to(device)
            criterion = nn.MSELoss()
        else:
            model = FCNet(
                input_dim=X.shape[1], hidden_sizes=hidden_sizes, n_classes=n_classes
            ).to(device)
            criterion = nn.CrossEntropyLoss()
        optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

        history = {"train_loss": [], "val_loss": [], "train_acc": [], "val_acc": []}

        for epoch in range(1, epochs + 1):
            # Training
            model.train()
            running_loss = 0.0
            correct = 0
            total = 0
            for xb, yb in train_loader:

                xb = xb.to(device)
                # if n_classes == 1:
                yb = yb.to(device)
                # else:
                #    yb = yb.to(device)
                optimizer.zero_grad()
                logits = model(xb)
                loss = criterion(logits, yb)
                loss.backward()
                optimizer.step()

                running_loss += loss.item() * xb.size(0)
                if n_classes > 1:
                    preds = logits.argmax(dim=1)
                    correct += (preds == yb).sum().item()
                total += xb.size(0)

            train_loss = running_loss / max(1, total)

            if n_classes == 1:
                train_acc = float(np.sqrt(train_loss))
            else:
                train_acc = correct / max(1, total)

            # Validation
            model.eval()
            val_loss = 0.0
            val_correct = 0
            val_total = 0
            with torch.no_grad():
                val_pred_y = []
                val_pred_proba_y = []
                val_actual_y = []
                for xb, yb in val_loader:
                    xb = xb.to(device)
                    # if n_classes == 1:
                    yb = yb.to(device)
                    # else:
                    #    yb = yb.to(device)
                    logits = model(xb)
                    loss = criterion(logits, yb)
                    val_loss += loss.item() * xb.size(0)
                    if n_classes > 1:
                        preds = logits.argmax(dim=1)
                        val_correct += (preds == yb).sum().item()
                    val_total += xb.size(0)
                    if n_classes == 1:
                        for l in logits.cpu().numpy():
                            val_pred_y.append(l)
                    else:
                        for l in logits.cpu().numpy():
                            val_pred_proba_y.append(l)
                        for p in preds.cpu().numpy():
                            val_pred_y.append(p)
                    for ay in yb.cpu().numpy():
                        val_actual_y.append(ay)

            val_loss = val_loss / max(1, val_total)

            if n_classes == 1:
                val_acc = float(np.sqrt(val_loss))
            else:
                val_acc = val_correct / max(1, total)

            history["train_loss"].append(train_loss)
            history["val_loss"].append(val_loss)
            history["train_acc"].append(train_acc)
            history["val_acc"].append(val_acc)
            histories.append(history)
            models.append(model)
            if epoch == epochs:
                print("saving preds!!!")
                all_preds.append(val_pred_y)
                all_preds_proba.append(val_pred_proba_y)
                all_y_vals.append(val_actual_y)
            # Simple progress print
            if epoch % max(1, epochs // 10) == 0 or epoch == 1 or epoch == epochs:
                print(
                    f"Epoch {epoch}/{epochs}  train_loss={train_loss:.4f} train_acc={train_acc:.3f}  val_loss={val_loss:.4f} val_acc={val_acc:.3f}"
                )

    return models, histories, all_preds, all_preds_proba, all_y_vals
