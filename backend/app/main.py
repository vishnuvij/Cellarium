from fastapi import FastAPI, HTTPException, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
import scanpy as sc
import numpy as np
import shutil

from app.io.loader import load_dataset

app = FastAPI(title="Cellarium API")

# --- CORS ---
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)


# --- Health ---
@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/")
def root():
    return {"message": "Cellarium backend is up. See /docs for the API."}


# --- Upload dataset ---
@app.post("/upload")
async def upload_dataset(file: UploadFile = File(...)):
    """
    Upload a new .h5ad dataset into the data/ folder.
    """
    if not file.filename.endswith(".h5ad"):
        raise HTTPException(status_code=400, detail="Only .h5ad files are supported")

    save_path = DATA_DIR / file.filename
    try:
        with open(save_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to save file: {e}")

    return {"status": "success", "filename": file.filename}


# --- Delete dataset ---
@app.delete("/dataset/{dataset_id}")
def delete_dataset(dataset_id: str):
    """
    Delete a dataset file from disk.
    """
    file_path = DATA_DIR / f"{dataset_id}.h5ad"
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="Dataset not found")

    try:
        file_path.unlink()
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to delete dataset: {e}")

    return {"status": "deleted", "dataset_id": dataset_id}


# --- List datasets ---
@app.get("/datasets")
def list_datasets():
    DATA_DIR.mkdir(exist_ok=True)
    files = sorted(DATA_DIR.glob("*.h5ad"))
    return [
        {"id": f.stem, "filename": f.name, "size_bytes": f.stat().st_size}
        for f in files
    ]


# --- Embedding route ---
@app.get("/embedding/{dataset_id}")
def embedding(dataset_id: str, n_pcs: int = 20):
    """
    Return PCA + UMAP (3D) embedding + cluster info for the dataset.
    """
    try:
        adata = load_dataset(dataset_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="Dataset not found")

    # Clean NaNs
    X = adata.X
    if isinstance(X, np.ndarray):
        X = np.nan_to_num(X, nan=0.0)
        adata.X = X

    # PCA
    if "X_pca" not in adata.obsm:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        try:
            sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
        except Exception as e:
            print("⚠️ HVG failed, skipping:", e)
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=n_pcs)

    # Neighbors
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata)

    # UMAP 3D
    if "X_umap" not in adata.obsm or adata.obsm["X_umap"].shape[1] < 3:
        sc.tl.umap(adata, n_components=3)

    # Leiden clustering
    if "leiden" not in adata.obs:
        sc.tl.leiden(adata, resolution=0.5)

    coords = adata.obsm["X_umap"][:, :3].tolist()
    cell_ids = adata.obs_names.tolist()
    clusters = adata.obs["leiden"].tolist()

    return {
        "dataset_id": dataset_id,
        "n_cells": len(cell_ids),
        "umap": coords,
        "cell_ids": cell_ids,
        "clusters": clusters,
    }


# --- Trajectory route ---
@app.get("/trajectory/{dataset_id}")
def trajectory(dataset_id: str):
    """
    Return pseudotime trajectory values for each cell.
    Currently creates simple linear pseudotime if none exists.
    """
    try:
        adata = load_dataset(dataset_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="Dataset not found")

    if "pseudotime" not in adata.obs:
        adata.obs["pseudotime"] = np.linspace(0, 1, adata.n_obs)

    pseudotime = adata.obs["pseudotime"].tolist()
    cell_ids = adata.obs_names.tolist()

    return {
        "dataset_id": dataset_id,
        "n_cells": int(adata.n_obs),
        "cell_ids": cell_ids,
        "pseudotime": pseudotime,
    }
