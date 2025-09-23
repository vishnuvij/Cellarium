from pathlib import Path
import scanpy as sc
import anndata as ad

DATA_DIR = Path("data")

def load_dataset(dataset_id: str) -> ad.AnnData:
    """
    Load a dataset by ID from the data/ folder.
    """
    path = DATA_DIR / f"{dataset_id}.h5ad"
    if not path.exists():
        raise FileNotFoundError(f"Dataset {dataset_id} not found")
    adata = sc.read_h5ad(path)
    return adata
