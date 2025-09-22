from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path

app = FastAPI(title="Cellarium API")

# CORS for future frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

DATA_DIR = Path("data")  # maps to /app/data inside the container via docker-compose

@app.get("/health")
def health():
    return {"status": "ok"}

@app.get("/")
def root():
    return {"message": "Cellarium backend is up. See /docs for the API."}

@app.get("/datasets")
def list_datasets():
    """
    Return a simple list of .h5ad files found in DATA_DIR.
    If empty, returns [].
    """
    DATA_DIR.mkdir(exist_ok=True)
    files = sorted(DATA_DIR.glob("*.h5ad"))
    return [{"id": f.stem, "filename": f.name, "size_bytes": f.stat().st_size} for f in files]
