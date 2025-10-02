# Cellarium

# Cellarium 

A lightweight web application for visualizing single-cell `.h5ad` datasets with **Scanpy** (backend) and **Plotly + React** (frontend).  
Upload your dataset, compute embeddings/trajectories, and interactively explore them in 3D.

---

##  Features
- Upload `.h5ad` datasets via the UI  
- Automatic preprocessing (PCA, neighbors, UMAP, Leiden clustering)  
- Interactive 3D plots with cluster and pseudotime views  
- Dark/Light mode toggle  
- Pseudotime slider for trajectory exploration  

---

##  Installation

### Requirements
- Docker & Docker Compose installed

### Clone Repository
```bash
git clone https://github.com/vishnuvij/Cellarium.git
cd Cellarium
```

### Build & Run
```bash
docker compose up --build
```

- Backend: http://localhost:8000  
- Frontend: http://localhost:5173  

---

##  Usage

1. Open **http://localhost:5173** in your browser.  
2. Upload a `.h5ad` dataset (e.g. PBMC3k).  
3. After upload, select your dataset from the dropdown.  
4. Explore:
   - **Cluster View** – Leiden clusters in 3D  
   - **Trajectory View** – pseudotime progression with a slider  
   - Highlight specific cells by ID  
   - Switch between dark/light mode  

---

##  Development Notes
- Backend: [FastAPI](https://fastapi.tiangolo.com/) + [Scanpy](https://scanpy.readthedocs.io/)  
- Frontend: [React](https://react.dev/) + [Plotly.js](https://plotly.com/javascript/)  
- Data directory (`/data`) is ignored via `.gitignore` – large datasets stay local, not in GitHub.  

---

##  Known Issues
- `leidenalg` must be installed for clustering:
  ```bash
  pip install leidenalg
  ```
  (already handled inside Dockerfile).  
- Large datasets may take time to embed.

---

##  License
MIT
