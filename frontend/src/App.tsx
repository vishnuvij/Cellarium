import React, { useEffect, useState, useRef } from "react";
import Plot from "react-plotly.js";
import axios from "axios";

const API_BASE = "http://localhost:8000";

const App = () => {
  const [datasets, setDatasets] = useState<{ id: string; filename: string }[]>([]);
  const [selectedDataset, setSelectedDataset] = useState<string>(""); // active dataset
  const [pendingDataset, setPendingDataset] = useState<string>("");   // selection waiting for OK
  const [points, setPoints] = useState<number[][]>([]);
  const [cellIds, setCellIds] = useState<string[]>([]);
  const [clusters, setClusters] = useState<string[]>([]);
  const [pseudotime, setPseudotime] = useState<number[]>([]);
  const [loading, setLoading] = useState(false);
  const [highlightCell, setHighlightCell] = useState<string>("");
  const [viewMode, setViewMode] = useState<"cluster" | "trajectory">("cluster");
  const [darkMode, setDarkMode] = useState(true);
  const [timeCutoff, setTimeCutoff] = useState<number>(1.0);

  const plotRef = useRef<any>(null);

  // --- Fetch available datasets ---
  const fetchDatasets = async () => {
    try {
      const res = await axios.get(`${API_BASE}/datasets`);
      setDatasets(res.data);
    } catch (err) {
      console.error("❌ Failed to fetch datasets:", err);
    }
  };

  // --- Fetch embedding + trajectory ---
  const fetchData = async (datasetId: string) => {
    if (!datasetId) return;
    setLoading(true);
    try {
      const [embedRes, trajRes] = await Promise.all([
        axios.get(`${API_BASE}/embedding/${datasetId}`),
        axios.get(`${API_BASE}/trajectory/${datasetId}`),
      ]);
      setPoints(embedRes.data.umap || []);
      setCellIds(embedRes.data.cell_ids || []);
      setClusters(embedRes.data.clusters || []);
      setPseudotime(trajRes.data.pseudotime || []);
    } catch (err) {
      console.error("❌ Failed to fetch embedding/trajectory:", err);
    } finally {
      setLoading(false);
    }
  };

  // --- File upload ---
  const handleUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    if (!e.target.files || e.target.files.length === 0) return;
    const file = e.target.files[0];
    const formData = new FormData();
    formData.append("file", file);
    try {
      await axios.post(`${API_BASE}/upload`, formData, {
        headers: { "Content-Type": "multipart/form-data" },
      });
      alert(`✅ Uploaded: ${file.name}`);
      fetchDatasets();
    } catch (err) {
      console.error("❌ Upload failed:", err);
      alert("Upload failed!");
    }
  };

  // --- Delete dataset ---
  const handleDelete = async (id: string) => {
    if (!window.confirm(`Delete dataset ${id}?`)) return;
    try {
      await axios.delete(`${API_BASE}/dataset/${id}`);
      alert(`🗑 Deleted: ${id}`);
      fetchDatasets();
      if (selectedDataset === id) {
        // Auto-clear if current dataset is deleted
        setSelectedDataset("");
        setPendingDataset("");
        setPoints([]);
        setClusters([]);
        setCellIds([]);
        setPseudotime([]);
      }
    } catch (err) {
      console.error("❌ Failed to delete:", err);
      alert("Delete failed!");
    }
  };

  useEffect(() => {
    fetchDatasets();
  }, []);

  // --- Startup view (no dataset loaded yet) ---
  if (!selectedDataset) {
    return (
      <div
        style={{
          width: "100%",
          height: "100vh",
          backgroundColor: darkMode ? "#111" : "#fff",
          color: darkMode ? "#eee" : "#111",
          display: "flex",
          flexDirection: "column",
          alignItems: "center",
          justifyContent: "center",
        }}
      >
        <h2>📂 Please upload or select a dataset to begin</h2>
        <input type="file" accept=".h5ad" onChange={handleUpload} />
        <div style={{ marginTop: "1rem" }}>
          {datasets.length > 0 && (
            <>
              <select
                value={pendingDataset}
                onChange={(e) => setPendingDataset(e.target.value)}
                style={{ padding: "0.5rem", marginRight: "1rem" }}
              >
                <option value="">-- Select dataset --</option>
                {datasets.map((d) => (
                  <option key={d.id} value={d.id}>
                    {d.filename}
                  </option>
                ))}
              </select>
              <button
                onClick={() => {
                  setSelectedDataset(pendingDataset);
                  fetchData(pendingDataset);
                }}
                disabled={!pendingDataset}
                style={{ padding: "0.5rem 1rem" }}
              >
                ✅ Load
              </button>
            </>
          )}
        </div>
      </div>
    );
  }

  // --- Loading state ---
  if (loading) return <div style={{ padding: "2rem" }}>⏳ Loading dataset...</div>;

  if (!points.length) {
    return <div style={{ padding: "2rem" }}>⚠️ No points found for {selectedDataset}</div>;
  }

  // --- Cluster grouping ---
  const clusterGroups: { [key: string]: { x: number[]; y: number[]; z: number[]; text: string[] } } = {};
  points.forEach((p, i) => {
    const cluster = clusters[i] ?? "Unassigned";
    if (!clusterGroups[cluster]) clusterGroups[cluster] = { x: [], y: [], z: [], text: [] };
    clusterGroups[cluster].x.push(p[0]);
    clusterGroups[cluster].y.push(p[1]);
    clusterGroups[cluster].z.push(p[2] ?? 0);
    clusterGroups[cluster].text.push(`Cell: ${cellIds[i]}<br>Cluster: ${cluster}`);
  });

  const clusterTraces = Object.keys(clusterGroups).map((cluster) => ({
    x: clusterGroups[cluster].x,
    y: clusterGroups[cluster].y,
    z: clusterGroups[cluster].z,
    mode: "markers",
    type: "scatter3d" as const,
    name: `Cluster ${cluster} (n=${clusterGroups[cluster].x.length})`,
    text: clusterGroups[cluster].text,
    hoverinfo: "text",
    marker: { size: 3 },
  }));

  // --- Trajectory trace ---
  const filteredIndices = pseudotime
    .map((pt, i) => ({ pt, i }))
    .filter(({ pt }) => pt <= timeCutoff)
    .map(({ i }) => i);

  const trajectoryTrace = {
    x: filteredIndices.map((i) => points[i][0]),
    y: filteredIndices.map((i) => points[i][1]),
    z: filteredIndices.map((i) => points[i][2] ?? 0),
    mode: "markers",
    type: "scatter3d" as const,
    name: "Trajectory (pseudotime)",
    text: filteredIndices.map((i) => `Cell: ${cellIds[i]}<br>Pseudotime: ${pseudotime[i]?.toFixed(3)}`),
    hoverinfo: "text",
    marker: {
      size: 3,
      color: filteredIndices.map((i) => pseudotime[i]),
      colorscale: "Plasma",
      colorbar: { title: "Pseudotime" },
    },
  };

  const traces = viewMode === "cluster" ? [...clusterTraces] : [trajectoryTrace];

  if (highlightCell && cellIds.includes(highlightCell)) {
    const i = cellIds.indexOf(highlightCell);
    traces.push({
      x: [points[i][0]],
      y: [points[i][1]],
      z: [points[i][2] ?? 0],
      mode: "markers+text",
      type: "scatter3d" as const,
      name: "Highlighted Cell",
      text: [`🔎 ${highlightCell}`],
      hoverinfo: "text",
      marker: { size: 8, color: "red", symbol: "diamond" },
    });
  }

  // --- Main UI ---
  return (
    <div
      style={{
        width: "100%",
        height: "100vh",
        backgroundColor: darkMode ? "#111" : "#fff",
        color: darkMode ? "#eee" : "#111",
      }}
    >
      <h2 style={{ textAlign: "center" }}>Cellarium 3D Embedding</h2>

      {/* Dataset controls */}
      <div style={{ textAlign: "center", marginBottom: "1rem" }}>
        <label style={{ marginRight: "1rem" }}>Dataset:</label>
        <select
          value={selectedDataset}
          onChange={(e) => {
            setSelectedDataset(e.target.value);
            fetchData(e.target.value);
          }}
          style={{ padding: "0.5rem", marginRight: "1rem" }}
        >
          {datasets.map((d) => (
            <option key={d.id} value={d.id}>
              {d.filename}
            </option>
          ))}
        </select>
        <input type="file" accept=".h5ad" onChange={handleUpload} style={{ marginRight: "1rem" }} />
        <button onClick={() => handleDelete(selectedDataset)} style={{ padding: "0.5rem 1rem" }}>
          🗑 Delete
        </button>
      </div>

      {/* Controls */}
      <div style={{ textAlign: "center", marginBottom: "1rem" }}>
        <input
          type="text"
          placeholder="Enter Cell ID..."
          value={highlightCell}
          onChange={(e) => setHighlightCell(e.target.value)}
          style={{ padding: "0.5rem", width: "250px", marginRight: "1rem" }}
        />
        <button onClick={() => setHighlightCell("")} style={{ padding: "0.5rem 1rem", marginRight: "1rem" }}>
          Reset View
        </button>
        <button
          onClick={() => setViewMode(viewMode === "cluster" ? "trajectory" : "cluster")}
          style={{ padding: "0.5rem 1rem", marginRight: "1rem" }}
        >
          {viewMode === "cluster" ? "Switch to Trajectory View" : "Switch to Cluster View"}
        </button>
        <button onClick={() => setDarkMode(!darkMode)} style={{ padding: "0.5rem 1rem", marginLeft: "1rem" }}>
          {darkMode ? "Switch to Light Mode" : "Switch to Dark Mode"}
        </button>
      </div>

      {/* Slider */}
      {viewMode === "trajectory" && (
        <div style={{ textAlign: "center", marginBottom: "1rem" }}>
          <label>Pseudotime cutoff: {timeCutoff.toFixed(2)}</label>
          <input
            type="range"
            min="0"
            max="1"
            step="0.01"
            value={timeCutoff}
            onChange={(e) => setTimeCutoff(parseFloat(e.target.value))}
            style={{ width: "60%", marginLeft: "1rem" }}
          />
        </div>
      )}

      {/* Plot */}
      <Plot
        ref={plotRef}
        data={traces}
        layout={{
          autosize: true,
          height: 700,
          paper_bgcolor: darkMode ? "#111" : "#fff",
          plot_bgcolor: darkMode ? "#111" : "#fff",
          font: { color: darkMode ? "#eee" : "#111" },
          legend: { title: { text: viewMode === "cluster" ? "Clusters" : "Trajectory" } },
          scene: {
            xaxis: { title: "UMAP-1" },
            yaxis: { title: "UMAP-2" },
            zaxis: { title: "UMAP-3" },
          },
        }}
        style={{ width: "100%", height: "100%" }}
        config={{ displayModeBar: true, scrollZoom: true }}
      />
    </div>
  );
};

export default App;
