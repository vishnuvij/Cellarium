import React, { useEffect, useState, useRef } from "react";
import ReactDOM from "react-dom/client";
import Plot from "react-plotly.js";
import axios from "axios";

const App = () => {
  const [points, setPoints] = useState<number[][]>([]);
  const [cellIds, setCellIds] = useState<string[]>([]);
  const [clusters, setClusters] = useState<string[]>([]);
  const [loading, setLoading] = useState(true);
  const [highlightCell, setHighlightCell] = useState<string>("");

  const plotRef = useRef<any>(null);

  // --- Fetch data ---
  useEffect(() => {
    const fetchData = async () => {
      try {
        const res = await axios.get("http://localhost:8000/embedding/scanpy-pbmc3k");
        console.log("✅ Data fetched:", res.data);

        setPoints(res.data.umap || []);
        setCellIds(res.data.cell_ids || []);
        setClusters(res.data.clusters || []);
      } catch (err) {
        console.error("❌ Failed to fetch embedding:", err);
      } finally {
        setLoading(false);
      }
    };
    fetchData();
  }, []);

  if (loading) return <div style={{ padding: "2rem" }}>Loading 3D scatter...</div>;

  if (!points.length) return <div style={{ padding: "2rem" }}>⚠️ No points found</div>;

  // --- Group by cluster ---
  const clusterGroups: {
    [key: string]: { x: number[]; y: number[]; z: number[]; text: string[] };
  } = {};

  points.forEach((p, i) => {
    const cluster = clusters[i] ?? "Unassigned";
    if (!clusterGroups[cluster]) {
      clusterGroups[cluster] = { x: [], y: [], z: [], text: [] };
    }
    clusterGroups[cluster].x.push(p[0]);
    clusterGroups[cluster].y.push(p[1]);
    clusterGroups[cluster].z.push(p[2] ?? 0);
    clusterGroups[cluster].text.push(`Cell: ${cellIds[i]}<br>Cluster: ${cluster}`);
  });

  const clusterNames = Object.keys(clusterGroups);

  // --- Build traces ---
  const traces = clusterNames.map((cluster) => ({
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

  console.log("📊 Traces built:", traces);

  // --- Highlight cell ---
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

  // --- Reset view ---
  const resetView = () => {
    if (plotRef.current) {
      plotRef.current.react(traces, {
        scene: {
          camera: {
            center: { x: 0, y: 0, z: 0 },
            eye: { x: 1.5, y: 1.5, z: 1.5 },
            up: { x: 0, y: 0, z: 1 },
          },
          xaxis: { title: "UMAP-1" },
          yaxis: { title: "UMAP-2" },
          zaxis: { title: "UMAP-3" },
        },
      });
    }
    setHighlightCell("");
  };

  return (
    <div style={{ width: "100%", height: "100vh" }}>
      <h2 style={{ textAlign: "center" }}>Cellarium 3D Embedding</h2>

      {/* Search + Reset */}
      <div style={{ textAlign: "center", marginBottom: "1rem" }}>
        <input
          type="text"
          placeholder="Enter Cell ID..."
          value={highlightCell}
          onChange={(e) => setHighlightCell(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && setHighlightCell(highlightCell.trim())}
          style={{ padding: "0.5rem", width: "250px", marginRight: "1rem" }}
        />
        <button onClick={resetView} style={{ padding: "0.5rem 1rem" }}>
          Reset View
        </button>
      </div>

      <Plot
        ref={plotRef}
        data={traces}
        layout={{
          autosize: true,
          height: 700,
          legend: { title: { text: "Clusters" } },
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

ReactDOM.createRoot(document.getElementById("root")!).render(<App />);
