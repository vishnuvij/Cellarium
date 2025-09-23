import React, { useEffect, useState } from "react";
import ReactDOM from "react-dom/client";
import Plot from "react-plotly.js";
import axios from "axios";

const App = () => {
  const [points, setPoints] = useState<number[][]>([]);
  const [cellIds, setCellIds] = useState<string[]>([]);
  const [clusters, setClusters] = useState<string[]>([]);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    const fetchData = async () => {
      try {
        const res = await axios.get("http://localhost:8000/embedding/scanpy-pbmc3k");
        setPoints(res.data.umap);
        setCellIds(res.data.cell_ids);
        setClusters(res.data.clusters);
      } catch (err) {
        console.error("❌ Failed to fetch embedding:", err);
      } finally {
        setLoading(false);
      }
    };
    fetchData();
  }, []);

  if (loading) return <div style={{ padding: "2rem" }}>Loading 3D scatter...</div>;

  return (
    <div style={{ width: "100%", height: "100vh" }}>
      <h2 style={{ textAlign: "center" }}>Cellarium 3D Embedding</h2>
      <Plot
        data={[
          {
            x: points.map(p => p[0]),
            y: points.map(p => p[1]),
            z: points.map(p => p[2]),
            mode: "markers",
            type: "scatter3d",
            text: cellIds.map((id, i) => `Cell: ${id}<br>Cluster: ${clusters[i]}`), // 👈 tooltip
            hoverinfo: "text",
            marker: {
              size: 3,
              color: clusters.map(c => parseInt(c)), // 👈 color by cluster
              colorscale: "Viridis"
            }
          }
        ]}
        layout={{
          autosize: true,
          height: 700,
          scene: {
            xaxis: { title: "UMAP-1" },
            yaxis: { title: "UMAP-2" },
            zaxis: { title: "UMAP-3" }
          }
        }}
        style={{ width: "100%", height: "100%" }}
      />
    </div>
  );
};

ReactDOM.createRoot(document.getElementById("root")!).render(<App />);
