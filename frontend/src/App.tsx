import React, { useEffect, useState } from 'react'
import axios from 'axios'

type Dataset = { id: string; filename: string; size_bytes: number }

export default function App() {
  const [datasets, setDatasets] = useState<Dataset[]>([])
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    const load = async () => {
      try {
        setLoading(true)
        setError(null)
        const r = await axios.get('http://localhost:8000/datasets')
        setDatasets(r.data)
      } catch (e: any) {
        setError(e?.message ?? 'Failed to load')
      } finally {
        setLoading(false)
      }
    }
    load()
  }, [])

  return (
    <div style={{ fontFamily: 'system-ui, sans-serif', padding: 16 }}>
      <h1>Cellarium</h1>
      <p>Backend: <code>http://localhost:8000</code></p>

      {loading && <p>Loading datasets…</p>}
      {error && <p style={{ color: 'crimson' }}>{error}</p>}

      <h2>Datasets in /data</h2>
      {datasets.length === 0 ? (
        <p>No .h5ad files found. Put files into the <code>data/</code> folder and refresh.</p>
      ) : (
        <ul>
          {datasets.map(d => (
            <li key={d.id}>
              <strong>{d.id}</strong> — {d.filename} ({d.size_bytes} bytes)
            </li>
          ))}
        </ul>
      )}
    </div>
  )
}
