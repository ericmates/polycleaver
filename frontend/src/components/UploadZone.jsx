import { useRef, useState } from 'react'
import { useStore } from '../store'
import { uploadCif } from '../lib/api'

export default function UploadZone() {
  const fileRef = useRef()
  const [dragging, setDragging] = useState(false)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState(null)
  const { setUploadResult, analysisJson, filename } = useStore()

  async function handle(file) {
    if (!file) return
    setLoading(true)
    setError(null)
    try {
      const data = await uploadCif(file)
      setUploadResult(data)
    } catch (e) {
      setError(e.message)
    } finally {
      setLoading(false)
    }
  }

  function onDrop(e) {
    e.preventDefault()
    setDragging(false)
    handle(e.dataTransfer.files[0])
  }

  return (
    <div>
      <p className="text-xs font-semibold uppercase tracking-widest text-muted mb-2">Structure</p>
      <div
        onDrop={onDrop}
        onDragOver={(e) => { e.preventDefault(); setDragging(true) }}
        onDragLeave={() => setDragging(false)}
        onClick={() => fileRef.current.click()}
        className={`
          relative border-2 border-dashed rounded-xl p-6 text-center cursor-pointer
          transition-all duration-200 select-none
          ${dragging ? 'border-accent bg-accent/10' : 'border-border hover:border-accent/60 hover:bg-white/[0.02]'}
          ${loading ? 'opacity-60 pointer-events-none' : ''}
        `}
      >
        <input
          ref={fileRef}
          type="file"
          accept=".cif"
          hidden
          onChange={(e) => handle(e.target.files[0])}
        />
        {loading ? (
          <div className="flex items-center justify-center gap-2 text-accent">
            <Spinner />
            <span className="text-sm">Parsing…</span>
          </div>
        ) : (
          <>
            <CifIcon />
            <p className="text-sm text-slate-300 mt-2">
              {filename ? filename : 'Drop CIF or click'}
            </p>
            <p className="text-xs text-muted mt-1">Crystallography Open Database format</p>
          </>
        )}
      </div>

      {error && (
        <p className="mt-2 text-xs text-red-400 font-mono">{error}</p>
      )}

      {analysisJson && (
        <div className="mt-3 flex flex-wrap gap-2">
          <Chip label={`${analysisJson.n_centers} scaffold`} />
          {analysisJson.n_cations > 0 && <Chip label={`${analysisJson.n_cations} cation${analysisJson.n_cations > 1 ? 's' : ''}`} />}
          {Object.entries(analysisJson.molecules).map(([f, n]) => (
            <Chip key={f} label={`${n}× ${f}`} />
          ))}
        </div>
      )}
    </div>
  )
}

function Chip({ label }) {
  return (
    <span className="px-2 py-0.5 rounded-full bg-accent/10 text-accent text-xs font-mono border border-accent/20">
      {label}
    </span>
  )
}

function CifIcon() {
  return (
    <svg className="mx-auto w-8 h-8 text-muted" fill="none" viewBox="0 0 24 24" stroke="currentColor">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
        d="M9 13h6m-3-3v6m5 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
    </svg>
  )
}

function Spinner() {
  return (
    <svg className="animate-spin h-4 w-4" viewBox="0 0 24 24" fill="none">
      <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
      <path className="opacity-75" fill="currentColor"
        d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
    </svg>
  )
}
