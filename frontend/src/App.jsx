import { useStore } from './store'
import UploadZone from './components/UploadZone'
import ControlPanel from './components/ControlPanel'
import OptimizePanel from './components/OptimizePanel'
import StructureViewer from './components/StructureViewer'
import ResultsStrip from './components/ResultsStrip'
import SlabInfoPanel from './components/SlabInfoPanel'

function ViewerOverlay() {
  const { bulkJson, slabs, viewerTarget, filename, viewerMode, wulffData } = useStore()

  if (viewerMode === 'wulff' && wulffData) {
    return (
      <div className="absolute inset-0 pointer-events-none p-4">
        <div className="bg-bg/70 backdrop-blur-sm border border-border rounded-lg px-3 py-1.5 inline-block">
          <p className="text-xs font-mono text-slate-300">Wulff construction</p>
        </div>
        <div className="absolute top-4 right-4 bg-bg/70 backdrop-blur-sm border border-border
                        rounded-lg px-3 py-1.5 text-right">
          <p className="text-xs font-mono text-slate-400">
            γ̄ = {wulffData.weighted_surface_energy.toFixed(3)} J/m²
          </p>
          <p className="text-xs font-mono text-slate-500">
            anisotropy {wulffData.anisotropy.toFixed(3)}
          </p>
        </div>
      </div>
    )
  }

  const data   = viewerTarget === null ? bulkJson : (slabs[viewerTarget] ?? bulkJson)
  const isSlab = viewerTarget !== null
  if (!data) return null

  const label = isSlab
    ? `Slab (${data.meta.miller_index.join(',')}) — ${data.meta.formula}`
    : `Bulk — ${filename ?? ''}`

  return (
    <div className="absolute inset-0 pointer-events-none p-4">
      {/* Top-left: current view label */}
      <div className="bg-bg/70 backdrop-blur-sm border border-border rounded-lg px-3 py-1.5 inline-block">
        <p className="text-xs font-mono text-slate-300">{label}</p>
      </div>

      {/* Top-right: atom count (only when no slab panel showing) */}
      {!isSlab && data.atoms && (
        <div className="absolute top-4 right-4 bg-bg/70 backdrop-blur-sm border border-border
                        rounded-lg px-3 py-1.5 text-right">
          <p className="text-xs font-mono text-slate-400">{data.atoms.length} atoms</p>
        </div>
      )}
    </div>
  )
}

export default function App() {
  return (
    <div className="flex h-screen bg-bg text-white overflow-hidden font-sans">
      {/* Sidebar */}
      <aside className="w-72 flex-shrink-0 flex flex-col gap-6 p-5
                        bg-panel border-r border-border overflow-y-auto">
        <div className="flex items-center gap-2 pt-1">
          <HexIcon />
          <span className="font-semibold text-base tracking-tight">Surface Lab</span>
        </div>
        <UploadZone />
        <div className="border-t border-border" />
        <ControlPanel />
        <div className="border-t border-border" />
        <OptimizePanel />
      </aside>

      {/* Main area */}
      <main className="flex-1 flex flex-col min-w-0">
        <div className="flex-1 relative min-h-0">
          <StructureViewer />
          <ViewerOverlay />
          <SlabInfoPanel />
        </div>
        <ResultsStrip />
      </main>
    </div>
  )
}

function HexIcon() {
  return (
    <svg className="w-6 h-6 text-accent" fill="currentColor" viewBox="0 0 24 24">
      <path d="M12 2l9 5v10l-9 5-9-5V7z" />
    </svg>
  )
}
