import { useStore } from '../store'
import { downloadUrl } from '../lib/api'

export default function SlabCard({ slab, index, jobId }) {
  const { selectSlab, viewerTarget } = useStore()
  const isActive = viewerTarget === index
  const m = slab.meta

  return (
    <div
      onClick={() => selectSlab(index)}
      className={`
        flex-shrink-0 w-48 cursor-pointer rounded-xl p-3 transition-all duration-200
        backdrop-blur-sm border animate-fade-in select-none
        ${isActive
          ? 'border-accent bg-accent/10 shadow-lg shadow-accent/20'
          : 'border-white/10 bg-white/[0.04] hover:border-accent/40 hover:bg-white/[0.07]'}
      `}
    >
      <div className="flex items-start justify-between gap-1 mb-2">
        <span className="font-mono text-sm text-accent font-semibold leading-none">
          ({m.miller_index.join(',')})
        </span>
        <a
          href={downloadUrl(jobId, index)}
          download
          onClick={e => e.stopPropagation()}
          className="text-muted hover:text-white transition-colors"
          title="Download CIF"
        >
          <DownloadIcon />
        </a>
      </div>

      <p className="text-xs text-white font-medium truncate mb-2">{m.formula}</p>

      <div className="grid grid-cols-2 gap-x-2 gap-y-1 text-[10px] font-mono text-slate-400">
        <span>d={m.thickness?.toFixed(1)} Å</span>
        <span>μ={m.dipole?.toFixed(3)}</span>
        <span>{m.is_symmetric ? '✓ sym' : '✗ asym'}</span>
        <span>{m.num_atoms} atoms</span>
      </div>

      {Object.keys(m.molecules ?? {}).length > 0 && (
        <div className="mt-2 flex flex-wrap gap-1">
          {Object.entries(m.molecules).map(([f, n]) => (
            <span key={f} className="text-[9px] font-mono px-1.5 py-0.5 rounded-full
                                     bg-white/5 border border-white/10 text-slate-400">
              {n}× {f}
            </span>
          ))}
        </div>
      )}
    </div>
  )
}

function DownloadIcon() {
  return (
    <svg className="w-3.5 h-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
        d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
    </svg>
  )
}
