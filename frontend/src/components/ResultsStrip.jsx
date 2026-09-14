import { useStore } from '../store'
import SlabCard from './SlabCard'

export default function ResultsStrip() {
  const { slabs, jobId, jobStatus } = useStore()

  if (slabs.length === 0 && jobStatus !== 'running') return null

  return (
    <div className="border-t border-border bg-panel/80 backdrop-blur-sm">
      <div className="flex items-center gap-2 px-4 pt-3 pb-1">
        <span className="text-xs font-semibold uppercase tracking-widest text-muted">
          Surfaces
        </span>
        {slabs.length > 0 && (
          <span className="text-xs font-mono text-accent">{slabs.length}</span>
        )}
      </div>
      <div className="flex gap-3 overflow-x-auto px-4 pb-4 pt-1 scrollbar-thin scrollbar-thumb-border">
        {slabs.map((slab, i) => (
          <SlabCard key={i} slab={slab} index={i} jobId={jobId} />
        ))}
        {jobStatus === 'running' && slabs.length === 0 && (
          <div className="flex-shrink-0 w-48 h-28 rounded-xl border border-dashed border-border
                          flex items-center justify-center text-muted text-xs animate-pulse">
            Searching…
          </div>
        )}
      </div>
    </div>
  )
}
