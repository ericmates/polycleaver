import { useEffect, useRef } from 'react'
import { useStore } from '../store'
import { startOptimize, pollOptStatus } from '../lib/api'

export default function OptimizePanel() {
  const {
    sessionId, jobId, jobStatus, slabs,
    optStatus, optPhaseLabel, optCompleted, optTotal,
    slabResults, wulffData, viewerMode,
    setOptStarted, updateOptStatus, setViewerMode,
  } = useStore()

  const pollRef = useRef(null)

  async function handleOptimize() {
    const { opt_job_id } = await startOptimize(sessionId, jobId)
    setOptStarted(opt_job_id)

    pollRef.current = setInterval(async () => {
      const data = await pollOptStatus(opt_job_id)
      updateOptStatus(data)
      if (data.status === 'done' || data.status === 'error') {
        clearInterval(pollRef.current)
      }
    }, 1000)
  }

  useEffect(() => () => clearInterval(pollRef.current), [])

  // Only show when generation is done and there are slabs
  if (jobStatus !== 'done' || slabs.length === 0) return null

  const isIdle    = optStatus === 'idle'
  const isRunning = optStatus === 'running'
  const isDone    = optStatus === 'done'
  const isError   = optStatus === 'error'

  const pct = optTotal > 0 ? Math.round((optCompleted / optTotal) * 100) : 0

  return (
    <div className="flex flex-col gap-3">
      <p className="text-xs font-semibold uppercase tracking-widest text-muted">
        Surface Energies
      </p>

      {isIdle && (
        <button
          onClick={handleOptimize}
          className="w-full text-xs bg-accent/15 hover:bg-accent/25 border border-accent/30
                     hover:border-accent/60 rounded-lg py-2 text-accent transition-colors"
        >
          Optimize with MACE-MP
        </button>
      )}

      {isRunning && (
        <div className="space-y-1.5">
          <p className="text-xs text-slate-400">{optPhaseLabel ?? 'Initializing…'}</p>
          <div className="w-full h-1 bg-white/10 rounded-full overflow-hidden">
            <div
              className="h-full bg-accent rounded-full transition-all duration-500"
              style={{ width: `${pct}%` }}
            />
          </div>
          {optTotal > 0 && (
            <p className="text-[10px] text-muted">{optCompleted}/{optTotal} slabs</p>
          )}
        </div>
      )}

      {(isDone || (isRunning && slabResults.length > 0)) && (
        <div className="space-y-0.5">
          {slabResults.map((r, i) => (
            <div key={i} className="flex justify-between items-baseline text-xs py-0.5
                                    border-b border-white/5">
              <span className="font-mono text-accent text-[11px]">
                ({r.hkl.join(',')})
              </span>
              <span className="font-mono text-slate-300">
                {r.e_surf_j_m2.toFixed(3)} J/m²
              </span>
            </div>
          ))}
        </div>
      )}

      {isDone && wulffData && (
        <button
          onClick={() => setViewerMode(viewerMode === 'wulff' ? 'structure' : 'wulff')}
          className="w-full text-xs border border-border hover:border-accent/50
                     rounded-lg py-1.5 text-muted hover:text-white transition-colors"
        >
          {viewerMode === 'wulff' ? '← Structure view' : 'View Wulff shape →'}
        </button>
      )}

      {isError && (
        <p className="text-xs text-red-400 font-mono break-all">
          Optimization failed. Check server logs.
        </p>
      )}
    </div>
  )
}
