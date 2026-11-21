import { useEffect, useRef, useState } from 'react'
import { useStore } from '../store'
import { startGenerate, pollStatus } from '../lib/api'

function parseMillerInput(raw) {
  const t = raw.trim()
  if (!t) return null
  if (/^\d+$/.test(t)) return { mode: 'max', value: parseInt(t) }

  // Single tuple: "2,1,0" or "(2,1,0)" or "[2,1,0]"
  const clean = t.replace(/[()[\]\s]/g, '')
  const nums = clean.split(',').map(Number)
  if (nums.length === 3 && nums.every(Number.isFinite))
    return { mode: 'single', value: nums }

  // List: semicolon-separated triples
  const triples = t.split(';').map(s =>
    s.trim().replace(/[()[\]\s]/g, '').split(',').map(Number)
  )
  if (triples.every(tr => tr.length === 3 && tr.every(Number.isFinite)))
    return { mode: 'list', value: triples }

  return null
}

export default function ControlPanel() {
  const { sessionId, jobStatus, setJobStarted, updateStatus,
          showBonds, showUnitCell, showEffects,
          toggleBonds, toggleCell, toggleEffects,
          progress, completedIndices, totalIndices, currentHkl, slabs, jobError } = useStore()

  const [thickness, setThickness] = useState(15)
  const [vacuum,    setVacuum]    = useState(15)
  const [dipole,    setDipole]    = useState(0.3)
  const [miller,    setMiller]    = useState('1')
  const [millerErr, setMillerErr] = useState(null)

  const pollRef = useRef(null)

  async function generate() {
    const mi = parseMillerInput(miller)
    if (!mi) { setMillerErr('Invalid format. Use: 1 · (2,1,0) · 1,0,0;0,1,0'); return }
    setMillerErr(null)

    const { job_id } = await startGenerate(sessionId, {
      thickness, vacuum, dipole_tol: dipole, miller_input: mi,
    })
    setJobStarted(job_id)

    pollRef.current = setInterval(async () => {
      const data = await pollStatus(job_id)
      updateStatus(data)
      if (data.status === 'done' || data.status === 'error') {
        clearInterval(pollRef.current)
      }
    }, 500)
  }

  useEffect(() => () => clearInterval(pollRef.current), [])

  const running    = jobStatus === 'running'
  const lastProg   = progress[progress.length - 1]
  const overallPct = totalIndices > 0 ? Math.round((completedIndices / totalIndices) * 100) : 0

  return (
    <div className="flex flex-col gap-5">
      {/* Parameters */}
      <div>
        <p className="text-xs font-semibold uppercase tracking-widest text-muted mb-3">Parameters</p>
        <div className="flex flex-col gap-2">
          <NumInput label="Thickness (Å)" value={thickness} onChange={setThickness} min={5}  step={0.5}  />
          <NumInput label="Vacuum (Å)"    value={vacuum}    onChange={setVacuum}    min={5}  step={0.5}  />
          <NumInput label="Dipole tol."   value={dipole}    onChange={setDipole}    min={0}  step={0.05} />

          <div className="flex flex-col gap-1">
            <label className="text-xs text-slate-400">Miller indices</label>
            <input
              value={miller}
              onChange={e => setMiller(e.target.value)}
              placeholder="1 · (2,1,0) · 1,0,0;0,1,0"
              className="bg-bg border border-border text-xs text-white font-mono rounded-lg px-2 py-1.5 w-full
                         focus:outline-none focus:border-accent placeholder:text-muted"
            />
            {millerErr && <p className="text-xs text-red-400">{millerErr}</p>}
            <p className="text-[10px] text-muted">int = max index · tuple = single · list = multiple (;)</p>
          </div>
        </div>
      </div>

      {/* Generate button */}
      <button
        onClick={generate}
        disabled={!sessionId || running}
        className={`
          w-full py-2.5 rounded-xl text-sm font-semibold transition-all duration-200
          ${sessionId && !running
            ? 'bg-accent hover:bg-accent-dim text-white shadow-lg shadow-accent/20'
            : 'bg-border text-muted cursor-not-allowed'}
        `}
      >
        {running ? 'Generating…' : 'Generate Surfaces'}
      </button>

      {/* Progress */}
      {(running || (jobStatus === 'done' && progress.length > 0)) && (
        <div className="space-y-2">
          {currentHkl && running && (
            <p className="text-xs font-mono text-slate-400">
              Processing ({currentHkl.join(',')})
              {lastProg && ` — Cut ${lastProg.cut}${lastProg.total_cuts ? `/${lastProg.total_cuts}` : ''} (${lastProg.percent.toFixed(0)}%)`}
            </p>
          )}
          {lastProg && <Bar pct={lastProg.percent} />}
          {totalIndices > 1 && (
            <>
              <p className="text-[10px] text-muted">{completedIndices}/{totalIndices} indices</p>
              <Bar pct={overallPct} dim />
            </>
          )}
        </div>
      )}

      {/* Results summary */}
      {jobStatus === 'done' && (
        <div className="text-xs">
          {slabs.length > 0
            ? <span className="text-green-400">{slabs.length} surface{slabs.length > 1 ? 's' : ''} found</span>
            : <span className="text-yellow-400">No surfaces passed the dipole filter. Try raising the tolerance.</span>}
        </div>
      )}
      {jobStatus === 'error' && (
        <p className="text-xs text-red-400 font-mono break-all">{jobError}</p>
      )}

      {/* Viewer toggles */}
      <div>
        <p className="text-xs font-semibold uppercase tracking-widest text-muted mb-2">Viewer</p>
        <div className="flex flex-col gap-1.5">
          <Toggle label="Unit cell" value={showUnitCell} onToggle={toggleCell}   />
          <Toggle label="Bonds"     value={showBonds}    onToggle={toggleBonds}  />
          <Toggle label="Effects"   value={showEffects}  onToggle={toggleEffects}/>
        </div>
      </div>
    </div>
  )
}

function NumInput({ label, value, onChange, min, step }) {
  return (
    <div className="flex items-center justify-between">
      <label className="text-xs text-slate-400">{label}</label>
      <input
        type="number"
        value={value}
        min={min}
        step={step}
        onChange={e => onChange(parseFloat(e.target.value))}
        className="w-24 bg-bg border border-border text-xs text-white font-mono text-right
                   rounded-lg px-2 py-1.5 focus:outline-none focus:border-accent"
      />
    </div>
  )
}

function Bar({ pct, dim }) {
  return (
    <div className={`w-full h-1 rounded-full ${dim ? 'bg-border' : 'bg-white/10'}`}>
      <div
        className={`h-1 rounded-full transition-all duration-300 ${dim ? 'bg-accent/40' : 'bg-accent'}`}
        style={{ width: `${pct}%` }}
      />
    </div>
  )
}

function Toggle({ label, value, onToggle }) {
  return (
    <button
      onClick={onToggle}
      className="flex items-center justify-between w-full px-3 py-1.5 rounded-lg
                 bg-white/[0.03] hover:bg-white/[0.06] transition-colors"
    >
      <span className="text-xs text-slate-400">{label}</span>
      <div className={`w-8 h-4 rounded-full transition-colors ${value ? 'bg-accent' : 'bg-border'}`}>
        <div className={`w-3 h-3 rounded-full bg-white m-0.5 transition-transform
                         ${value ? 'translate-x-4' : 'translate-x-0'}`} />
      </div>
    </button>
  )
}
