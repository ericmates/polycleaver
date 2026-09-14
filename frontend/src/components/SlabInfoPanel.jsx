import { useStore } from '../store'

function Row({ label, value, accent }) {
  return (
    <div className="flex items-baseline justify-between gap-3 py-1 border-b border-white/5">
      <span className="text-[10px] text-muted uppercase tracking-wider flex-shrink-0">{label}</span>
      <span className={`text-xs font-mono text-right ${accent ? 'text-accent' : 'text-slate-300'}`}>
        {value}
      </span>
    </div>
  )
}

function Section({ title, children }) {
  return (
    <div className="mb-3">
      <p className="text-[9px] font-semibold uppercase tracking-widest text-muted/60 mb-1">{title}</p>
      {children}
    </div>
  )
}

export default function SlabInfoPanel() {
  const { slabs, viewerTarget, viewerMode, showBulk } = useStore()

  if (viewerMode === 'wulff' || viewerTarget === null || !slabs[viewerTarget]) return null
  const m = slabs[viewerTarget].meta

  const hkl        = `(${m.miller_index.join(', ')})`
  const symLabel   = m.is_symmetric ? 'Symmetric' : 'Asymmetric'
  const molEntries = Object.entries(m.molecules ?? {})

  return (
    <div className="absolute right-4 top-4 bottom-4 w-52 pointer-events-auto flex flex-col gap-0
                    bg-bg/80 backdrop-blur-md border border-border rounded-xl p-3 overflow-y-auto
                    text-white shadow-xl shadow-black/40">

      {/* Header */}
      <div className="mb-3">
        <p className="font-mono text-base text-accent font-semibold">{hkl}</p>
        <p className="text-xs text-slate-300 mt-0.5 break-all leading-tight">{m.formula}</p>
      </div>

      <Section title="Geometry">
        <Row label="Thickness"    value={`${m.thickness?.toFixed(2)} Å`} />
        <Row label="Surface area" value={`${m.surface_area?.toFixed(1)} Å²`} />
        <Row label="Symmetry"     value={symLabel} accent={m.is_symmetric} />
      </Section>

      <Section title="Polarity">
        <Row label="Dipole (Z)"   value={`${m.dipole?.toFixed(4)} eÅ`}
          accent={m.dipole != null && m.dipole < 0.5} />
      </Section>

      <Section title="Composition">
        <Row label="Atoms"        value={m.num_atoms} />
        <Row label="Multiplicity" value={m.multiplicity != null ? `${m.multiplicity.toFixed(2)}× bulk` : '—'} />
        <Row label="Scaffold"     value={m.n_centers} />
        <Row label="Cations"      value={m.n_cations} />
        {molEntries.length > 0 && molEntries.map(([f, n]) => (
          <Row key={f} label={f} value={n} />
        ))}
      </Section>

      <div className="mt-auto pt-2">
        <button
          onClick={showBulk}
          className="w-full text-xs text-muted hover:text-white border border-border
                     hover:border-accent/40 rounded-lg py-1.5 transition-colors"
        >
          ← Back to bulk
        </button>
      </div>
    </div>
  )
}
