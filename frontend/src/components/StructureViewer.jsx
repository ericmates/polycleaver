import { useEffect, useMemo, useRef } from 'react'
import { Canvas, useFrame, useThree } from '@react-three/fiber'
import { OrbitControls, Environment } from '@react-three/drei'
import { EffectComposer, SSAO, Bloom, ToneMapping } from '@react-three/postprocessing'
import * as THREE from 'three'
import { ConvexGeometry } from 'three/addons/geometries/ConvexGeometry.js'
import { useStore } from '../store'
import { getColor, getRadius } from '../lib/elementData'
import WulffViewer from './WulffViewer'

// ---------------------------------------------------------------------------
// Shared: half-coloured bond cylinder (used for water AND scaffold/molecule bonds)
// ---------------------------------------------------------------------------

function ColoredBond({ posA, posB, colorA, colorB, radius = 0.09 }) {
  const vA   = useMemo(() => new THREE.Vector3(...posA), [posA])
  const vB   = useMemo(() => new THREE.Vector3(...posB), [posB])
  const len  = useMemo(() => vA.distanceTo(vB), [vA, vB])
  const dir  = useMemo(() => vB.clone().sub(vA).normalize(), [vA, vB])
  const quat = useMemo(
    () => new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir),
    [dir],
  )
  const midA = useMemo(() => vA.clone().lerp(vB, 0.25), [vA, vB])
  const midB = useMemo(() => vA.clone().lerp(vB, 0.75), [vA, vB])
  const half = len / 2

  return (
    <>
      <mesh position={midA} quaternion={quat}>
        <cylinderGeometry args={[radius, radius, half, 10]} />
        <meshPhysicalMaterial color={colorA} roughness={0.35} metalness={0.05} />
      </mesh>
      <mesh position={midB} quaternion={quat}>
        <cylinderGeometry args={[radius, radius, half, 10]} />
        <meshPhysicalMaterial color={colorB} roughness={0.35} metalness={0.05} />
      </mesh>
    </>
  )
}

// ---------------------------------------------------------------------------
// Atom cloud — one InstancedMesh per element, excludes water atoms
// ---------------------------------------------------------------------------

function ElementInstances({ element, atoms }) {
  const ref    = useRef()
  const color  = useMemo(() => new THREE.Color(getColor(element)), [element])
  const radius = getRadius(element)

  useEffect(() => {
    if (!ref.current) return
    const m = new THREE.Matrix4()
    atoms.forEach((a, i) => {
      m.setPosition(a.x, a.y, a.z)
      ref.current.setMatrixAt(i, m)
    })
    ref.current.instanceMatrix.needsUpdate = true
  }, [atoms])

  return (
    <instancedMesh ref={ref} args={[null, null, atoms.length]} castShadow receiveShadow>
      <icosahedronGeometry args={[radius, 3]} />
      <meshPhysicalMaterial
        color={color}
        roughness={0.30}
        metalness={0.05}
        clearcoat={0.7}
        clearcoatRoughness={0.15}
      />
    </instancedMesh>
  )
}

function AtomCloud({ atoms, waterSet }) {
  const byElement = useMemo(() => {
    const map = {}
    atoms.forEach((a, i) => {
      if (waterSet.has(i)) return
      if (!map[a.element]) map[a.element] = []
      map[a.element].push(a)
    })
    return map
  }, [atoms, waterSet])

  return Object.entries(byElement).map(([el, group]) => (
    <ElementInstances key={el} element={el} atoms={group} />
  ))
}

// ---------------------------------------------------------------------------
// Water sticks — small O sphere + half-coloured O-H cylinders, no H sphere
// ---------------------------------------------------------------------------

function WaterSticks({ atoms, waterBonds }) {
  const oxygens = useMemo(() => {
    if (!waterBonds?.length) return []
    const seen = new Set()
    const list = []
    waterBonds.forEach(b => {
      const idx = b.elem_a === 'O' ? b.a : b.b
      if (!seen.has(idx)) { seen.add(idx); list.push(atoms[idx]) }
    })
    return list
  }, [atoms, waterBonds])

  if (!waterBonds?.length) return null

  return (
    <>
      {/* Junction sphere at each oxygen */}
      {oxygens.map((a, i) => (
        <mesh key={i} position={[a.x, a.y, a.z]}>
          <icosahedronGeometry args={[0.18, 2]} />
          <meshPhysicalMaterial color={getColor('O')} roughness={0.30}
            metalness={0.05} clearcoat={0.6} />
        </mesh>
      ))}

      {/* Half-coloured O-H cylinders — use PBC-corrected positions when present */}
      {waterBonds.map((b, i) => (
        <ColoredBond
          key={i}
          posA={b.pos_a ?? [atoms[b.a].x, atoms[b.a].y, atoms[b.a].z]}
          posB={b.pos_b ?? [atoms[b.b].x, atoms[b.b].y, atoms[b.b].z]}
          colorA={getColor(b.elem_a)}
          colorB={getColor(b.elem_b)}
          radius={0.08}
        />
      ))}
    </>
  )
}

// ---------------------------------------------------------------------------
// Chemical bonds (scaffold + non-water molecules) — toggle only
// ---------------------------------------------------------------------------

function ChemBonds({ bonds, atoms }) {
  if (!bonds?.length) return null
  return bonds.map((b, i) => (
    <ColoredBond
      key={i}
      posA={[atoms[b.a].x, atoms[b.a].y, atoms[b.a].z]}
      posB={[atoms[b.b].x, atoms[b.b].y, atoms[b.b].z]}
      colorA={getColor(b.elem_a)}
      colorB={getColor(b.elem_b)}
      radius={0.08}
    />
  ))
}

// ---------------------------------------------------------------------------
// Coordination polyhedra — convex hull of O vertices, PBC-corrected
// ---------------------------------------------------------------------------

function SinglePolyhedron({ tet }) {
  const geometry = useMemo(() => {
    const pts = tet.vertices.map(v => new THREE.Vector3(...v))
    try { return new ConvexGeometry(pts) } catch { return null }
  }, [tet])

  if (!geometry) return null

  return (
    <mesh geometry={geometry}>
      <meshPhysicalMaterial
        color={tet.color}
        transparent
        opacity={0.22}
        side={THREE.DoubleSide}
        depthWrite={false}
        roughness={0.2}
      />
    </mesh>
  )
}

function CoordinationPolyhedra({ polyhedra }) {
  if (!polyhedra?.length) return null
  return polyhedra.map((tet, i) => <SinglePolyhedron key={i} tet={tet} />)
}

// ---------------------------------------------------------------------------
// Unit cell parallelepiped
// ---------------------------------------------------------------------------

const CELL_EDGES = [
  [0,1],[0,2],[0,4],[1,3],[1,5],[2,3],[2,6],[4,5],[4,6],[3,7],[5,7],[6,7],
]

function UnitCell({ matrix }) {
  const geometry = useMemo(() => {
    const [a, b, c] = matrix
    const corners = []
    for (let h = 0; h <= 1; h++)
      for (let k = 0; k <= 1; k++)
        for (let l = 0; l <= 1; l++)
          corners.push([h*a[0]+k*b[0]+l*c[0], h*a[1]+k*b[1]+l*c[1], h*a[2]+k*b[2]+l*c[2]])
    const pts = []
    CELL_EDGES.forEach(([i, j]) => { pts.push(...corners[i], ...corners[j]) })
    const geo = new THREE.BufferGeometry()
    geo.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3))
    return geo
  }, [matrix])

  return (
    <lineSegments geometry={geometry}>
      <lineBasicMaterial color="#3b82f6" transparent opacity={0.7} />
    </lineSegments>
  )
}

// ---------------------------------------------------------------------------
// Camera — auto zoom-to-fit; slab view from +X axis, bulk from default angle
// ---------------------------------------------------------------------------

function CameraRig({ structureData }) {
  const { camera } = useThree()
  const controlsRef = useRef()
  const target  = useRef({ pos: new THREE.Vector3(), look: new THREE.Vector3() })
  const lerping = useRef(false)
  const isSlab  = !!structureData?.meta

  useEffect(() => {
    if (!structureData?.atoms?.length) return
    const atoms = structureData.atoms
    const cx = atoms.reduce((s, a) => s + a.x, 0) / atoms.length
    const cy = atoms.reduce((s, a) => s + a.y, 0) / atoms.length
    const cz = atoms.reduce((s, a) => s + a.z, 0) / atoms.length
    const rSq = atoms.reduce((m, a) => {
      const d = (a.x-cx)**2 + (a.y-cy)**2 + (a.z-cz)**2
      return d > m ? d : m
    }, 0)
    const r = Math.sqrt(rSq) || 10
    target.current.look.set(cx, cy, cz)
    if (isSlab) {
      // Side view from +X: shows slab layers (Z is surface normal)
      target.current.pos.set(cx + r * 3.0, cy, cz)
    } else {
      target.current.pos.set(cx, cy + r * 0.3, cz + r * 2.5)
    }
    lerping.current = true
  }, [structureData])

  useFrame(() => {
    if (!lerping.current) return
    camera.position.lerp(target.current.pos, 0.08)
    if (controlsRef.current) {
      controlsRef.current.target.lerp(target.current.look, 0.08)
      controlsRef.current.update()
    }
    if (camera.position.distanceTo(target.current.pos) < 0.05) lerping.current = false
  })

  return <OrbitControls ref={controlsRef} makeDefault enableDamping dampingFactor={0.08} />
}

// ---------------------------------------------------------------------------
// Main viewer
// ---------------------------------------------------------------------------

export default function StructureViewer() {
  const { bulkJson, slabs, viewerTarget, showBonds, showUnitCell, showEffects,
          viewerMode, wulffData } = useStore()

  const structureData = viewerTarget === null
    ? bulkJson
    : (slabs[viewerTarget] ?? bulkJson)

  const waterSet = useMemo(
    () => new Set(structureData?.water_atom_indices ?? []),
    [structureData],
  )

  const showWulff = viewerMode === 'wulff' && wulffData

  if (!structureData && !showWulff) {
    return (
      <div className="w-full h-full flex items-center justify-center text-muted select-none">
        <div className="text-center space-y-2">
          <svg className="w-16 h-16 mx-auto opacity-20" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1}
              d="M20 7l-8-4-8 4m16 0l-8 4m8-4v10l-8 4m0-10L4 7m8 10V7" />
          </svg>
          <p className="text-sm">Upload a CIF file to begin</p>
        </div>
      </div>
    )
  }

  return (
    <Canvas
      dpr={[1, 2]}
      camera={{ position: [0, 0, 50], fov: 45 }}
      shadows
      gl={{ antialias: true, toneMapping: THREE.ACESFilmicToneMapping }}
    >
      <color attach="background" args={['#080d1a']} />

      {showWulff ? (
        <WulffViewer wulffData={wulffData} />
      ) : (
        <>
          <CameraRig structureData={structureData} />
          <Environment preset="studio" />
          <directionalLight position={[15, 25, 20]} intensity={1.5} castShadow />
          <ambientLight intensity={0.3} />

          <AtomCloud atoms={structureData.atoms} waterSet={waterSet} />
          <WaterSticks atoms={structureData.atoms} waterBonds={structureData.water_bonds} />
          <CoordinationPolyhedra polyhedra={structureData.polyhedra} />

          {showUnitCell && <UnitCell matrix={structureData.lattice.matrix} />}

          {/* Bonds: scaffold + molecule bonds (chemical only, not all-atom contacts) */}
          {showBonds && <ChemBonds bonds={structureData.bonds} atoms={structureData.atoms} />}

          {showEffects && (
            <EffectComposer>
              <SSAO radius={0.04} intensity={20} luminanceInfluence={0.6} />
              <Bloom luminanceThreshold={0.6} intensity={0.25} mipmapBlur />
              <ToneMapping />
            </EffectComposer>
          )}
        </>
      )}
    </Canvas>
  )
}
