import { useMemo, useEffect, useRef } from 'react'
import { useFrame, useThree } from '@react-three/fiber'
import { Html, OrbitControls } from '@react-three/drei'
import * as THREE from 'three'

// ---------------------------------------------------------------------------
// Single Wulff face — triangle fan, colored by energy
// ---------------------------------------------------------------------------

function WulffFace({ vertices, energy, minE, maxE }) {
  const geometry = useMemo(() => {
    if (vertices.length < 3) return null
    const pts = vertices.map(v => new THREE.Vector3(...v))
    const center = new THREE.Vector3()
    pts.forEach(p => center.add(p))
    center.divideScalar(pts.length)

    const positions = []
    for (let i = 0; i < pts.length; i++) {
      positions.push(
        ...center.toArray(),
        ...pts[i].toArray(),
        ...pts[(i + 1) % pts.length].toArray(),
      )
    }
    const geo = new THREE.BufferGeometry()
    geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3))
    geo.computeVertexNormals()
    return geo
  }, [vertices])

  const color = useMemo(() => {
    const t = maxE > minE ? (energy - minE) / (maxE - minE) : 0
    const hue = Math.round((1 - t) * 220)   // blue (220°) → red (0°)
    return `hsl(${hue}, 75%, 55%)`
  }, [energy, minE, maxE])

  if (!geometry) return null
  return (
    <mesh geometry={geometry}>
      <meshPhysicalMaterial
        color={color}
        transparent
        opacity={0.65}
        side={THREE.DoubleSide}
        depthWrite={false}
        roughness={0.25}
        metalness={0.05}
      />
    </mesh>
  )
}

// ---------------------------------------------------------------------------
// Wireframe edges for all faces
// ---------------------------------------------------------------------------

function WulffEdges({ faces }) {
  const geometry = useMemo(() => {
    const pts = []
    faces.forEach(face => {
      const v = face.vertices
      for (let i = 0; i < v.length; i++) {
        pts.push(...v[i], ...v[(i + 1) % v.length])
      }
    })
    const geo = new THREE.BufferGeometry()
    geo.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3))
    return geo
  }, [faces])

  return (
    <lineSegments geometry={geometry}>
      <lineBasicMaterial color="#94a3b8" transparent opacity={0.5} />
    </lineSegments>
  )
}

// ---------------------------------------------------------------------------
// Camera auto-fit for Wulff shape
// ---------------------------------------------------------------------------

function WulffCameraRig({ faces }) {
  const { camera } = useThree()
  const controlsRef = useRef()
  const target = useRef({ pos: new THREE.Vector3(), look: new THREE.Vector3() })
  const lerping = useRef(false)

  useEffect(() => {
    if (!faces.length) return
    const allVerts = faces.flatMap(f => f.vertices.map(v => new THREE.Vector3(...v)))
    const cx = allVerts.reduce((s, v) => s + v.x, 0) / allVerts.length
    const cy = allVerts.reduce((s, v) => s + v.y, 0) / allVerts.length
    const cz = allVerts.reduce((s, v) => s + v.z, 0) / allVerts.length
    const r = allVerts.reduce((m, v) => {
      const d = Math.sqrt((v.x - cx) ** 2 + (v.y - cy) ** 2 + (v.z - cz) ** 2)
      return d > m ? d : m
    }, 0) || 5
    target.current.look.set(cx, cy, cz)
    target.current.pos.set(cx, cy + r * 0.3, cz + r * 2.8)
    lerping.current = true
  }, [faces])

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
// Labels: one per unique (hkl), placed at the outermost face centroid
// ---------------------------------------------------------------------------

function WulffLabels({ faces, minE, maxE }) {
  const labels = useMemo(() => {
    const byHkl = {}
    faces.forEach(face => {
      const key = face.hkl.join(',')
      const verts = face.vertices
      const cx = verts.reduce((s, v) => s + v[0], 0) / verts.length
      const cy = verts.reduce((s, v) => s + v[1], 0) / verts.length
      const cz = verts.reduce((s, v) => s + v[2], 0) / verts.length
      const r = Math.sqrt(cx * cx + cy * cy + cz * cz)
      if (!byHkl[key] || r > byHkl[key].r) {
        byHkl[key] = { hkl: face.hkl, energy: face.energy, center: [cx, cy, cz], r }
      }
    })
    return Object.values(byHkl)
  }, [faces])

  return (
    <>
      {labels.map((lbl, i) => {
        const t = maxE > minE ? (lbl.energy - minE) / (maxE - minE) : 0
        const hue = Math.round((1 - t) * 220)
        const borderColor = `hsl(${hue}, 60%, 45%)`
        return (
          <Html key={i} position={lbl.center} center distanceFactor={10}>
            <div style={{
              background: 'rgba(8,13,26,0.88)',
              border: `1px solid ${borderColor}`,
              borderRadius: 6,
              padding: '3px 7px',
              pointerEvents: 'none',
              whiteSpace: 'nowrap',
              userSelect: 'none',
            }}>
              <div style={{ color: `hsl(${hue}, 80%, 65%)`, fontFamily: 'monospace', fontSize: 11, fontWeight: 600 }}>
                ({lbl.hkl.join(',')})
              </div>
              <div style={{ color: '#94a3b8', fontFamily: 'monospace', fontSize: 10 }}>
                {lbl.energy.toFixed(3)} J/m²
              </div>
            </div>
          </Html>
        )
      })}
    </>
  )
}

// ---------------------------------------------------------------------------
// Main export — renders inside an existing Canvas
// ---------------------------------------------------------------------------

export default function WulffViewer({ wulffData }) {
  const { faces } = wulffData
  const energies = faces.map(f => f.energy)
  const minE = Math.min(...energies)
  const maxE = Math.max(...energies)

  return (
    <>
      <WulffCameraRig faces={faces} />
      <ambientLight intensity={0.6} />
      <directionalLight position={[10, 15, 10]} intensity={1.0} />

      {faces.map((face, i) => (
        <WulffFace key={i} vertices={face.vertices} energy={face.energy} minE={minE} maxE={maxE} />
      ))}
      <WulffEdges faces={faces} />
      <WulffLabels faces={faces} minE={minE} maxE={maxE} />
    </>
  )
}
