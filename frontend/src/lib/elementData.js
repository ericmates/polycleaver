export const CPK = {
  H:  '#FFFFFF', C:  '#909090', N:  '#3050F8', O:  '#FF0D0D',
  F:  '#90E050', P:  '#FF8000', S:  '#FFFF30', Cl: '#1FF01F',
  Ca: '#3DFF00', Na: '#AB5CF2', K:  '#8F40D4', Mg: '#8AFF00',
  Al: '#BFA6A6', Si: '#F0C8A0', Fe: '#E06633', Mn: '#9C7AC7',
  Zn: '#7D80B0', Cu: '#C88033', B:  '#FFB5B5', Br: '#A62929',
  I:  '#940094', As: '#BD80E3', Se: '#FFA100',
}

export const VDW = {
  H:  1.20, C:  1.70, N:  1.55, O:  1.52, F:  1.47,
  P:  1.80, S:  1.80, Cl: 1.75, Ca: 2.31, Na: 2.27,
  K:  2.75, Mg: 1.73, Al: 1.84, Si: 2.10, Fe: 2.04,
  Mn: 2.05, Zn: 2.10, Cu: 1.96, B:  1.92, Br: 1.85,
  I:  1.98, As: 1.85, Se: 1.90,
}

// Per-element display radius overrides (smaller for light atoms)
const DISPLAY_RADIUS = {
  H: 0.14,
  O: 0.34,
  N: 0.36,
}

export const getColor  = sym => CPK[sym] ?? '#FF69B4'
export const getRadius = sym => DISPLAY_RADIUS[sym] ?? (VDW[sym] ?? 2.0) * 0.4
