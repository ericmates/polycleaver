import { create } from 'zustand'

export const useStore = create((set) => ({
  // Upload / session
  sessionId: null,
  bulkJson: null,
  analysisJson: null,
  filename: null,

  // Slab generation job
  jobId: null,
  jobStatus: 'idle',     // 'idle' | 'running' | 'done' | 'error'
  jobError: null,
  progress: [],
  completedIndices: 0,
  totalIndices: 0,
  currentHkl: null,

  // Results
  slabs: [],

  // Optimization job
  optJobId: null,
  optStatus: 'idle',     // 'idle' | 'running' | 'done' | 'error'
  optError: null,
  optPhase: null,
  optPhaseLabel: null,
  optCompleted: 0,
  optTotal: 0,
  slabResults: [],
  wulffData: null,

  // Viewer
  viewerTarget: null,    // null = bulk, number = slab index
  viewerMode: 'structure',  // 'structure' | 'wulff'
  showBonds: false,
  showUnitCell: true,
  showEffects: true,

  // Actions
  setUploadResult: ({ session_id, bulk, analysis, filename }) =>
    set({
      sessionId: session_id, bulkJson: bulk, analysisJson: analysis, filename,
      slabs: [], jobId: null, jobStatus: 'idle', jobError: null,
      progress: [], viewerTarget: null,
      optJobId: null, optStatus: 'idle', optError: null,
      slabResults: [], wulffData: null, viewerMode: 'structure',
    }),

  setJobStarted: (jobId) =>
    set({ jobId, jobStatus: 'running', jobError: null, progress: [], slabs: [],
          optJobId: null, optStatus: 'idle', slabResults: [], wulffData: null }),

  updateStatus: (data) =>
    set({
      jobStatus: data.status,
      progress: data.progress ?? [],
      slabs: data.results ?? [],
      completedIndices: data.completed_indices ?? 0,
      totalIndices: data.total_indices ?? 0,
      currentHkl: data.current_hkl ?? null,
      jobError: data.error ?? null,
    }),

  setOptStarted: (optJobId) =>
    set({ optJobId, optStatus: 'running', optError: null,
          optCompleted: 0, optTotal: 0, slabResults: [], wulffData: null }),

  updateOptStatus: (data) =>
    set({
      optStatus:    data.status,
      optPhase:     data.phase ?? null,
      optPhaseLabel: data.phase_label ?? null,
      optCompleted: data.completed ?? 0,
      optTotal:     data.total ?? 0,
      slabResults:  data.slab_results ?? [],
      wulffData:    data.wulff_data ?? null,
      optError:     data.error ?? null,
    }),

  selectSlab:     (i)    => set({ viewerTarget: i, viewerMode: 'structure' }),
  showBulk:       ()     => set({ viewerTarget: null, viewerMode: 'structure' }),
  setViewerMode:  (mode) => set({ viewerMode: mode }),
  toggleBonds:    ()     => set((s) => ({ showBonds: !s.showBonds })),
  toggleCell:     ()     => set((s) => ({ showUnitCell: !s.showUnitCell })),
  toggleEffects:  ()     => set((s) => ({ showEffects: !s.showEffects })),
}))
