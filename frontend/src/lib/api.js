const BASE = '/api'

async function request(path, options) {
  const res = await fetch(BASE + path, options)
  if (!res.ok) {
    const text = await res.text()
    throw new Error(text || `HTTP ${res.status}`)
  }
  return res.json()
}

export function uploadCif(file) {
  const fd = new FormData()
  fd.append('file', file)
  return request('/upload', { method: 'POST', body: fd })
}

export function startGenerate(sessionId, params) {
  return request('/generate', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ session_id: sessionId, ...params }),
  })
}

export function pollStatus(jobId) {
  return request(`/status/${jobId}`)
}

export function downloadUrl(jobId, slabIndex) {
  return `${BASE}/download/${jobId}/${slabIndex}`
}

export function startOptimize(sessionId, jobId) {
  return request('/optimize', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ session_id: sessionId, job_id: jobId }),
  })
}

export function pollOptStatus(optJobId) {
  return request(`/opt_status/${optJobId}`)
}
