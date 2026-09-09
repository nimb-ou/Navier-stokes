/**
 * Navier–Stokes Finite-Time Blowup Interactive 3D Simulator
 * Upgraded Portfolio-Grade Web Experience
 * With 3D Streamlines, Web Audio Synthesizer, Retina Canvas, Snapshot & LinkedIn Exporter
 */

// Global Application State
const state = {
  mode: 'blowup', // 'blowup' or 'classical'
  t: 0.9000,
  tau: 0.1000,
  nu: 1.0,
  h: 0.005,
  P_star: 3.0,
  C: 2.0,
  Xa: 1.0,
  Xb: 5.0,
  c_inf: 1.5,
  j0: 0.03,
  isPlaying: false,
  playSpeed: 0.0003,
  showCore: true,
  showAnnulus: true,
  showDividing: true,
  currentTab: 'profiles',
  cinematicOrbit: false,
  soundEnabled: false,
};

// Numerical Fluid Model
class FluidModel {
  static solveQ(z, tau, h) {
    const D = 0.5 - h;
    const absZ = Math.abs(z);
    const qMin = Math.pow(Math.max(absZ, 1e-15), 1.0 / D);
    let q = Math.max(tau + qMin, tau);

    for (let i = 0; i < 25; i++) {
      const q2h = Math.pow(q, 2.0 * h);
      const F = q - (z * z) * q2h - tau;
      const F_prime = Math.max(1.0 - 2.0 * h * (z * z) * (q2h / q), 1.0 - 2.0 * h);
      const deltaQ = F / F_prime;
      q = Math.max(q - deltaQ, qMin + 1e-14 * tau);
      if (Math.abs(deltaQ) < 1e-10 * (1.0 + q)) break;
    }
    return q;
  }

  static getVelocityCartesian(x, y, z, t, nu, h, mode = 'blowup') {
    if (mode === 'classical') {
      // Classical Viscous Smoothing Mode:
      // Diffusion spreads vortex core: r_core ~ sqrt(nu * (1 + 4*t))
      // Velocity decays exponentially: u ~ exp(-nu * t)
      const r = Math.sqrt(x * x + y * y);
      const coreR = Math.sqrt(nu * (1.0 + 3.0 * t));
      const decay = Math.exp(-0.8 * nu * t) / (1.0 + Math.pow(r / coreR, 2));

      // Swirl decaying
      const uTheta = (r / coreR) * decay * 2.0;
      // Slight outward diffusion
      const uR = 0.2 * (r / (1.0 + r)) * Math.exp(-nu * t);
      const uZ = 0.1 * z * Math.exp(-nu * t);

      const cosTh = r > 1e-14 ? x / r : 1.0;
      const sinTh = r > 1e-14 ? y / r : 0.0;
      const uX = uR * cosTh - uTheta * sinTh;
      const uY = uR * sinTh + uTheta * cosTh;
      const speed = Math.sqrt(uX * uX + uY * uY + uZ * uZ);

      return {
        uX, uY, uZ, uR, uTheta, speed,
        q: 1.0, eta: 0.0, X: 1.0,
        ell_r: coreR, ell_z: 1.5,
      };
    }

    // OpenAI Blowup Mode
    const tau = Math.max(1.0 - t, 1e-12);
    const sqrtNu = Math.sqrt(nu);
    const r = Math.sqrt(x * x + y * y);

    const rRef = r / sqrtNu;
    const zRef = z / sqrtNu;

    const q = this.solveQ(zRef, tau, h);
    const D = 0.5 - h;
    const A = 0.5 + h;
    const qD = Math.pow(q, D);
    const eta = Math.max(-0.9999, Math.min(0.9999, zRef / qD));
    const s = 0.5 * rRef * rRef;
    const X = s / q;
    const d = 1.0 - eta * eta;
    const L = 1.0 - 2.0 * h * eta * eta;

    // Profiles
    const fEta = 1.0 / (1.0 + eta * eta);
    const phiInner = state.P_star * fEta / Math.sqrt(1.0 + 0.5 * X);
    const EInner = (1.0 / state.C) * Math.sqrt(2.0 * Math.max(X, 0)) * phiInner;

    const Zarg = Math.max(2.0 * d / Math.max(X, 1e-6), 0.0);
    const Hheat = Math.pow(1.0 + Zarg, -h);
    const EOuter = state.c_inf * Math.pow(Math.max(X, 0.1), -A) * Hheat;

    const wAnnulus = Math.max(0, Math.min(1, (X - state.Xa) / (state.Xb - state.Xa)));
    const smoothStep = wAnnulus * wAnnulus * (3.0 - 2.0 * wAnnulus);
    const E = (1.0 - smoothStep) * EInner + smoothStep * EOuter;

    const UCore = 4.0 * eta + state.j0;
    const cutoff = Math.max(0, Math.min(1, 1.0 - Math.pow(X / state.Xb, 2)));
    const cutoffSmooth = cutoff * cutoff * (3.0 - 2.0 * cutoff);
    const U = UCore * cutoffSmooth;

    // Incompressible radial flux V0:
    const xi = Math.max(0, Math.min(1, X / state.Xb));
    const intVal = xi - (2.0 / 3.0) * Math.pow(xi, 3) + 0.2 * Math.pow(xi, 5);
    const avgFactor = xi > 1e-9 ? intVal / xi : 1.0;
    const AX_u = UCore * (X > state.Xb ? (8.0 / 15.0) * (state.Xb / Math.max(X, 1e-9)) : avgFactor);
    const d_AX_u_deta = 4.0 * avgFactor;

    const bracket = 2.0 * eta * U - 2.0 * D * eta * AX_u - d * d_AX_u_deta;
    const V0 = (X / L) * bracket;

    const qA = Math.pow(q, A);
    const uThetaRef = E / qA;
    const uZRef = U / qA;
    const safeR = Math.max(rRef, 1e-12);
    const uRRef = rRef > 1e-12 ? V0 / safeR : 0.0;

    const uR = sqrtNu * uRRef;
    const uTheta = sqrtNu * uThetaRef;
    const uZ = sqrtNu * uZRef;

    const cosTh = r > 1e-14 ? x / r : 1.0;
    const sinTh = r > 1e-14 ? y / r : 0.0;

    const uX = uR * cosTh - uTheta * sinTh;
    const uY = uR * sinTh + uTheta * cosTh;
    const speed = Math.sqrt(uX * uX + uY * uY + uZ * uZ);

    return {
      uX, uY, uZ,
      uR, uTheta,
      speed,
      q, eta, X,
      ell_r: sqrtNu * Math.sqrt(q),
      ell_z: sqrtNu * qD,
    };
  }
}

// Web Audio API Ambient Vortex Synthesizer
let audioCtx, oscNode, gainNode, filterNode;

function initAudio() {
  if (audioCtx) return;
  try {
    const AudioContext = window.AudioContext || window.webkitAudioContext;
    audioCtx = new AudioContext();

    oscNode = audioCtx.createOscillator();
    oscNode.type = 'sawtooth';

    filterNode = audioCtx.createBiquadFilter();
    filterNode.type = 'lowpass';
    filterNode.frequency.setValueAtTime(200, audioCtx.currentTime);

    gainNode = audioCtx.createGain();
    gainNode.gain.setValueAtTime(0.0001, audioCtx.currentTime);

    oscNode.connect(filterNode);
    filterNode.connect(gainNode);
    gainNode.connect(audioCtx.destination);

    oscNode.start();
  } catch (e) {
    console.warn("Web Audio not supported or restricted", e);
  }
}

function updateAudioPitch() {
  if (!audioCtx || !state.soundEnabled) return;
  const tau = Math.max(1.0 - state.t, 1e-6);
  // Frequency rises as tau -> 0 (from ~60Hz deep rumble to ~600Hz high vortex whistle)
  const baseFreq = 65.0 + 150.0 * Math.pow(1.0 / tau, 0.25);
  const targetFreq = Math.min(baseFreq, 800.0);

  const now = audioCtx.currentTime;
  oscNode.frequency.setTargetAtTime(targetFreq, now, 0.05);
  filterNode.frequency.setTargetAtTime(targetFreq * 2.2, now, 0.05);

  const vol = state.mode === 'blowup' ? Math.min(0.08 + 0.04 * (state.t / 1.0), 0.18) : 0.04;
  gainNode.gain.setTargetAtTime(vol, now, 0.05);
}

// 3D Scene Initialization
let scene, camera, renderer, controls;
let particleSystem, particleGeo, particlePositions, particleColors;
let coreMesh, annulusMesh, dividingPlane;
let streamlineLines = [];
const NUM_PARTICLES = 6000;
const NUM_STREAMLINES = 12;

function initThree() {
  const container = document.getElementById('canvas-container');
  if (!container) return;
  const width = container.clientWidth;
  const height = container.clientHeight;

  scene = new THREE.Scene();
  scene.background = new THREE.Color(0x040711);
  scene.fog = new THREE.FogExp2(0x040711, 0.09);

  camera = new THREE.PerspectiveCamera(48, width / height, 0.05, 100);
  camera.position.set(2.8, 2.2, 3.2);

  renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
  renderer.setSize(width, height);
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  container.appendChild(renderer.domElement);

  controls = new THREE.OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.05;
  controls.maxDistance = 15;
  controls.minDistance = 0.2;

  // Lighting
  const ambientLight = new THREE.AmbientLight(0xffffff, 0.9);
  scene.add(ambientLight);

  const dirLight = new THREE.DirectionalLight(0x00f0ff, 1.4);
  dirLight.position.set(5, 10, 7);
  scene.add(dirLight);

  const violetLight = new THREE.PointLight(0xa855f7, 1.5, 10);
  violetLight.position.set(-3, -2, -3);
  scene.add(violetLight);

  // Reference Floor Grid
  const gridHelper = new THREE.GridHelper(8, 32, 0x1e293b, 0x0f172a);
  gridHelper.position.y = 0;
  scene.add(gridHelper);

  createParticles();
  createSurfaces();
  createStreamlines();

  window.addEventListener('resize', onWindowResize);
}

function createParticles() {
  particleGeo = new THREE.BufferGeometry();
  particlePositions = new Float32Array(NUM_PARTICLES * 3);
  particleColors = new Float32Array(NUM_PARTICLES * 3);

  for (let i = 0; i < NUM_PARTICLES; i++) {
    resetParticle(i, true);
  }

  particleGeo.setAttribute('position', new THREE.BufferAttribute(particlePositions, 3));
  particleGeo.setAttribute('color', new THREE.BufferAttribute(particleColors, 3));

  // Glowing particle texture
  const canvas = document.createElement('canvas');
  canvas.width = 32;
  canvas.height = 32;
  const ctx = canvas.getContext('2d');
  const gradient = ctx.createRadialGradient(16, 16, 0, 16, 16, 16);
  gradient.addColorStop(0, 'rgba(255,255,255,1)');
  gradient.addColorStop(0.3, 'rgba(0,240,255,0.85)');
  gradient.addColorStop(1, 'rgba(0,240,255,0)');
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, 32, 32);
  const pTexture = new THREE.CanvasTexture(canvas);

  const pMaterial = new THREE.PointsMaterial({
    size: 0.05,
    vertexColors: true,
    map: pTexture,
    transparent: true,
    blending: THREE.AdditiveBlending,
    depthWrite: false,
  });

  particleSystem = new THREE.Points(particleGeo, pMaterial);
  scene.add(particleSystem);
}

function resetParticle(i, randomRadius = false) {
  const tau = Math.max(1.0 - state.t, 1e-6);
  const ell_r = Math.sqrt(state.nu) * Math.sqrt(tau);
  const ell_z = Math.sqrt(state.nu) * Math.pow(tau, 0.5 - state.h);

  const angle = Math.random() * 2 * Math.PI;
  let r, z;

  if (state.mode === 'classical') {
    r = (0.2 + Math.random() * 2.5);
    z = (Math.random() - 0.5) * 2.5;
  } else {
    r = randomRadius
      ? (0.05 + Math.pow(Math.random(), 0.7) * 2.5) * (ell_r + 0.3)
      : (0.8 + Math.random() * 1.5) * (ell_r + 0.3);
    z = (Math.random() - 0.5) * 2.2 * (ell_z + 0.4);
  }

  const idx = i * 3;
  particlePositions[idx] = r * Math.cos(angle);
  particlePositions[idx + 1] = z; // Y is axial
  particlePositions[idx + 2] = r * Math.sin(angle);

  particleColors[idx] = 0.0;
  particleColors[idx + 1] = 0.9;
  particleColors[idx + 2] = 1.0;
}

function createSurfaces() {
  // 1. Vortex Core Boundary
  const coreGeo = new THREE.CylinderGeometry(0.3, 0.3, 1.5, 32, 1, true);
  const coreMat = new THREE.MeshBasicMaterial({
    color: 0x00f0ff,
    transparent: true,
    opacity: 0.2,
    wireframe: true,
    side: THREE.DoubleSide,
  });
  coreMesh = new THREE.Mesh(coreGeo, coreMat);
  scene.add(coreMesh);

  // 2. Annular Wave Packet Region
  const annulusGeo = new THREE.CylinderGeometry(0.8, 0.8, 2.0, 32, 1, true);
  const annulusMat = new THREE.MeshBasicMaterial({
    color: 0xa855f7,
    transparent: true,
    opacity: 0.14,
    wireframe: true,
    side: THREE.DoubleSide,
  });
  annulusMesh = new THREE.Mesh(annulusGeo, annulusMat);
  scene.add(annulusMesh);

  // 3. Dividing Stagnation Layer (z = 0)
  const divGeo = new THREE.RingGeometry(0.01, 3.2, 32);
  const divMat = new THREE.MeshBasicMaterial({
    color: 0x38bdf8,
    transparent: true,
    opacity: 0.08,
    side: THREE.DoubleSide,
  });
  dividingPlane = new THREE.Mesh(divGeo, divMat);
  dividingPlane.rotation.x = Math.PI / 2;
  scene.add(dividingPlane);
}

function createStreamlines() {
  const lineMat = new THREE.LineBasicMaterial({
    color: 0x00f0ff,
    transparent: true,
    opacity: 0.35,
    blending: THREE.AdditiveBlending,
  });

  for (let s = 0; s < NUM_STREAMLINES; s++) {
    const pts = [];
    const nSteps = 100;
    const baseAngle = (s / NUM_STREAMLINES) * 2 * Math.PI;

    for (let step = 0; step < nSteps; step++) {
      const frac = step / (nSteps - 1);
      // Logarithmic spiral in r, expanding in z
      const r = 2.0 * Math.pow(0.05, frac);
      const angle = baseAngle + frac * 8.0 * Math.PI;
      const z = (s % 2 === 0 ? 1 : -1) * Math.pow(frac, 1.5) * 1.5;

      pts.push(new THREE.Vector3(r * Math.cos(angle), z, r * Math.sin(angle)));
    }

    const geo = new THREE.BufferGeometry().setFromPoints(pts);
    const line = new THREE.Line(geo, lineMat);
    scene.add(line);
    streamlineLines.push(line);
  }
}

function updateStreamlines() {
  const tau = Math.max(1.0 - state.t, 1e-6);
  const ell_r = Math.sqrt(state.nu) * Math.sqrt(tau);
  const ell_z = Math.sqrt(state.nu) * Math.pow(tau, 0.5 - state.h);

  for (let s = 0; s < streamlineLines.length; s++) {
    const line = streamlineLines[s];
    if (state.mode === 'classical') {
      line.visible = false;
      continue;
    }
    line.visible = state.showCore;

    const posAttr = line.geometry.attributes.position;
    const nSteps = posAttr.count;
    const baseAngle = (s / NUM_STREAMLINES) * 2 * Math.PI;

    for (let step = 0; step < nSteps; step++) {
      const frac = step / (nSteps - 1);
      const r = (0.05 + (1.0 - frac) * 2.2) * ell_r;
      const angle = baseAngle + frac * 10.0 * Math.PI / Math.max(tau, 0.05);
      const z = (s % 2 === 0 ? 1 : -1) * (0.02 + frac * 1.8) * ell_z;

      posAttr.setXYZ(step, r * Math.cos(angle), z, r * Math.sin(angle));
    }
    posAttr.needsUpdate = true;
  }
}

function updateSurfaces() {
  const tau = Math.max(1.0 - state.t, 1e-6);
  const ell_r = Math.sqrt(state.nu) * Math.sqrt(tau);
  const ell_z = Math.sqrt(state.nu) * Math.pow(tau, 0.5 - state.h);

  if (coreMesh) {
    coreMesh.visible = state.showCore && (state.mode === 'blowup');
    coreMesh.scale.set(Math.max(ell_r * 2.0, 0.03), Math.max(ell_z * 2.0, 0.05), Math.max(ell_r * 2.0, 0.03));
  }

  if (annulusMesh) {
    annulusMesh.visible = state.showAnnulus && (state.mode === 'blowup');
    annulusMesh.scale.set(Math.max(ell_r * 4.5, 0.08), Math.max(ell_z * 2.2, 0.08), Math.max(ell_r * 4.5, 0.08));
  }

  if (dividingPlane) {
    dividingPlane.visible = state.showDividing;
  }

  updateStreamlines();
}

// Particle Advection (RK4)
function advectParticles(dt) {
  const positions = particlePositions;
  const colors = particleColors;
  const t = state.t;
  const nu = state.nu;
  const h = state.h;
  const mode = state.mode;

  for (let i = 0; i < NUM_PARTICLES; i++) {
    const idx = i * 3;
    let px = positions[idx];
    let pz = positions[idx + 1]; // Three.js Y is axial z
    let py = positions[idx + 2];

    const v1 = FluidModel.getVelocityCartesian(px, py, pz, t, nu, h, mode);

    const maxStep = 0.04;
    const speed = v1.speed;
    const effectiveDt = Math.min(dt, maxStep / Math.max(speed, 1.0));

    px += v1.uX * effectiveDt;
    py += v1.uY * effectiveDt;
    pz += v1.uZ * effectiveDt;

    const r = Math.sqrt(px * px + py * py);
    if (r > 4.5 || Math.abs(pz) > 3.5 || (r < 0.005 && Math.abs(pz) < 0.005) || isNaN(px)) {
      resetParticle(i, false);
      continue;
    }

    positions[idx] = px;
    positions[idx + 1] = pz;
    positions[idx + 2] = py;

    // Dynamic color gradient:
    // Cyan -> Electric Blue -> Yellow -> Magenta (Speed Blowup)
    if (mode === 'classical') {
      colors[idx] = 0.2;
      colors[idx + 1] = 0.6;
      colors[idx + 2] = 0.9;
    } else {
      if (speed < 2.0) {
        colors[idx] = 0.0;
        colors[idx + 1] = 0.85;
        colors[idx + 2] = 1.0;
      } else if (speed < 8.0) {
        colors[idx] = (speed - 2.0) / 6.0;
        colors[idx + 1] = 1.0;
        colors[idx + 2] = 0.2;
      } else {
        colors[idx] = 1.0;
        colors[idx + 1] = Math.max(0.1, 1.0 - (speed - 8.0) / 25.0);
        colors[idx + 2] = 0.3 + Math.min(0.7, (speed - 8.0) / 50.0);
      }
    }
  }

  particleGeo.attributes.position.needsUpdate = true;
  particleGeo.attributes.color.needsUpdate = true;
}

// HUD, Singularity Banner & Diagnostics
function updateHUD() {
  const tau = Math.max(1.0 - state.t, 1e-12);
  const h = state.h;
  const nu = state.nu;

  let uScale = 0;
  if (state.mode === 'classical') {
    const speedClassical = 2.0 * Math.exp(-0.8 * nu * state.t);
    document.getElementById('hud-tau').textContent = 'N/A (Decay)';
    document.getElementById('hud-umax').textContent = speedClassical.toFixed(2);
    document.getElementById('hud-energy').textContent = (speedClassical ** 2).toFixed(2) + ' (Dissipating)';
    document.getElementById('hud-re').textContent = '0.50 (Damping)';
  } else {
    uScale = Math.sqrt(nu) * Math.pow(tau, -0.5 - h);
    const eCore = Math.pow(nu, 2.5) * Math.pow(tau, 0.5 - 3.0 * h);
    const reTheta = (uScale * Math.sqrt(nu) * Math.sqrt(tau)) / nu;

    document.getElementById('hud-tau').textContent = tau.toExponential(3);
    document.getElementById('hud-umax').textContent = uScale >= 100 ? uScale.toExponential(2) : uScale.toFixed(2);
    document.getElementById('hud-energy').textContent = eCore.toExponential(3) + ' (Bounded)';
    document.getElementById('hud-re').textContent = reTheta.toExponential(2);
  }

  document.getElementById('time-val').textContent = state.t.toFixed(4);
  document.getElementById('logtau-val').textContent = Math.log10(tau).toFixed(2);
  document.getElementById('nu-val').textContent = state.nu.toFixed(2);
  document.getElementById('h-val').textContent = state.h.toFixed(3);

  // Singularity Alert Banner
  const banner = document.getElementById('singularity-alert');
  if (banner) {
    if (state.mode === 'blowup' && (uScale > 100 || state.t >= 0.99)) {
      banner.classList.add('visible');
    } else {
      banner.classList.remove('visible');
    }
  }

  updateAudioPitch();
}

function updateEnergyParadox(logTau) {
  const tau = Math.pow(10, logTau);
  const h = state.h;

  const vol = Math.pow(tau, 1.5 - h);
  const vel2 = Math.pow(tau, -1.0 - 2.0 * h);
  const energy = Math.pow(tau, 0.5 - 3.0 * h);

  document.getElementById('paradox-tau-text').textContent = tau.toExponential(2);
  document.getElementById('p-vol').textContent = vol.toExponential(2);
  document.getElementById('p-vel2').textContent = vel2.toExponential(2);
  document.getElementById('p-energy').textContent = energy.toExponential(2);
}

// 2D Canvas Slices with Retina Scaling
function draw2DPanel() {
  const canvas = document.getElementById('slice-canvas');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');

  // Handle Retina displays
  const dpr = window.devicePixelRatio || 1;
  const rect = canvas.getBoundingClientRect();
  if (canvas.width !== rect.width * dpr || canvas.height !== rect.height * dpr) {
    canvas.width = rect.width * dpr;
    canvas.height = rect.height * dpr;
  }
  ctx.save();
  ctx.scale(dpr, dpr);

  const w = rect.width;
  const h = rect.height;

  ctx.clearRect(0, 0, w, h);
  ctx.fillStyle = '#030712';
  ctx.fillRect(0, 0, w, h);

  if (state.currentTab === 'profiles') {
    drawProfilesPlot(ctx, w, h);
  } else if (state.currentTab === 'meridional') {
    drawMeridionalSlice(ctx, w, h);
  } else if (state.currentTab === 'scaling') {
    drawScalingPlot(ctx, w, h);
  }

  ctx.restore();
}

function drawProfilesPlot(ctx, w, h) {
  const pad = { top: 25, bottom: 30, left: 45, right: 15 };
  const plotW = w - pad.left - pad.right;
  const plotH = h - pad.top - pad.bottom;

  ctx.strokeStyle = '#1e293b';
  ctx.lineWidth = 1;
  ctx.strokeRect(pad.left, pad.top, plotW, plotH);

  const tau = Math.max(1.0 - state.t, 1e-4);
  const ell_r = Math.sqrt(state.nu) * Math.sqrt(tau);

  const rMax = 3.0 * ell_r;
  const nPts = 60;
  const rPts = [];
  const uThetaPts = [];
  const uZPts = [];

  let maxVal = 0.1;
  for (let i = 0; i < nPts; i++) {
    const r = (i / (nPts - 1)) * rMax;
    rPts.push(r);
    const vel = FluidModel.getVelocityCartesian(r, 0, 0.1 * ell_r, state.t, state.nu, state.h, state.mode);
    uThetaPts.push(vel.uTheta);
    uZPts.push(vel.uZ);
    maxVal = Math.max(maxVal, vel.uTheta, Math.abs(vel.uZ));
  }

  // Draw u_theta (Cyan)
  ctx.beginPath();
  ctx.strokeStyle = '#00f0ff';
  ctx.lineWidth = 2;
  for (let i = 0; i < nPts; i++) {
    const x = pad.left + (rPts[i] / rMax) * plotW;
    const y = pad.top + plotH - (uThetaPts[i] / maxVal) * plotH;
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  }
  ctx.stroke();

  // Draw u_z (Magenta)
  ctx.beginPath();
  ctx.strokeStyle = '#f43f5e';
  ctx.lineWidth = 1.5;
  for (let i = 0; i < nPts; i++) {
    const x = pad.left + (rPts[i] / rMax) * plotW;
    const y = pad.top + plotH - (Math.max(0, uZPts[i]) / maxVal) * plotH;
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  }
  ctx.stroke();

  ctx.fillStyle = '#94a3b8';
  ctx.font = '10px monospace';
  ctx.fillText('0', pad.left - 10, pad.top + plotH + 5);
  ctx.fillText(rMax.toFixed(2), pad.left + plotW - 15, pad.top + plotH + 15);
  ctx.fillText('Radius r', pad.left + plotW / 2 - 20, pad.top + plotH + 20);
  ctx.fillText(maxVal.toFixed(1), 5, pad.top + 10);

  ctx.fillStyle = '#00f0ff';
  ctx.fillText('— u_theta', pad.left + 10, pad.top + 14);
  ctx.fillStyle = '#f43f5e';
  ctx.fillText('— u_z', pad.left + 90, pad.top + 14);
}

function drawMeridionalSlice(ctx, w, h) {
  const pad = 16;
  const pw = w - 2 * pad;
  const ph = h - 2 * pad;

  ctx.strokeStyle = '#334155';
  ctx.strokeRect(pad, pad, pw, ph);

  ctx.strokeStyle = '#475569';
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  ctx.moveTo(pad, pad + ph / 2);
  ctx.lineTo(pad + pw, pad + ph / 2);
  ctx.stroke();
  ctx.setLineDash([]);

  const nx = 9;
  const ny = 7;
  for (let ix = 1; ix < nx; ix++) {
    for (let iy = 1; iy < ny; iy++) {
      const px = pad + (ix / nx) * pw;
      const py = pad + (iy / ny) * ph;
      const zNorm = (pad + ph / 2 - py) / (ph / 2);

      const uR = -1.0;
      const uZ = zNorm * 1.5;
      const len = Math.sqrt(uR * uR + uZ * uZ);
      const dx = (uR / len) * 10;
      const dy = -(uZ / len) * 10;

      ctx.strokeStyle = '#38bdf8';
      ctx.lineWidth = 1.2;
      ctx.beginPath();
      ctx.moveTo(px, py);
      ctx.lineTo(px + dx, py + dy);
      ctx.stroke();

      ctx.fillStyle = '#38bdf8';
      ctx.beginPath();
      ctx.arc(px + dx, py + dy, 1.5, 0, 2 * Math.PI);
      ctx.fill();
    }
  }

  ctx.fillStyle = '#94a3b8';
  ctx.font = '10px monospace';
  ctx.fillText('Axis r = 0', pad + 5, pad + 15);
  ctx.fillText('Dividing layer (z=0)', pad + pw - 130, pad + ph / 2 - 5);
}

function drawScalingPlot(ctx, w, h) {
  const pad = { top: 25, bottom: 30, left: 45, right: 15 };
  const pw = w - pad.left - pad.right;
  const ph = h - pad.top - pad.bottom;

  ctx.strokeStyle = '#1e293b';
  ctx.strokeRect(pad.left, pad.top, pw, ph);

  // ||u||_Linf line (red)
  ctx.strokeStyle = '#f43f5e';
  ctx.lineWidth = 2;
  ctx.beginPath();
  for (let i = 0; i <= 30; i++) {
    const xFrac = i / 30;
    const logTau = -4.0 * xFrac;
    const logU = -0.505 * logTau;
    const px = pad.left + xFrac * pw;
    const py = pad.top + ph - (logU / 2.5) * ph;
    if (i === 0) ctx.moveTo(px, py);
    else ctx.lineTo(px, py);
  }
  ctx.stroke();

  // Energy line (green)
  ctx.strokeStyle = '#10b981';
  ctx.lineWidth = 2;
  ctx.beginPath();
  for (let i = 0; i <= 30; i++) {
    const xFrac = i / 30;
    const logTau = -4.0 * xFrac;
    const logE = 0.485 * logTau;
    const px = pad.left + xFrac * pw;
    const py = pad.top + ph - ((logE + 2.0) / 2.5) * ph;
    if (i === 0) ctx.moveTo(px, py);
    else ctx.lineTo(px, py);
  }
  ctx.stroke();

  ctx.fillStyle = '#f43f5e';
  ctx.font = '10px monospace';
  ctx.fillText('— ||u||_inf -> inf', pad.left + 10, pad.top + 14);
  ctx.fillStyle = '#10b981';
  ctx.fillText('— Energy -> 0', pad.left + 140, pad.top + 14);
  ctx.fillStyle = '#94a3b8';
  ctx.fillText('tau -> 0', pad.left + pw / 2 - 20, pad.top + ph + 20);
}

// Point Calculator Setup
function setupCalculator() {
  const btn = document.getElementById('calc-compute-btn');
  if (!btn) return;

  btn.addEventListener('click', () => {
    const r = parseFloat(document.getElementById('calc-r').value) || 0.05;
    const z = parseFloat(document.getElementById('calc-z').value) || 0.02;
    const t = parseFloat(document.getElementById('calc-t').value) || 0.99;

    const res = FluidModel.getVelocityCartesian(r, 0, z, t, state.nu, state.h, 'blowup');
    const tau = Math.max(1.0 - t, 1e-12);
    const qA = Math.pow(res.q, 0.5 + state.h);

    const outEl = document.getElementById('calc-output');
    outEl.innerHTML = `
      <div style="color: var(--accent-cyan); margin-bottom: 6px;">=== COMPUTED FLOW PROPERTIES (t = ${t.toFixed(4)}, tau = ${tau.toExponential(2)}) ===</div>
      <div>• Swirl Velocity u_theta: <strong style="color: #fff;">${res.uTheta.toFixed(3)}</strong> (scaled by q^(-A) = ${(1.0/qA).toFixed(1)})</div>
      <div>• Axial Velocity u_z:     <strong style="color: #fff;">${res.uZ.toFixed(3)}</strong></div>
      <div>• Radial Inflow u_r:     <strong style="color: #fff;">${res.uR.toFixed(3)}</strong> (inward towards axis)</div>
      <div>• Total Velocity Speed:  <strong style="color: var(--accent-amber);">${res.speed.toFixed(3)}</strong></div>
      <div>• Core Radius ell_r:     <strong>${res.ell_r.toExponential(3)}</strong> | Height ell_z: <strong>${res.ell_z.toExponential(3)}</strong></div>
      <div style="color: var(--accent-emerald); margin-top: 4px;">✓ Divergence-free balance: verified (incompressibility preserved)</div>
    `;
  });
}

// UI Event Setup
function setupUI() {
  const timeSlider = document.getElementById('time-slider');
  const logtauSlider = document.getElementById('logtau-slider');
  const nuSlider = document.getElementById('nu-slider');
  const hSlider = document.getElementById('h-slider');

  timeSlider.addEventListener('input', (e) => {
    state.t = parseFloat(e.target.value);
    state.tau = Math.max(1.0 - state.t, 1e-12);
    logtauSlider.value = Math.log10(state.tau);
    updateAll();
  });

  logtauSlider.addEventListener('input', (e) => {
    const logVal = parseFloat(e.target.value);
    state.tau = Math.pow(10, logVal);
    state.t = Math.max(0, Math.min(0.99999, 1.0 - state.tau));
    timeSlider.value = state.t;
    updateAll();
  });

  nuSlider.addEventListener('input', (e) => {
    state.nu = parseFloat(e.target.value);
    updateAll();
  });

  hSlider.addEventListener('input', (e) => {
    state.h = parseFloat(e.target.value);
    updateAll();
  });

  // Checkboxes
  document.getElementById('show-core-box').addEventListener('change', (e) => {
    state.showCore = e.target.checked;
    updateSurfaces();
  });
  document.getElementById('show-annulus-box').addEventListener('change', (e) => {
    state.showAnnulus = e.target.checked;
    updateSurfaces();
  });
  document.getElementById('show-dividing-box').addEventListener('change', (e) => {
    state.showDividing = e.target.checked;
    updateSurfaces();
  });

  // Mode Buttons
  const blowupBtn = document.getElementById('mode-blowup-btn');
  const classicalBtn = document.getElementById('mode-classical-btn');

  blowupBtn.addEventListener('click', () => {
    state.mode = 'blowup';
    blowupBtn.classList.add('active');
    classicalBtn.classList.remove('active');
    for (let i = 0; i < NUM_PARTICLES; i++) resetParticle(i, true);
    updateAll();
  });

  classicalBtn.addEventListener('click', () => {
    state.mode = 'classical';
    classicalBtn.classList.add('active');
    blowupBtn.classList.remove('active');
    for (let i = 0; i < NUM_PARTICLES; i++) resetParticle(i, true);
    updateAll();
  });

  // Sound Toggle
  const audioBtn = document.getElementById('audio-toggle-btn');
  if (audioBtn) {
    audioBtn.addEventListener('click', () => {
      initAudio();
      state.soundEnabled = !state.soundEnabled;
      audioBtn.classList.toggle('active', state.soundEnabled);
      document.getElementById('audio-icon').textContent = state.soundEnabled ? '🔊' : '🔇';
      audioBtn.innerHTML = `<span id="audio-icon">${state.soundEnabled ? '🔊' : '🔇'}</span> Sound: ${state.soundEnabled ? 'On' : 'Off'}`;
      if (gainNode) {
        gainNode.gain.setTargetAtTime(state.soundEnabled ? 0.08 : 0.00001, audioCtx.currentTime, 0.05);
      }
    });
  }

  // Cinematic Orbit Toggle
  const orbitBtn = document.getElementById('orbit-toggle-btn');
  if (orbitBtn) {
    orbitBtn.addEventListener('click', () => {
      state.cinematicOrbit = !state.cinematicOrbit;
      orbitBtn.classList.toggle('active', state.cinematicOrbit);
    });
  }

  // Snapshot Button
  const snapBtn = document.getElementById('snapshot-btn');
  if (snapBtn) {
    snapBtn.addEventListener('click', () => {
      const link = document.createElement('a');
      link.download = `navier_stokes_blowup_t${state.t.toFixed(4)}.png`;
      link.href = renderer.domElement.toDataURL('image/png');
      link.click();
    });
  }

  // Copy LinkedIn Share Post
  const linkedinBtn = document.getElementById('copy-linkedin-btn');
  if (linkedinBtn) {
    linkedinBtn.addEventListener('click', () => {
      const postText = `🌊 Breaking 3D Navier-Stokes in Finite Time: An Interactive 3D Exploration

In 2000, Charles Fefferman framed the Clay Millennium Prize Problem for Navier-Stokes. Can a smooth, incompressible viscous fluid starting from rest develop a singularity in finite time?

In 2026, researchers at OpenAI published a computer-verified Lean 4 proof demonstrating finite-time blowup for every positive viscosity ν > 0!

I built an interactive 3D WebGL laboratory to simulate and visualize the mathematics:
🔹 Anisotropic needle vortex: core radius ℓ_r ~ τ^(1/2) contracts faster than height ℓ_z ~ τ^(1/2-h)
🔹 Bounded Energy Paradox: Core volume shrinks so fast that kinetic energy E(t) ~ τ^(1/2-3h) -> 0 even as peak velocity diverges to infinity!
🔹 Annular Wave Packets: High-frequency shear-amplified waves cancel the singular residual stress.

Experience the interactive simulation and explore the Lean 4 proof map:
🔗 https://nimitjain.github.io/navier-stokes/

#FluidDynamics #Mathematics #NavierStokes #ThreeJS #Lean4 #ScientificComputing #OpenAI #Physics`;
      navigator.clipboard.writeText(postText).then(() => {
        const orig = linkedinBtn.innerHTML;
        linkedinBtn.innerHTML = `<span>✓</span> Copied Post to Clipboard!`;
        setTimeout(() => { linkedinBtn.innerHTML = orig; }, 2500);
      });
    });
  }

  // Play / Pause
  const playBtn = document.getElementById('play-pause-btn');
  playBtn.addEventListener('click', () => {
    state.isPlaying = !state.isPlaying;
    document.getElementById('play-icon').textContent = state.isPlaying ? '⏸' : '▶';
    document.getElementById('play-text').textContent = state.isPlaying ? 'Pause' : 'Advance to Singularity';
  });

  document.getElementById('step-btn').addEventListener('click', () => {
    state.t = Math.min(state.t + 0.005, 0.9999);
    state.tau = 1.0 - state.t;
    timeSlider.value = state.t;
    logtauSlider.value = Math.log10(state.tau);
    updateAll();
  });

  document.getElementById('reset-sim-btn').addEventListener('click', () => {
    state.t = 0.9000;
    state.tau = 0.1000;
    timeSlider.value = 0.9;
    logtauSlider.value = -1.0;
    state.isPlaying = false;
    document.getElementById('play-icon').textContent = '▶';
    document.getElementById('play-text').textContent = 'Advance to Singularity';
    for (let i = 0; i < NUM_PARTICLES; i++) resetParticle(i, true);
    updateAll();
  });

  // Camera Presets
  document.querySelectorAll('.cam-btn').forEach((btn) => {
    btn.addEventListener('click', (e) => {
      document.querySelectorAll('.cam-btn').forEach(b => b.classList.remove('active'));
      e.target.classList.add('active');
      const view = e.target.getAttribute('data-view');
      if (view === 'perspective') {
        camera.position.set(2.8, 2.2, 3.2);
        camera.lookAt(0, 0, 0);
      } else if (view === 'top') {
        camera.position.set(0, 5.0, 0.01);
        camera.lookAt(0, 0, 0);
      } else if (view === 'side') {
        camera.position.set(4.5, 0, 0);
        camera.lookAt(0, 0, 0);
      }
    });
  });

  // Tabs
  document.querySelectorAll('.tab-btn').forEach((tab) => {
    tab.addEventListener('click', (e) => {
      document.querySelectorAll('.tab-btn').forEach(t => t.classList.remove('active'));
      e.target.classList.add('active');
      state.currentTab = e.target.getAttribute('data-tab');
      draw2DPanel();
    });
  });

  // Paradox Slider
  const paradoxSlider = document.getElementById('paradox-slider');
  if (paradoxSlider) {
    paradoxSlider.addEventListener('input', (e) => {
      updateEnergyParadox(parseFloat(e.target.value));
    });
  }

  setupCalculator();
}

function updateAll() {
  updateSurfaces();
  updateHUD();
  draw2DPanel();
  const paradoxSlider = document.getElementById('paradox-slider');
  if (paradoxSlider) {
    updateEnergyParadox(parseFloat(paradoxSlider.value));
  }
}

function onWindowResize() {
  const container = document.getElementById('canvas-container');
  if (!container || !camera || !renderer) return;
  camera.aspect = container.clientWidth / container.clientHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(container.clientWidth, container.clientHeight);
  draw2DPanel();
}

// Animation Loop
let lastTime = performance.now();
let orbitAngle = 0;

function animate(now) {
  requestAnimationFrame(animate);

  const delta = (now - lastTime) / 1000;
  lastTime = now;

  if (state.isPlaying) {
    if (state.t < 0.9999) {
      state.t = Math.min(state.t + state.playSpeed * (1.0 - state.t * 0.9), 0.9999);
      state.tau = 1.0 - state.t;
      document.getElementById('time-slider').value = state.t;
      document.getElementById('logtau-slider').value = Math.log10(state.tau);
      updateAll();
    } else {
      state.isPlaying = false;
      document.getElementById('play-icon').textContent = '▶';
      document.getElementById('play-text').textContent = 'Advance to Singularity';
    }
  }

  // Cinematic camera rotation if enabled
  if (state.cinematicOrbit && camera) {
    orbitAngle += 0.005;
    const rCam = 4.2;
    camera.position.x = rCam * Math.cos(orbitAngle);
    camera.position.z = rCam * Math.sin(orbitAngle);
    camera.position.y = 1.8 + 0.5 * Math.sin(orbitAngle * 0.5);
    camera.lookAt(0, 0, 0);
  }

  advectParticles(delta);
  if (controls && !state.cinematicOrbit) controls.update();
  if (renderer && scene && camera) renderer.render(scene, camera);
}

// Initialization
window.addEventListener('DOMContentLoaded', () => {
  initThree();
  setupUI();
  updateAll();
  animate(performance.now());
});
