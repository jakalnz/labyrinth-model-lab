// Builds ear-data.json for the Labyrinth Model Lab: every IE-Map piece in one shared
// frame (HeadFrame mm: +X anterior, +Y left, +Z superior; origin = utricle centroid,
// right ear), plus a pre-built "unified skin" (voxel union of the membranous pieces ->
// blur -> marching cubes).
//
// Usage: node scripts/build-data.mjs <out.json> [voxelMm=0.06] [blurVoxels=1.1]
// Set IEMAP_DATA to the IE-Map dataset's data/ folder (default: ./IEMap_data_v_1_0/data/).
import { parseVtkPolydata } from './vtk.mjs';
import { edgeTable, triTable } from 'three/examples/jsm/objects/MarchingCubes.js';
import { writeFileSync } from 'node:fs';

const D = (process.env.IEMAP_DATA ?? 'IEMap_data_v_1_0/data').replace(/\/?$/, '/');
const OUT = process.argv[2];
// RAS (3D Slicer: +X right, +Y anterior, +Z superior) -> HeadFrame, in mm.
const rasToHead = (x, y, z) => [y, -x, z];

const ut = parseVtkPolydata(D + 'ls_Hsapiens_utricle.vtk').points;
const anchor = [0, 0, 0];
for (let i = 0; i < ut.length; i += 3) {
  const h = rasToHead(ut[i], ut[i + 1], ut[i + 2]);
  for (let k = 0; k < 3; k++) anchor[k] += h[k];
}
for (let k = 0; k < 3; k++) anchor[k] /= ut.length / 3;

function load(file) {
  const { points, indices } = parseVtkPolydata(D + file);
  const p = new Float32Array(points.length);
  for (let i = 0; i < points.length; i += 3) {
    const h = rasToHead(points[i], points[i + 1], points[i + 2]);
    p[i] = h[0] - anchor[0];
    p[i + 1] = h[1] - anchor[1];
    p[i + 2] = h[2] - anchor[2];
  }
  return { points: p, indices: Uint32Array.from(indices) };
}

// id, file, label, group, inUnion
const PIECES = [
  ['post-duct', 'ls_Hsapiens_Sp.vtk', 'Posterior duct', 'posterior', true],
  ['ant-duct', 'ls_Hsapiens_Sa.vtk', 'Anterior duct', 'anterior', true],
  ['hor-duct', 'ls_Hsapiens_Sl.vtk', 'Horizontal duct', 'horizontal', true],
  ['post-ampulla', 'ls_Hsapiens_Ap.vtk', 'Posterior ampulla', 'posterior', true],
  ['ant-ampulla', 'ls_Hsapiens_Aa.vtk', 'Anterior ampulla', 'anterior', true],
  ['hor-ampulla', 'ls_Hsapiens_Al.vtk', 'Horizontal ampulla', 'horizontal', true],
  ['post-cupula', 'ls_Hsapiens_Post_Cup_Wall.vtk', 'Posterior cupula', 'cupula', false],
  ['ant-cupula', 'ls_Hsapiens_Ant_Cup_Wall.vtk', 'Anterior cupula', 'cupula', false],
  ['hor-cupula', 'ls_Hsapiens_Lat_Cup_Wall.vtk', 'Horizontal cupula', 'cupula', false],
  ['post-cup-ut', 'ls_Hsapiens_Post_Cup_Ut.vtk', 'Posterior ampulla to utricle wall', 'posterior', true],
  ['ant-cup-ut', 'ls_Hsapiens_Ant_Cup_Ut.vtk', 'Anterior ampulla to utricle wall', 'anterior', true],
  ['hor-cup-ut', 'ls_Hsapiens_Lat_Cup_Ut.vtk', 'Horizontal ampulla to utricle wall', 'horizontal', true],
  ['post-cup-duct', 'ls_Hsapiens_Post_Cup_Duct.vtk', 'Posterior ampulla to duct wall', 'posterior', true],
  ['ant-cup-duct', 'ls_Hsapiens_Ant_Cup_Duct.vtk', 'Anterior ampulla to duct wall', 'anterior', true],
  ['hor-cup-duct', 'ls_Hsapiens_Lat_Cup_Duct.vtk', 'Horizontal ampulla to duct wall', 'horizontal', true],
  ['crus-a', 'ls_Hsapiens_CCa.vtk', 'Common crus, anterior half', 'crus', true],
  ['crus-p', 'ls_Hsapiens_CCp.vtk', 'Common crus, posterior half', 'crus', true],
  ['ucp', 'ls_Hsapiens_UCp.vtk', 'Utricle to crus junction (UCp)', 'crus', true],
  ['sinus', 'ls_Hsapiens_SC.vtk', 'Horizontal canal utricular sinus', 'horizontal', true],
  ['utricle', 'ls_Hsapiens_utricle.vtk', 'Utricle', 'otolith', true],
  ['saccule', 'mesh_david_sacculus.vtk', 'Saccule', 'otolith', true],
  ['cochlea', 'mesh_david_cochlea.vtk', 'Cochlear duct', 'cochlea', false],
  ['envelope', 'segT2_th1.07_Otsu.vtk', 'Fluid envelope (T2 MRI)', 'housing', false],
];

const meshes = {};
for (const [id, file] of PIECES) meshes[id] = load(file);

// --- unified skin ---------------------------------------------------------------------
const H = Number(process.argv[3] ?? 0.06); // voxel size mm
const SIGMA = Number(process.argv[4] ?? 1.1); // blur, voxels
const unionIds = PIECES.filter((p) => p[4]).map((p) => p[0]);
const mn = [1e9, 1e9, 1e9];
const mx = [-1e9, -1e9, -1e9];
for (const id of unionIds) {
  const p = meshes[id].points;
  for (let i = 0; i < p.length; i += 3)
    for (let k = 0; k < 3; k++) {
      mn[k] = Math.min(mn[k], p[i + k]);
      mx[k] = Math.max(mx[k], p[i + k]);
    }
}
for (let k = 0; k < 3; k++) {
  mn[k] -= 6 * H;
  mx[k] += 6 * H;
}
const [NX, NY, NZ] = mn.map((v, k) => Math.ceil((mx[k] - v) / H) + 1);
const idx3 = (x, y, z) => x + NX * (y + NY * z);
const occ = new Float32Array(NX * NY * NZ);
console.log('grid', [NX, NY, NZ], ((NX * NY * NZ) / 1e6).toFixed(1) + 'M voxels');

// Per-mesh x-ray parity fill (every piece is a closed surface), OR'd into occ.
for (const id of unionIds) {
  const { points: p, indices: t } = meshes[id];
  const cols = new Map();
  for (let f = 0; f < t.length; f += 3) {
    const a = t[f] * 3;
    const b = t[f + 1] * 3;
    const c = t[f + 2] * 3;
    const ax = (p[a] - mn[0]) / H, ay = (p[a + 1] - mn[1]) / H, az = (p[a + 2] - mn[2]) / H;
    const bx = (p[b] - mn[0]) / H, by = (p[b + 1] - mn[1]) / H, bz = (p[b + 2] - mn[2]) / H;
    const cx = (p[c] - mn[0]) / H, cy = (p[c + 1] - mn[1]) / H, cz = (p[c + 2] - mn[2]) / H;
    const y0 = Math.ceil(Math.min(ay, by, cy)), y1 = Math.floor(Math.max(ay, by, cy));
    const z0 = Math.ceil(Math.min(az, bz, cz)), z1 = Math.floor(Math.max(az, bz, cz));
    const det = (by - ay) * (cz - az) - (cy - ay) * (bz - az);
    if (Math.abs(det) < 1e-12) continue;
    for (let y = y0; y <= y1; y++)
      for (let z = z0; z <= z1; z++) {
        const py = y + 1e-4, pz = z + 2e-4; // jitter avoids exact edge/vertex hits
        const u = ((py - ay) * (cz - az) - (cy - ay) * (pz - az)) / det;
        const v = ((by - ay) * (pz - az) - (py - ay) * (bz - az)) / det;
        if (u < 0 || v < 0 || u + v > 1) continue;
        const x = ax + u * (bx - ax) + v * (cx - ax);
        const key = y * NZ + z;
        let arr = cols.get(key);
        if (!arr) cols.set(key, (arr = []));
        arr.push(x);
      }
  }
  for (const [key, xs] of cols) {
    if (xs.length < 2) continue;
    xs.sort((m, n) => m - n);
    const y = Math.floor(key / NZ), z = key % NZ;
    for (let i = 0; i + 1 < xs.length; i += 2)
      for (let x = Math.ceil(xs[i]); x <= Math.floor(xs[i + 1]); x++) occ[idx3(x, y, z)] = 1;
  }
}

function blur(src, sigma) {
  const r = Math.ceil(sigma * 2.5);
  const w = [];
  let s = 0;
  for (let i = -r; i <= r; i++) {
    const v = Math.exp(-(i * i) / (2 * sigma * sigma));
    w.push(v);
    s += v;
  }
  for (let i = 0; i < w.length; i++) w[i] /= s;
  let a = src;
  let b = new Float32Array(src.length);
  const strides = [1, NX, NX * NY];
  const lens = [NX, NY, NZ];
  for (let axis = 0; axis < 3; axis++) {
    const st = strides[axis], L = lens[axis];
    for (let i = 0; i < a.length; i++) {
      const coord = Math.floor(i / st) % L;
      let acc = 0;
      for (let k = -r; k <= r; k++) {
        const c = coord + k;
        if (c < 0 || c >= L) continue;
        acc += a[i + k * st] * w[k + r];
      }
      b[i] = acc;
    }
    [a, b] = [b, a];
  }
  return a;
}

function marchingCubes(field, iso) {
  const verts = [];
  const tris = [];
  const edgeVert = new Map();
  const corner = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]];
  const edges = [[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]];
  const vals = new Float32Array(8);
  const ev = new Int32Array(12);
  for (let z = 0; z < NZ - 1; z++)
    for (let y = 0; y < NY - 1; y++)
      for (let x = 0; x < NX - 1; x++) {
        let ci = 0;
        for (let c = 0; c < 8; c++) {
          vals[c] = field[idx3(x + corner[c][0], y + corner[c][1], z + corner[c][2])];
          if (vals[c] < iso) ci |= 1 << c;
        }
        const em = edgeTable[ci];
        if (!em) continue;
        for (let e = 0; e < 12; e++) {
          if (!(em & (1 << e))) continue;
          const [c0, c1] = edges[e];
          const g0 = [x + corner[c0][0], y + corner[c0][1], z + corner[c0][2]];
          const g1 = [x + corner[c1][0], y + corner[c1][1], z + corner[c1][2]];
          const k0 = idx3(g0[0], g0[1], g0[2]);
          const k1 = idx3(g1[0], g1[1], g1[2]);
          const lo = Math.min(k0, k1), d = Math.abs(k1 - k0);
          const key = lo * 4 + (d === 1 ? 0 : d === NX ? 1 : 2);
          let vi = edgeVert.get(key);
          if (vi === undefined) {
            const t = (iso - vals[c0]) / (vals[c1] - vals[c0]);
            vi = verts.length / 3;
            for (let k = 0; k < 3; k++) verts.push(mn[k] + (g0[k] + t * (g1[k] - g0[k])) * H);
            edgeVert.set(key, vi);
          }
          ev[e] = vi;
        }
        const o = ci << 4;
        for (let i = 0; triTable[o + i] !== -1; i += 3)
          tris.push(ev[triTable[o + i]], ev[triTable[o + i + 2]], ev[triTable[o + i + 1]]);
      }
  return { points: Float32Array.from(verts), indices: Uint32Array.from(tris) };
}

meshes.skin = marchingCubes(blur(occ, SIGMA), 0.5);
console.log('skin verts', meshes.skin.points.length / 3, 'tris', meshes.skin.indices.length / 3);

function minDist(a, b) {
  let best = 1e9;
  for (let i = 0; i < a.length; i += 3)
    for (let j = 0; j < b.length; j += 3) {
      const d = (a[i] - b[j]) ** 2 + (a[i + 1] - b[j + 1]) ** 2 + (a[i + 2] - b[j + 2]) ** 2;
      if (d < best) best = d;
    }
  return Math.sqrt(best);
}
const report = {
  sacculeToUtricleGapMm: minDist(meshes.saccule.points, meshes.utricle.points),
  ucpToUtricleGapMm: minDist(meshes.ucp.points, meshes.utricle.points),
  ucpToCrusGapMm: minDist(meshes.ucp.points, meshes['crus-p'].points),
  crusToUtricleGapMm: minDist(meshes['crus-p'].points, meshes.utricle.points),
};
console.log(report);

let maxAbs = 0;
for (const m of Object.values(meshes)) for (const v of m.points) maxAbs = Math.max(maxAbs, Math.abs(v));
const scale = maxAbs / 32000;
const b64 = (ta) => Buffer.from(ta.buffer, ta.byteOffset, ta.byteLength).toString('base64');
const pieces = [
  ...PIECES.map(([id, , label, group, inUnion]) => ({ id, label, group, inUnion })),
  { id: 'skin', label: 'Unified smooth skin (merged)', group: 'housing', inUnion: false },
];
for (const pc of pieces) {
  const m = meshes[pc.id];
  const q = new Int16Array(m.points.length);
  for (let i = 0; i < q.length; i++) q[i] = Math.round(m.points[i] / scale);
  const big = m.points.length / 3 > 65535;
  pc.vertexCount = m.points.length / 3;
  pc.pos = b64(q);
  pc.idx = b64(big ? m.indices : Uint16Array.from(m.indices));
  pc.idx32 = big;
}
writeFileSync(
  OUT,
  JSON.stringify({
    units: 'mm',
    frame: 'HeadFrame (+X anterior, +Y left, +Z superior), right ear, origin = utricle centroid',
    scale,
    report,
    pieces,
  })
);
console.log('wrote', OUT);
