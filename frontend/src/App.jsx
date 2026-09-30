import { useEffect, useMemo, useRef, useState } from 'react';
import './App.css';

const NGL_SOURCES = [
  'https://unpkg.com/ngl@2.0.0-dev.37/dist/ngl.js',
  'https://cdn.jsdelivr.net/npm/ngl@2.0.0-dev.37/dist/ngl.js',
];

let nglPromise;

function loadNGL() {
  if (window.NGL) return Promise.resolve(window.NGL);
  if (nglPromise) return nglPromise;

  nglPromise = new Promise((resolve, reject) => {
    let index = 0;

    const tryNext = () => {
      if (window.NGL) {
        resolve(window.NGL);
        return;
      }
      if (index >= NGL_SOURCES.length) {
        reject(new Error('3D 查看器加载失败，请检查网络后刷新页面。'));
        return;
      }

      const script = document.createElement('script');
      script.src = NGL_SOURCES[index++];
      script.async = true;
      script.onload = () => (window.NGL ? resolve(window.NGL) : tryNext());
      script.onerror = tryNext;
      document.head.appendChild(script);
    };

    tryNext();
  });

  return nglPromise;
}

function getResidueScores(pdbText) {
  const residues = new Map();

  for (const line of pdbText.split(/\r?\n/)) {
    if (!line.startsWith('ATOM')) continue;

    const chain = line.slice(21, 22).trim();
    const resid = line.slice(22, 26).trim() + line.slice(26, 27).trim();
    const key = `${chain}:${resid}`;
    const value = Number.parseFloat(line.slice(60, 66).trim());

    if (Number.isFinite(value)) residues.set(key, value);
  }

  return [...residues.values()];
}

function parseSecondaryStructure(pdbText) {
  const ranges = [];

  for (const line of pdbText.split(/\r?\n/)) {
    if (line.startsWith('HELIX')) {
      const chain = line.slice(19, 20).trim();
      const start = Number.parseInt(line.slice(21, 25).trim(), 10);
      const endChain = line.slice(31, 32).trim();
      const end = Number.parseInt(line.slice(33, 37).trim(), 10);
      if (Number.isFinite(start) && Number.isFinite(end)) {
        ranges.push({ type: 'H', chain, start, endChain: endChain || chain, end });
      }
    }

    if (line.startsWith('SHEET')) {
      const chain = line.slice(21, 22).trim();
      const start = Number.parseInt(line.slice(22, 26).trim(), 10);
      const endChain = line.slice(32, 33).trim();
      const end = Number.parseInt(line.slice(33, 37).trim(), 10);
      if (Number.isFinite(start) && Number.isFinite(end)) {
        ranges.push({ type: 'E', chain, start, endChain: endChain || chain, end });
      }
    }
  }

  return ranges;
}

function residueTypeFromRanges(atom, ranges) {
  const chain = String(atom.chainname || '').trim();
  const resno = Number(atom.resno);

  for (const range of ranges) {
    if (chain !== range.chain && chain !== range.endChain) continue;
    if (resno >= range.start && resno <= range.end) return range.type;
  }

  const fallback = String(atom.sstruc || '').toLowerCase();
  if (fallback === 'h') return 'H';
  if (fallback === 's') return 'E';
  return 'C';
}

function mix(a, b, t) {
  return Math.round(a + (b - a) * t);
}

function rgb(r, g, b) {
  return (r << 16) | (g << 8) | b;
}

function createScoreScheme(NGL, maxAbs, ranges) {
  const scale = Math.max(Number(maxAbs) || 1, 0.0001);

  return NGL.ColormakerRegistry.addScheme(function () {
    this.atomColor = function (atom) {
      const value = Number(atom.bfactor) || 0;
      const ss = residueTypeFromRanges(atom, ranges);

      // 没有参与评分的 Coil：中性灰色。
      if (ss === 'C') return 0x9fa6a2;

      // 接近 0：高饱和黄色。
      if (Math.abs(value) <= 0.05) return 0xffd400;

      // 负值：高饱和红色，数值越负越深。
      if (value < 0) {
        const t = Math.min(1, Math.abs(value) / scale);
        return rgb(mix(238, 155, t), mix(35, 0, t), mix(35, 18, t));
      }

      // 正值：高饱和绿色，数值越正越深。
      const t = Math.min(1, value / scale);
      return rgb(mix(20, 0, t), mix(190, 105, t), mix(70, 38, t));
    };
  });
}

function App() {
  const [file, setFile] = useState(null);
  const [loading, setLoading] = useState(false);
  const [message, setMessage] = useState('');
  const [error, setError] = useState('');
  const [pdbText, setPdbText] = useState('');
  const [downloadUrl, setDownloadUrl] = useState('');
  const [stats, setStats] = useState(null);
  const [viewerStatus, setViewerStatus] = useState('正在准备 3D 查看器…');
  const [viewerReady, setViewerReady] = useState(false);

  const viewerRef = useRef(null);
  const stageRef = useRef(null);
  const componentRef = useRef(null);

  const residueScores = useMemo(() => getResidueScores(pdbText), [pdbText]);
  const secondaryStructure = useMemo(() => parseSecondaryStructure(pdbText), [pdbText]);

  useEffect(() => {
    let disposed = false;

    loadNGL()
      .then((NGL) => {
        if (disposed || !viewerRef.current) return;
        const stage = new NGL.Stage(viewerRef.current, { backgroundColor: 'white' });
        stageRef.current = stage;
        setViewerReady(true);
        setViewerStatus('3D 查看器已就绪');
        window.addEventListener('resize', () => stage.handleResize());
      })
      .catch((err) => {
        console.error(err);
        if (!disposed) setViewerStatus(err.message || '3D 查看器加载失败');
      });

    return () => {
      disposed = true;
      if (componentRef.current) componentRef.current.dispose();
      componentRef.current = null;
      if (stageRef.current) stageRef.current.dispose();
      stageRef.current = null;
    };
  }, []);

  useEffect(() => {
    if (!pdbText || !viewerReady || !stageRef.current || !window.NGL) return;

    let cancelled = false;
    const stage = stageRef.current;
    const maxAbs = residueScores.reduce((m, value) => Math.max(m, Math.abs(value)), 0) || 1;

    async function loadStructure() {
      setViewerStatus('正在加载蛋白质 3D 结构…');
      try {
        if (componentRef.current) {
          componentRef.current.dispose();
          componentRef.current = null;
        }

        const blob = new Blob([pdbText], { type: 'text/plain' });
        const fileObject = new File([blob], 'alphabetteR-result.pdb', { type: 'text/plain' });
        const component = await stage.loadFile(fileObject, { ext: 'pdb' });
        if (cancelled) {
          component.dispose();
          return;
        }

        componentRef.current = component;
        const schemeId = createScoreScheme(window.NGL, maxAbs, secondaryStructure);

        component.addRepresentation('cartoon', {
          color: schemeId,
          quality: 'high',
          smoothSheet: true,
          aspectRatio: 4,
        });

        component.autoView(700);
        stage.handleResize();
        setViewerStatus('3D 结构加载完成');
      } catch (err) {
        console.error('3D 加载失败:', err);
        if (!cancelled) setViewerStatus(`3D 加载失败：${err.message || '未知错误'}`);
      }
    }

    loadStructure();

    return () => {
      cancelled = true;
    };
  }, [pdbText, residueScores, secondaryStructure, viewerReady]);

  function chooseFile(nextFile) {
    setError('');
    setMessage('');
    setStats(null);
    setPdbText('');
    setViewerStatus('正在准备 3D 查看器…');
    if (downloadUrl) URL.revokeObjectURL(downloadUrl);
    setDownloadUrl('');

    if (!nextFile) {
      setFile(null);
      return;
    }

    if (!nextFile.name.toLowerCase().endsWith('.pdb')) {
      setFile(null);
      setError('请选择 .pdb 格式的蛋白质结构文件。');
      return;
    }

    setFile(nextFile);
  }

  async function handleUpload() {
    if (!file || loading) return;

    setLoading(true);
    setError('');
    setMessage('正在上传并计算稳定性，请稍候…');
    setPdbText('');
    setStats(null);

    try {
      const formData = new FormData();
      formData.append('file', file);

      const response = await fetch('/api/upload', {
        method: 'POST',
        body: formData,
      });

      if (!response.ok) {
        let detail = `服务器返回 ${response.status}`;
        try {
          const data = await response.json();
          if (data.detail) detail = data.detail;
        } catch {
          // 保留默认错误信息。
        }
        throw new Error(detail);
      }

      const text = await response.text();
      const url = URL.createObjectURL(new Blob([text], { type: 'chemical/x-pdb' }));
      setDownloadUrl(url);
      setPdbText(text);

      const scores = getResidueScores(text);
      const positive = scores.filter((v) => v > 0.05).length;
      const negative = scores.filter((v) => v < -0.05).length;
      const neutral = scores.length - positive - negative;
      const nonZero = scores.filter((v) => Math.abs(v) > 0.05);

      setStats({
        residues: scores.length,
        positive,
        negative,
        neutral,
        min: nonZero.length ? Math.min(...nonZero) : 0,
        max: nonZero.length ? Math.max(...nonZero) : 0,
      });
      setMessage('分析完成！下面可以直接查看 3D 结构。');
    } catch (err) {
      console.error(err);
      setError(err.message || '分析失败，请稍后重试。');
      setMessage('');
    } finally {
      setLoading(false);
    }
  }

  function handleDownload() {
    if (!downloadUrl) return;
    const a = document.createElement('a');
    a.href = downloadUrl;
    a.download = `colored_${file?.name || 'protein.pdb'}`;
    document.body.appendChild(a);
    a.click();
    a.remove();
  }

  return (
    <main className="page-shell">
      <header className="hero">
        <div className="brand-mark">AB</div>
        <div>
          <div className="eyebrow">PROTEIN STABILITY ANALYSIS</div>
          <h1>AlphaBetter</h1>
          <p>蛋白质结构稳定性分析与 3D 可视化平台</p>
        </div>
      </header>

      <section className="upload-card">
        <div className="section-title">
          <span className="step">01</span>
          <div>
            <h2>上传 PDB 文件</h2>
            <p>选择蛋白质结构文件，系统将在云端自动完成稳定性分析。</p>
          </div>
        </div>

        <label className={`drop-zone ${file ? 'has-file' : ''}`} htmlFor="pdb-input">
          <input
            id="pdb-input"
            type="file"
            accept=".pdb"
            onChange={(event) => chooseFile(event.target.files?.[0] || null)}
          />
          <div className="upload-icon">↑</div>
          <strong>{file ? file.name : '点击选择 PDB 文件'}</strong>
          <span>{file ? '文件已准备好，可以开始分析' : '支持 .pdb 格式，单个文件不超过 20 MB'}</span>
        </label>

        <div className="actions">
          <button className="primary-button" onClick={handleUpload} disabled={!file || loading}>
            {loading ? '分析中…' : '开始分析'}
          </button>
          {downloadUrl && (
            <button className="secondary-button" onClick={handleDownload}>下载彩色 PDB</button>
          )}
        </div>

        {loading && <div className="progress-line" aria-label="正在分析"><span /></div>}
        {message && <div className="success-message">✓ {message}</div>}
        {error && <div className="error-message">{error}</div>}
      </section>

      {pdbText && (
        <section className="result-section">
          <div className="section-title result-heading">
            <span className="step">02</span>
            <div>
              <h2>稳定性结果</h2>
              <p>H/E 区域按稳定性评分着色：负值红色、接近 0 黄色、正值绿色；Coil 为灰色。</p>
            </div>
          </div>

          {stats && (
            <div className="stats-grid">
              <div className="stat-card"><span>残基数</span><strong>{stats.residues}</strong></div>
              <div className="stat-card green"><span>稳定区域</span><strong>{stats.positive}</strong></div>
              <div className="stat-card red"><span>不稳定区域</span><strong>{stats.negative}</strong></div>
              <div className="stat-card yellow"><span>中性区域</span><strong>{stats.neutral}</strong></div>
            </div>
          )}

          <div className="viewer-card">
            <div className="viewer-toolbar">
              <div className="legend">
                <span><i className="dot red-dot" />负值 / 不稳定</span>
                <span><i className="dot yellow-dot" />接近 0</span>
                <span><i className="dot green-dot" />正值 / 稳定</span>
                <span><i className="dot gray-dot" />Coil</span>
              </div>
              <span className="viewer-tip">拖动旋转 · 滚轮缩放 · 右键平移</span>
            </div>
            <div className="viewer-wrap">
              <div ref={viewerRef} className="viewer-canvas" />
              <div className="viewer-status">{viewerStatus}</div>
            </div>
          </div>

          <div className="score-note">
            <span>评分范围</span>
            <strong>{stats ? `${stats.min.toFixed(2)} ～ ${stats.max.toFixed(2)}` : '—'}</strong>
          </div>
        </section>
      )}

      <footer>AlphaBetter · Protein Stability Analysis & Visualization</footer>
    </main>
  );
}

export default App;
