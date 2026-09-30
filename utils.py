import csv
import os
import subprocess
import sys
from pathlib import Path


def run_calc_stability(script_path: Path, pdb_path: Path, di_db: Path, penta_db: Path, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(script_path),
        "--pdb", str(pdb_path),
        "--di", str(di_db),
        "--penta", str(penta_db),
        "--out", str(out_dir),
    ]
    env = os.environ.copy()
    env.setdefault("PYTHONUNBUFFERED", "1")
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip() or "未知错误"
        raise RuntimeError(f"稳定性计算失败：{detail}")

    score_csv = out_dir / "residue_scores.csv"
    if not score_csv.exists():
        raise FileNotFoundError("稳定性计算完成，但没有生成 residue_scores.csv")
    return score_csv


def read_residue_scores(csv_path: Path):
    scores = {}
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            chain = (row.get("chain") or "").strip()
            resid = (row.get("resid") or "").strip()
            if not chain or not resid:
                continue
            scores[(chain, resid)] = float(row.get("score") or 0.0)
    return scores


def write_b_factor(input_pdb_path: Path, output_pdb_path: Path, res_score_map):
    """将每个残基的稳定性分数写入 PDB ATOM 记录的 B-factor 列。"""
    match_count = 0
    with input_pdb_path.open("r", encoding="utf-8", errors="ignore") as f_in, \
            output_pdb_path.open("w", encoding="utf-8") as f_out:
        for line in f_in:
            if line.startswith("ATOM"):
                chain_id = line[21:22].strip()
                res_seq = line[22:26].strip()
                icode = line[26:27].strip()
                resid = f"{res_seq}{icode}".strip()
                score = res_score_map.get((chain_id, resid), 0.0)
                if (chain_id, resid) in res_score_map:
                    match_count += 1

                # PDB B-factor 固定宽度 6 列（61-66）。
                new_b = f"{score:6.2f}"
                f_out.write(line[:60] + new_b + line[66:])
            else:
                f_out.write(line)
    return match_count
