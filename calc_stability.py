import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np


AA3_TO_AA1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}
STANDARD_AA3 = set(AA3_TO_AA1.keys())


def aa3_to_aa1(resname):
    return AA3_TO_AA1.get((resname or "").upper(), "X")


def normalize_ss(ss):
    ss = (ss or "").strip().upper()
    if ss in {"H", "G", "I"}:
        return ss
    if ss in {"E", "B"}:
        return ss
    if ss in {"S", "T"}:
        return ss
    return "C"


def load_database(csv_path, expected_length):
    """
    读取数据库文件。
    支持列：
        peptide, ss, db_value
    额外列会自动忽略。
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"数据库文件不存在: {csv_path}")

    db = {}
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        required = {"peptide", "ss", "db_value"}
        if not reader.fieldnames:
            raise ValueError(f"{csv_path} 为空或缺少表头")

        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{csv_path} 缺少列: {sorted(missing)}")

        for row in reader:
            peptide = (row.get("peptide") or "").strip().upper()
            ss = normalize_ss(row.get("ss"))

            if len(peptide) != expected_length:
                continue
            if ss not in {"H", "E"}:
                continue

            try:
                value = float(row.get("db_value"))
            except Exception:
                continue

            if not math.isfinite(value):
                continue

            db[(peptide, ss)] = value

    if not db:
        raise ValueError(f"{csv_path} 中未读取到有效数据库记录")

    return db


def parse_pdb_atom_line(line):
    """
    解析 PDB ATOM/HETATM 行（按固定列）
    返回：
      chain_id, resname, resseq, icode, atom_name, coord(x,y,z)
    """
    atom_name = line[12:16].strip()
    resname = line[17:20].strip().upper()
    chain_id = line[21].strip()
    resseq = line[22:26].strip()
    icode = line[26].strip()

    x = line[30:38].strip()
    y = line[38:46].strip()
    z = line[46:54].strip()

    if not resseq:
        return None
    try:
        resseq = int(resseq)
        coord = np.array([float(x), float(y), float(z)], dtype=float)
    except Exception:
        return None

    return chain_id, resname, resseq, icode, atom_name, coord


def signed_dihedral(p1, p2, p3, p4):
    p1 = np.array(p1, dtype=float)
    p2 = np.array(p2, dtype=float)
    p3 = np.array(p3, dtype=float)
    p4 = np.array(p4, dtype=float)

    b0 = -1.0 * (p2 - p1)
    b1 = p3 - p2
    b2 = p4 - p3

    b1 /= (np.linalg.norm(b1) + 1e-12)

    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def angle_deg(v1, v2):
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < 1e-12 or n2 < 1e-12:
        return None
    cosv = float(np.dot(v1, v2) / (n1 * n2))
    cosv = max(-1.0, min(1.0, cosv))
    return float(np.degrees(np.arccos(cosv)))


def vec_norm(v):
    v = np.array(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-12:
        return v * 0.0
    return v / n


@dataclass
class ResidueRecord:
    chain: str
    resnum: int
    icode: str
    aa: str
    resname: str
    atoms: dict
    ss: str = "C"


def extract_first_chain_records(pdb_path):
    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB 文件不存在: {pdb_path}")

    first_chain = None
    with pdb_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            parsed = parse_pdb_atom_line(line)
            if parsed is None:
                continue
            chain_id = parsed[0]
            if chain_id:
                first_chain = chain_id
                break

    if first_chain is None:
        raise ValueError("PDB 中未找到有效的 ATOM/HETATM 记录或链 ID")

    residues_map = {}
    order = []

    with pdb_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            parsed = parse_pdb_atom_line(line)
            if parsed is None:
                continue

            chain_id, resname, resseq, icode, atom_name, coord = parsed
            if chain_id != first_chain:
                continue

            if resname not in STANDARD_AA3:
                continue

            key = (chain_id, resseq, icode, resname)
            if key not in residues_map:
                residues_map[key] = ResidueRecord(
                    chain=chain_id,
                    resnum=resseq,
                    icode=icode,
                    aa=aa3_to_aa1(resname),
                    resname=resname,
                    atoms={},
                )
                order.append(key)

            residues_map[key].atoms[atom_name] = coord

    residues = [residues_map[k] for k in order]
    if not residues:
        raise ValueError("第一条链中未解析到标准氨基酸残基")

    return residues


def residue_key(res):
    return (res.chain, res.resnum, res.icode)


def are_consecutive(r1, r2):
    return (
        r2.resnum == r1.resnum + 1
        and (r2.icode or "") == (r1.icode or "")
    )


def get_atom(res, name):
    return res.atoms.get(name)


def has_backbone_atoms(res):
    return all(a in res.atoms for a in ("N", "CA", "C", "O"))


def build_peptide_backbone(residues):
    """
    为缺氢残基补 N-H。
    返回每个 residue 的氢键供体位置字典:
      donor_h[res_idx] = np.array([x,y,z]) or None
    """
    donor_h = [None] * len(residues)

    for i, res in enumerate(residues):
        if res.resname == "PRO":
            continue
        if "H" in res.atoms or "HN" in res.atoms or "1H" in res.atoms:
            for hname in ("H", "HN", "1H"):
                if hname in res.atoms:
                    donor_h[i] = res.atoms[hname]
                    break
            continue

        if not has_backbone_atoms(res):
            continue

        n = res.atoms["N"]
        ca = res.atoms["CA"]
        c_prev = None
        if i > 0 and has_backbone_atoms(residues[i - 1]):
            c_prev = residues[i - 1].atoms["C"]
        elif "C" in res.atoms:
            c_prev = res.atoms["C"]

        if c_prev is None:
            continue

        v1 = vec_norm(n - c_prev)
        v2 = vec_norm(ca - n)
        bis = vec_norm(v1 + v2)
        if np.linalg.norm(bis) < 1e-8:
            bis = vec_norm(ca - c_prev)

        h_pos = n + 1.01 * bis
        donor_h[i] = h_pos

    return donor_h


def compute_hydrogen_bonds(residues, donor_h, cutoff=-0.5):
    """
    计算氢键：E = 27.888*(1/ON + 1/CH - 1/OH - 1/CN)
    返回:
      hbonds = list of dict
    """
    hbonds = []

    for i, donor in enumerate(residues):
        if donor_h[i] is None:
            continue
        if not has_backbone_atoms(donor):
            continue

        n = donor.atoms["N"]
        h = donor_h[i]
        c = donor.atoms["C"]

        for j, acceptor in enumerate(residues):
            if i == j:
                continue
            if "O" not in acceptor.atoms or "C" not in acceptor.atoms:
                continue

            o = acceptor.atoms["O"]
            c_acc = acceptor.atoms["C"]

            ON = np.linalg.norm(o - n)
            CH = np.linalg.norm(c - h)
            OH = np.linalg.norm(o - h)
            CN = np.linalg.norm(c_acc - n)

            if min(ON, CH, OH, CN) < 1e-6:
                continue

            e = 27.888 * (1.0 / ON + 1.0 / CH - 1.0 / OH - 1.0 / CN)
            if e < cutoff:
                hbonds.append({
                    "donor": i,
                    "acceptor": j,
                    "energy": float(e),
                    "ON": float(ON),
                    "CH": float(CH),
                    "OH": float(OH),
                    "CN": float(CN),
                })

    return hbonds


def assign_helices(residues, hbonds):
    """
    识别 H/G/I:
      alpha: i,i+4
      3/10: i,i+3
      pi: i,i+5
    采用简单稳定规则：连续模式成段后标注。
    """
    n = len(residues)
    ss = ["C"] * n

    bond_pairs = {(hb["donor"], hb["acceptor"]) for hb in hbonds}

    def mark_pattern(offset, label, min_run=2):
        runs = []
        i = 0
        while i + offset < n:
            if (i, i + offset) in bond_pairs:
                start = i
                end = i + offset
                run_pairs = 1
                k = i + 1
                while k + offset < n and (k, k + offset) in bond_pairs:
                    run_pairs += 1
                    end = k + offset
                    k += 1
                if run_pairs >= min_run:
                    runs.append((start, end, run_pairs))
                i = k
            else:
                i += 1

        for start, end, _ in runs:
            for idx in range(start, end + 1):
                if label == "H":
                    ss[idx] = "H"
                elif label == "G":
                    if ss[idx] == "C":
                        ss[idx] = "G"
                elif label == "I":
                    if ss[idx] == "C":
                        ss[idx] = "I"

    mark_pattern(4, "H", min_run=2)
    mark_pattern(3, "G", min_run=2)
    mark_pattern(5, "I", min_run=1)

    return ss


def ca_bend_angle(residues, center_idx):
    if center_idx - 3 < 0 or center_idx + 3 >= len(residues):
        return None
    r1 = residues[center_idx - 3]
    r2 = residues[center_idx]
    r3 = residues[center_idx + 3]
    if "CA" not in r1.atoms or "CA" not in r2.atoms or "CA" not in r3.atoms:
        return None
    return angle_deg(r1.atoms["CA"] - r2.atoms["CA"], r3.atoms["CA"] - r2.atoms["CA"])


def assign_bends(residues, ss, bend_threshold=75.0):
    """
    按报告：内部夹角 < 75° 判为 bend (S)
    优先级高于所有元素，但如果落在螺旋第2残基，螺旋自此开始。
    """
    n = len(residues)
    bend_flags = [False] * n

    for i in range(n):
        ang = ca_bend_angle(residues, i)
        if ang is None:
            continue
        if ang < bend_threshold:
            bend_flags[i] = True

    for i in range(n):
        if bend_flags[i]:
            if ss[i] in {"H", "G", "I"}:
                if i > 0 and ss[i - 1] in {"H", "G", "I"}:
                    continue
            ss[i] = "S"

    return ss


def assign_beta_bridges(residues, hbonds):
    """
    简化版 beta bridge / ladder:
    - 使用主链氢键 i->j
    - 若形成成对/连续配对，则为 ladder
    - ladder >= 2 桥 => E
    - 孤立桥 => B
    """
    n = len(residues)
    ss = ["C"] * n
    bridge_pairs = []

    for hb in hbonds:
        i = hb["donor"]
        j = hb["acceptor"]
        if abs(i - j) < 2:
            continue
        bridge_pairs.append(tuple(sorted((i, j))))

    bridge_pairs = sorted(set(bridge_pairs))
    adj = defaultdict(set)
    for i, j in bridge_pairs:
        adj[i].add(j)
        adj[j].add(i)

    visited = set()
    for i, j in bridge_pairs:
        if (i, j) in visited:
            continue
        ladder = [(i, j)]
        visited.add((i, j))

        cur_i, cur_j = i, j
        while True:
            next_pair = None
            for ni in (cur_i + 1, cur_i - 1):
                for nj in (cur_j + 1, cur_j - 1):
                    if (min(ni, nj), max(ni, nj)) in bridge_pairs:
                        next_pair = (min(ni, nj), max(ni, nj))
                        break
                if next_pair:
                    break
            if next_pair and next_pair not in visited:
                ladder.append(next_pair)
                visited.add(next_pair)
                cur_i, cur_j = next_pair
            else:
                break

        if len(ladder) >= 2:
            for a, b in ladder:
                ss[a] = "E"
                ss[b] = "E"
        else:
            a, b = ladder[0]
            if ss[a] == "C":
                ss[a] = "B"
            if ss[b] == "C":
                ss[b] = "B"

    return ss


def combine_ss(residues, ss_helices, ss_beta, ss_bend):
    """
    按优先级合并：
    S > G > H > I > E > B > C
    这里为了最终输出兼容打分，映射为 H/E/C,
    但保留中间状态用于调试。
    """
    n = len(residues)
    final = ["C"] * n

    for i in range(n):
        candidates = [ss_beta[i], ss_helices[i], ss_bend[i]]

        if ss_bend[i] == "S":
            final[i] = "C"  # bend 不直接用于打分，保留为 C
        elif ss_helices[i] in {"G", "H", "I"}:
            final[i] = ss_helices[i]
        elif ss_beta[i] in {"E", "B"}:
            final[i] = ss_beta[i]
        else:
            final[i] = "C"

    return final


def run_ussa_from_residues(residues):
    """
    USSA 风格结构分配：
    1) 补氢
    2) 氢键
    3) 螺旋
    4) beta bridge / ladder
    5) bend
    最终输出 H/E/C
    """
    donor_h = build_peptide_backbone(residues)
    hbonds = compute_hydrogen_bonds(residues, donor_h, cutoff=-0.5)

    ss_helices = assign_helices(residues, hbonds)
    ss_beta = assign_beta_bridges(residues, hbonds)
    ss_bend = ["C"] * len(residues)
    ss_bend = assign_bends(residues, ss_helices[:], bend_threshold=75.0)

    final = combine_ss(residues, ss_helices, ss_beta, ss_bend)

    # 最终映射为 H/E/C
    out = []
    for x in final:
        if x in {"H", "G", "I"}:
            out.append("H")
        elif x in {"E", "B"}:
            out.append("E")
        else:
            out.append("C")
    return out


def parse_pdb_secondary_structure(pdb_path, residues):
    """
    读取 PDB 文件自带的 HELIX / SHEET 注释，并映射到已解析的残基。
    对本项目的稳定性数据库只需要 H / E / C 三类：
      HELIX -> H
      SHEET -> E
      其它 -> C

    注意：PDB 的残基编号可能带 insertion code，因此同时比较
    (chain, resseq, icode)。
    """
    annotations = []

    with Path(pdb_path).open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            record = line[0:6].strip()

            if record == "HELIX":
                try:
                    chain = line[19].strip()
                    start_res = int(line[21:25].strip())
                    start_icode = line[25].strip()
                    end_res = int(line[33:37].strip())
                    end_icode = line[37].strip()
                    annotations.append(("H", chain, start_res, start_icode,
                                        end_res, end_icode))
                except (ValueError, IndexError):
                    continue

            elif record == "SHEET":
                try:
                    chain = line[21].strip()
                    start_res = int(line[22:26].strip())
                    start_icode = line[26].strip()
                    end_res = int(line[33:37].strip())
                    end_icode = line[37].strip()
                    annotations.append(("E", chain, start_res, start_icode,
                                        end_res, end_icode))
                except (ValueError, IndexError):
                    continue

    if not annotations:
        return None

    def in_range(res, start_res, start_icode, end_res, end_icode):
        # 当前 test PDB 没有复杂 insertion code；这里仍保留基本支持。
        start_key = (start_res, start_icode or "")
        end_key = (end_res, end_icode or "")
        res_key = (res.resnum, res.icode or "")
        return start_key <= res_key <= end_key

    ss = ["C"] * len(residues)
    assigned = 0

    # SHEET/HELIX 可能存在重叠。优先保留先出现的明确注释，避免
    # 后面的记录把已有结构覆盖掉。
    for i, res in enumerate(residues):
        for label, chain, start_res, start_icode, end_res, end_icode in annotations:
            if res.chain != chain:
                continue
            if in_range(res, start_res, start_icode, end_res, end_icode):
                ss[i] = label
                assigned += 1
                break

    return ss if assigned else None


def extract_first_chain_records_with_ss(pdb_path):
    residues = extract_first_chain_records(pdb_path)

    # 先保留原有的几何/氢键方式计算；如果它在当前 PDB 上无法识别出
    # H/E，则使用 PDB 自带的 HELIX/SHEET 注释作为可靠的结构标签。
    # 这样不会把“结构识别失败”误当成“所有残基都是 C”。
    ss_list = run_ussa_from_residues(residues)
    detected_he = sum(ss in {"H", "E"} for ss in ss_list)

    pdb_ss = parse_pdb_secondary_structure(pdb_path, residues)
    if pdb_ss is not None and (detected_he == 0 or detected_he < max(2, len(residues) // 20)):
        ss_list = pdb_ss
        print(f"[INFO] USSA 未可靠识别二级结构（H/E={detected_he}），改用 PDB HELIX/SHEET 注释。")

    for res, ss in zip(residues, ss_list):
        res.ss = ss
    return residues


def are_consecutive_rec(r1, r2):
    return (
        r2.resnum == r1.resnum + 1
        and (r2.icode or "") == (r1.icode or "")
    )


def contiguous(window):
    for i in range(1, len(window)):
        if not are_consecutive_rec(window[i - 1], window[i]):
            return False
    return True


def calc_stability(pdb_file, di_csv, penta_csv, out_dir):
    di_db = load_database(di_csv, expected_length=2)
    penta_db = load_database(penta_csv, expected_length=5)

    records = extract_first_chain_records_with_ss(pdb_file)
    chain_id = records[0].chain

    n = len(records)
    residue_scores = [0.0] * n
    hits = []

    # 二肽打分
    for i in range(n - 1):
        w = records[i:i + 2]
        if not contiguous(w):
            continue

        ss0 = normalize_ss(w[0].ss)
        ss1 = normalize_ss(w[1].ss)

        if ss0 != ss1 or ss0 not in {"H", "E"}:
            continue

        pep = w[0].aa + w[1].aa
        val = di_db.get((pep, ss0))
        if val is not None:
            residue_scores[i] += val
            hits.append({
                "chain": chain_id,
                "type": "di",
                "start": i + 1,
                "end": i + 2,
                "center": "",
                "peptide": pep,
                "ss": ss0,
                "db_value": val,
            })

    # 五肽打分
    for center in range(2, n - 2):
        w = records[center - 2:center + 3]
        if not contiguous(w):
            continue

        ss_list = [normalize_ss(x.ss) for x in w]
        if len(set(ss_list)) != 1 or ss_list[0] not in {"H", "E"}:
            continue

        ss = ss_list[0]
        pep = "".join(x.aa for x in w)
        val = penta_db.get((pep, ss))
        if val is not None:
            residue_scores[center] += val
            hits.append({
                "chain": chain_id,
                "type": "penta",
                "start": center - 1,
                "end": center + 3,
                "center": center + 1,
                "peptide": pep,
                "ss": ss,
                "db_value": val,
            })

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    residue_csv = out_dir / "residue_scores.csv"
    hits_csv = out_dir / "hits.csv"

    with residue_csv.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["chain", "resid", "score"])
        writer.writeheader()
        for rec, score in zip(records, residue_scores):
            resid = f"{rec.resnum}{rec.icode}".strip()
            writer.writerow({
                "chain": rec.chain,
                "resid": resid,
                "score": score,
            })

    with hits_csv.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["chain", "type", "start", "end", "center", "peptide", "ss", "db_value"]
        )
        writer.writeheader()
        writer.writerows(hits)

    print(f"[DONE] residue_scores: {residue_csv}")
    print(f"[DONE] hits: {hits_csv}")


def main():
    parser = argparse.ArgumentParser(description="calc_stability")
    parser.add_argument("--pdb", required=True, help="原始 PDB 文件")
    parser.add_argument("--di", required=True, help="二肽数据库 CSV")
    parser.add_argument("--penta", required=True, help="五肽数据库 CSV")
    parser.add_argument("--out", required=True, help="输出目录")
    args = parser.parse_args()

    calc_stability(
        pdb_file=args.pdb,
        di_csv=args.di,
        penta_csv=args.penta,
        out_dir=args.out,
    )


if __name__ == "__main__":
    main()