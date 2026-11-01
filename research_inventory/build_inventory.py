"""Rebuild the read-only experiment indexes from catalog.py and existing files.

Usage: python research_inventory/build_inventory.py
No solver is invoked. Only files inside research_inventory are written.
"""

from __future__ import annotations

import csv
import json
import subprocess
from collections import Counter
from pathlib import Path

from catalog import CASES, LEG


ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
SCAN = ("experiments", "examples", "outputs", "tests", "docs", "manuscript", "src")
DATA_EXT = {".csv", ".json", ".npz", ".log", ".txt"}
FIG_EXT = {".pdf", ".svg", ".png", ".jpg", ".jpeg"}
ORIG_FIG = Path(LEG) / "paper_draft_arxiv_cn/figures"


def paths(value: str) -> list[str]:
    return [p.strip() for p in value.split(";") if p.strip()]


def show(path: str) -> str:
    p = Path(path)
    if p.is_absolute():
        return f"`{p}`"
    return f"[`{p}`](../{p})"


def cell(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def all_paths() -> list[Path]:
    out = []
    for folder in SCAN:
        base = ROOT / folder
        if base.exists():
            out.extend(p for p in base.rglob("*") if p.is_file() and "__pycache__" not in p.parts)
    out.extend(p for name in ("README.md", "pyproject.toml") if (p := ROOT / name).is_file())
    return sorted(out)


def classify(path: Path) -> list[str]:
    s = str(path).lower()
    n = path.stem.lower()
    ids = []
    for c in CASES:
        for field in ("code", "data", "figure"):
            if any(str(path) == str((ROOT / p).resolve()) if not Path(p).is_absolute()
                   else str(path) == p for p in paths(c[field])):
                ids.append(c["id"])
                break
    if ids:
        return sorted(set(ids))
    if "engineering_inverse" in s:
        if "aircraft" in n:
            return ["P04"] if "tangent" in n else (["P03"] if any(t in n for t in ("screen", "shared", "reuse", "efficiency", "cost")) else ["I05"])
        if "automotive" in n:
            return ["P05"] if any(t in n for t in ("efficiency", "optimization_cost", "forward_cost", "convergence")) else ["I07"]
        return ["I05", "I07", "P03", "P04", "P05"]
    if "cmame_round6_applications" in s:
        return ["I04"] if "aircraft" in n else (["I06"] if "automotive" in n else ["I04", "I06"])
    if "cmame_round7_aircraft" in s:
        return ["I04"]
    if "cmame_round5" in s:
        if "shape" in n:
            return ["I02"]
        if "rotation" in n or "large" in n:
            return ["I03"]
        return ["I01"]
    if "cmame_positioning" in s:
        return ["M09"]
    if "cmame_relative_value" in s:
        if any(t in n for t in ("completion", "ablation", "spectrum")):
            return ["M04", "M06"]
        if any(t in n for t in ("cylinder", "isometry", "thickness", "locking")):
            return ["M05"]
        return ["P02"]
    if "plate" in n:
        return ["F01"]
    if "scordelis" in n or "roof" in n:
        return ["F02"]
    if "m3" in n or "hemisphere" in n:
        return ["F03", "F05"] if "nonlinear" in n else ["F03"]
    if "cantilever" in n:
        return ["F04"]
    if "cylinder" in n:
        return ["F06"]
    if "mesh" in n or "distortion" in n:
        return ["M08"]
    if "beta" in n:
        return ["M06"]
    if "spectrum" in n or "mode" in n:
        return ["M07"]
    if "gauge" in n:
        return ["M02"]
    if "geometry" in n or "curvature" in n:
        return ["M01"]
    return []


def master() -> None:
    labels = {"M": "几何与理论机制", "F": "正向壳力学", "I": "逆向设计", "P": "性能与实现", "X": "历史与独立支线"}
    lines = ["# 全项目独立数值实验总清单", "",
             "本清单以**研究问题/物理案例**而非网格、厚度、优化方法或轮次计数。每项下的系列见 [CASE_DETAILS.md](CASE_DETAILS.md)，逐系列机器索引见 [EXPERIMENT_MASTER.csv](EXPERIMENT_MASTER.csv)。",
             "", "检索范围：当前仓库 `experiments/`、`examples/`、`outputs/`、`tests/`、`docs/`、`manuscript/`；原始 arXiv 历史仓库 `/home/song/projects/CB-shell` 的 `v1.0-arxiv` 标签及输出；外部 4D 打印仓库仅作独立支线登记。当前仓库 Git 仅有 `0699638` 基线提交，绝大多数实验目录未被该 Git 追踪，不能从当前 Git 精确恢复每轮运行版本。", "",
             "**计数：32 项**＝机制 9＋正向 7＋逆向 7＋性能 5＋历史/外部 4。主仓库或原始论文相关内容为 28 项（含历史失败 F07）；X01–X04 是不同程度的独立支线，**不计作本篇壳单元论文的四个新验证**。其中非工程逆向为 I01–I03 共 3 项，工程/工程前驱为 I04–I07 共 4 项。", ""]
    for prefix, label in labels.items():
        lines += [f"## {label}", "", "| ID | 案例 | 来源 | 系列 | 可信状态 | 证据作用 |", "|---|---|---|---|---|---|"]
        for c in CASES:
            if c["id"].startswith(prefix):
                lines.append("| " + " | ".join(map(cell, (f"[{c['id']}](CASE_DETAILS.md#{c['id'].lower()})", c["name"], c["origin"], c["series"], c["status"], c["role"]))) + " |")
        lines.append("")
    OUT.joinpath("MASTER_CASE_INDEX.md").write_text("\n".join(lines), encoding="utf-8")


def details() -> None:
    lines = ["# 独立案例详细记录", "", "数值均引用已有记录；未重新运行求解器。`源程序`、`原始数据`、`代表图表`是已核实存在的**代表路径**，全部文件见 `ARTIFACT_MANIFEST.csv` 和 `FIGURE_AND_DATA_INDEX.md`。一项案例可以横跨多轮，系列名不另计案例。", ""]
    keys = [("原始研究目的", "purpose"), ("物理对象与设定", "setting"), ("实验系列", "series"), ("方法", "methods"), ("主要数值结果", "result"), ("参考与比较基准", "reference"), ("版本/轮次", "origin"), ("验证状态", "status"), ("已知问题", "issue"), ("能够支持", "supports"), ("不能支持", "cannot"), ("论文证据作用", "role")]
    for c in CASES:
        lines += [f"<a id=\"{c['id'].lower()}\"></a>", "", f"## {c['id']} {c['name']}", ""]
        lines.extend(f"- **{label}：** {c[key]}" for label, key in keys)
        lines += [f"- **源程序：** {'；'.join(map(show, paths(c['code']))) or '未找到'}",
                  f"- **原始数据：** {'；'.join(map(show, paths(c['data']))) or '未找到'}",
                  f"- **代表图表：** {'；'.join(map(show, paths(c['figure']))) or '未找到图表；请按数据/报告核实'}", ""]
    OUT.joinpath("CASE_DETAILS.md").write_text("\n".join(lines), encoding="utf-8")


def csv_index() -> None:
    columns = ["Case ID", "实验名称", "实验系列", "研究类别", "来源/版本", "物理对象与设定", "计算方法", "主要指标", "结果摘要", "比较基准", "代码路径", "数据路径", "对应图表", "可信状态", "已知问题", "可支持结论", "不能支持结论", "证据作用"]
    with OUT.joinpath("EXPERIMENT_MASTER.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for c in CASES:
            series = [s.strip() for s in c["series"].split("；") if s.strip()]
            for s in series:
                writer.writerow(dict(zip(columns, [c["id"], c["name"], s, c["category"], c["origin"], c["setting"], c["methods"], c["purpose"], c["result"], c["reference"], c["code"], c["data"], c["figure"], c["status"], c["issue"], c["supports"], c["cannot"], c["role"]])))


def artifacts() -> list[Path]:
    files = all_paths()
    with OUT.joinpath("ARTIFACT_MANIFEST.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["path", "extension", "size_bytes", "case_candidates", "attribution_status"])
        for p in files:
            candidates = classify(p)
            writer.writerow([str(p.relative_to(ROOT)), p.suffix.lower(), p.stat().st_size, ";".join(candidates), "exact/candidate by name" if candidates else "unassigned; see directory audit"])
    return files


def historical_artifacts() -> set[str]:
    hist = Path(LEG)
    tag = subprocess.check_output(
        ["git", "-C", str(hist), "ls-tree", "-r", "--name-only", "v1.0-arxiv"],
        text=True).splitlines()
    tagged = set(tag)
    files = []
    for folder in ("experiments", "outputs", "paper_draft_arxiv_cn/figures", "docs", "src"):
        base = hist / folder
        if base.exists():
            files.extend(p for p in base.rglob("*") if p.is_file() and "__pycache__" not in p.parts)
    with OUT.joinpath("HISTORICAL_ARTIFACT_MANIFEST.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["absolute_path", "relative_path", "size_bytes", "present_in_v1_arxiv_tag", "historical_scope"])
        for p in sorted(files):
            rel = str(p.relative_to(hist))
            scope = "original-paper figure/source" if rel.startswith("paper_draft_arxiv_cn/figures/") and rel in tagged else "historical worktree; tag status recorded separately"
            writer.writerow([str(p), rel, p.stat().st_size, int(rel in tagged), scope])
    return tagged


def historical_case(path: Path) -> str:
    s = str(path).lower()
    n = path.stem.lower()
    if "two_port" in s or "triangular_shell_validation" in s or "fpc" in s:
        return "X03"
    if "matrix_free" in s or "multiscenario_gpu" in s or "cross_method_gpu" in s:
        return "X02"
    if any(k in s for k in ("membrane_final", "locking_free", "intrinsic_m", "membrane_locking", "phase_completion")):
        return "X01"
    if "pinched_cylinder" in s:
        return "F07"
    if "open_cylinder" in s:
        return "F06"
    if "cantilever" in n:
        return "F04"
    if "nonlinear" in s and "hemisphere" in n:
        return "F05"
    if "plate" in n:
        return "F01"
    if "roof" in n or "scordelis" in n:
        return "F02"
    if "hemisphere" in n or "m3" in n:
        return "F03"
    if "geometry" in n or "bform" in s:
        return "M01"
    if "completion" in n or "blind" in n:
        return "M03;M04"
    return "历史待核实"


def figure_index(files: list[Path], tagged: set[str]) -> None:
    figures = [p for p in files if p.suffix.lower() in FIG_EXT]
    figures.extend(p for base in (Path(LEG) / "outputs", Path(LEG) / "experiments", ORIG_FIG) if base.exists()
                   for p in base.rglob("*") if p.is_file() and p.suffix.lower() in FIG_EXT)
    figures = sorted(set(figures))
    pyfiles = [p for p in files if p.suffix.lower() == ".py"]
    pytexts = {}
    for p in pyfiles:
        try:
            pytexts[p] = p.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            pass
    rows = []
    for fig in figures:
        stem = fig.stem
        local = fig.is_relative_to(ROOT)
        tagged_original = not local and str(fig.relative_to(LEG)) in tagged
        scope = "当前仓库" if local else ("arXiv v1 标签文件" if tagged_original else "历史仓库工作树；非标签文件")
        owner = ";".join(classify(fig)) if local else historical_case(fig)
        matches = [p for p, content in pytexts.items() if fig.name in content or stem in content]
        if local:
            same_root = fig.parent.parent if fig.parent.name in ("figures", "plots") else fig.parent
            data_candidates = [p for p in same_root.rglob("*") if p.is_file() and p.suffix.lower() in DATA_EXT and "__pycache__" not in p.parts]
            exact = [p for p in data_candidates if p.stem == stem]
            data = ";".join(str(p.relative_to(ROOT)) for p in exact[:10]) if exact else (str(same_root.relative_to(ROOT)) + "/*.{csv,json,npz,log}（目录级候选；未逐图验证）" if data_candidates else "未找到")
            script = ";".join(str(p.relative_to(ROOT)) for p in matches[:5]) if matches else "未定位精确生成脚本"
            fpath = str(fig.relative_to(ROOT))
        else:
            data = str(fig.parent.parent) + "/*.{csv,json,npz}（历史目录级候选；未逐图验证）"
            script = "未在当前仓库定位精确生成脚本；见历史仓库"
            fpath = str(fig)
        paper_status = ("原稿图目录且标签存在" if tagged_original and fig.is_relative_to(ORIG_FIG)
                        else "标签存在；原稿是否引用待核实" if tagged_original
                        else "后续或历史版本" if not local
                        else "后续新增/待与原稿核对")
        rows.append((fpath, scope, owner, script, data, paper_status))
    with OUT.joinpath("FIGURE_AND_DATA_INDEX.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["figure_path", "source", "case_candidates", "generating_script_evidence", "raw_data_evidence", "paper_status"])
        writer.writerows(rows)
    counts = Counter(r[1] for r in rows)
    lines = ["# 现有图表与原始数据索引", "",
             f"已逐文件枚举 {len(rows)} 个 PDF/SVG/PNG/JPEG 图表文件：" + "，".join(f"{k} {v}" for k, v in counts.items()) + "。同一图的多种格式分别保留，不能把文件数当作独立实验数。逐行机器索引见 [FIGURE_AND_DATA_INDEX.csv](FIGURE_AND_DATA_INDEX.csv)。",
             "", "**归属规则：** 代表图在 `CASE_DETAILS.md` 中经存在性核查。下表其他文件按名称/目录给出候选案例；`未定位精确生成脚本` 与 `目录级候选` 不表示已证实逐图数据链。原稿发布图与后续翻绘图可能同名，应以路径和版本判别。", "",
             "| 图表文件 | 来源 | 案例候选 | 生成脚本证据 | 原始数据证据 |", "|---|---|---|---|---|"]
    for fpath, scope, owner, script, data, _ in rows:
        lines.append("| " + " | ".join(cell(x) for x in (fpath, scope, owner or "待核实", script, data)) + " |")
    OUT.joinpath("FIGURE_AND_DATA_INDEX.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def audit() -> None:
    missing = []
    for c in CASES:
        for field in ("code", "data", "figure"):
            for p in paths(c[field]):
                target = Path(p) if Path(p).is_absolute() else ROOT / p
                if not target.exists():
                    missing.append({"case": c["id"], "field": field, "path": p})
    observed = sorted(p.name for p in (ROOT / "experiments").iterdir() if p.is_dir())
    coverage = {d: [c["id"] for c in CASES if f"experiments/{d}/" in " ".join((c["code"], c["data"], c["figure"]))] for d in observed}
    report = {"case_count": len(CASES), "by_prefix": dict(Counter(c["id"][0] for c in CASES)), "experiment_directory_coverage": coverage, "missing_representative_paths": missing}
    OUT.joinpath("COVERAGE_AUDIT.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    master()
    details()
    csv_index()
    scanned = artifacts()
    tagged = historical_artifacts()
    figure_index(scanned, tagged)
    audit()
    print(f"indexed {len(CASES)} cases, {len(scanned)} current-repo artifacts")
