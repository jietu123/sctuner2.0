from __future__ import annotations

from pathlib import Path
from typing import Iterable


def sample_dir_candidates(project_root: Path, sample: str, sim_group: str = "real_brca") -> list[Path]:
    root = Path(project_root).resolve()
    s = str(sample).strip()
    cands: list[Path] = [root / "data" / "sim" / sim_group / s]
    sim_root = root / "data" / "sim"
    if sim_root.exists():
        for grp in sorted(sim_root.iterdir()):
            if not grp.is_dir():
                continue
            cand = grp / s
            if cand not in cands:
                cands.append(cand)
    for extra in [root / "data" / "sim" / s, root / "data" / "raw" / s]:
        if extra not in cands:
            cands.append(extra)
    return cands


def resolve_sample_dir(
    project_root: Path,
    sample: str,
    *,
    sim_group: str = "real_brca",
    must_exist: bool = True,
) -> Path:
    cands = sample_dir_candidates(project_root, sample, sim_group=sim_group)
    for c in cands:
        if c.exists():
            return c
    if must_exist:
        raise FileNotFoundError(
            f"sample dir not found for '{sample}'. tried: "
            + ", ".join(str(x) for x in cands)
        )
    return cands[0]


def resolve_from_roots(roots: Iterable[Path], rel_path: str) -> Path | None:
    p = Path(rel_path)
    for r in roots:
        cand = Path(r) / p
        if cand.exists():
            return cand
    return None


def infer_sim_group(
    project_root: Path,
    sample: str,
    *,
    sim_group: str = "real_brca",
) -> str | None:
    """
    Infer sample group under data/sim/<group>/<sample>.
    Returns None when sample resolves outside data/sim/<group>/ (e.g. data/raw).
    """
    root = Path(project_root).resolve()
    s = str(sample).strip()
    sim_root = root / "data" / "sim"
    if sim_root.exists():
        for grp in sorted(sim_root.iterdir()):
            if not grp.is_dir():
                continue
            if (grp / s).exists():
                return str(grp.name)

    sample_dir = resolve_sample_dir(root, s, sim_group=sim_group, must_exist=True)
    try:
        rel = sample_dir.resolve().relative_to(sim_root.resolve())
    except Exception:
        return None
    parts = rel.parts
    if len(parts) >= 2:
        return str(parts[0])

    name = s.lower()
    if name.startswith("adult_mouse_kidney") or name.startswith("amk_"):
        return "adult_mouse_kidney"
    if name.startswith("real_brca"):
        return "real_brca"
    return None
