from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract Xenium gene panel names to a text file")
    p.add_argument("--xenium_dir", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    xenium_dir = Path(args.xenium_dir).resolve()
    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    panel = json.loads((xenium_dir / "gene_panel.json").read_text(encoding="utf-8"))
    genes: list[str] = []
    for target in panel.get("payload", {}).get("targets", []):
        t = target.get("type", {})
        if t.get("descriptor") == "gene":
            name = t.get("data", {}).get("name")
            if name:
                genes.append(str(name))
    ordered: list[str] = []
    seen = set()
    for gene in genes:
        if gene not in seen:
            ordered.append(gene)
            seen.add(gene)
    out_path.write_text("\n".join(ordered) + "\n", encoding="utf-8")
    print(f"[panel] genes={len(ordered)} -> {out_path}")


if __name__ == "__main__":
    main()
