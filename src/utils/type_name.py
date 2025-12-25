from __future__ import annotations

from pathlib import Path
from typing import Dict

try:
    import yaml  # type: ignore
except Exception:  # pragma: no cover
    yaml = None


def normalize_type_name(name: str | None) -> str:
    """
    基础规范化：
    - 转成字符串
    - 去首尾空格
    - 折叠中间多余空格
    不做大小写变化，避免和现有标签完全不匹配。
    """
    if name is None:
        return ""
    # 转成字符串并按空白切分再用单空格拼回
    parts = str(name).split()
    return " ".join(parts)


def load_alias_map(path: Path) -> Dict[str, str]:
    """
    从 YAML 加载别名映射:

    aliases:
      "B cell": "B cells"
      "CD4 T": "T cells CD4"

    key / value 都会先做 normalize，再存入映射。
    """
    alias_map: Dict[str, str] = {}
    if not path.exists() or yaml is None:
        return alias_map

    try:
        with path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
    except Exception:
        return alias_map

    raw_aliases = data.get("aliases", {}) if isinstance(data, dict) else {}
    if not isinstance(raw_aliases, dict):
        return alias_map

    for k, v in raw_aliases.items():
        nk = normalize_type_name(k)
        nv = normalize_type_name(v)
        if nk and nv:
            alias_map[nk] = nv

    return alias_map


def canonicalize_type_name(name: str | None, alias_map: Dict[str, str]) -> str:
    """
    先 normalize，再通过 alias_map 做一次映射。
    若无映射，则返回 normalize 后的结果。
    """
    norm = normalize_type_name(name)
    if not norm:
        return ""
    return alias_map.get(norm, norm)


