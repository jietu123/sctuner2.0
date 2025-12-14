from __future__ import annotations

import ast
import re
import sys
from pathlib import Path


def _die(msg: str) -> None:
    raise SystemExit(msg)


def _count_defs(src: str, name: str) -> int:
    return len(re.findall(rf"^def {re.escape(name)}\b", src, flags=re.M))


def _get_run_plus_chunk(src: str) -> str:
    tree = ast.parse(src)
    run_plus = None
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == "CytoSPACEBackend":
            for m in node.body:
                if isinstance(m, ast.FunctionDef) and m.name == "run_plus":
                    run_plus = m
                    break
    if run_plus is None:
        _die("run_plus not found in CytoSPACEBackend")
    if not hasattr(run_plus, "end_lineno") or run_plus.end_lineno is None:
        _die("Python AST missing end_lineno; please run with Python>=3.8")
    lines = src.splitlines(True)
    return "".join(lines[run_plus.lineno - 1 : run_plus.end_lineno])


def _get_method_chunk(src: str, class_name: str, method_name: str) -> str:
    tree = ast.parse(src)
    target = None
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for m in node.body:
                if isinstance(m, ast.FunctionDef) and m.name == method_name:
                    target = m
                    break
    if target is None:
        _die(f"{method_name} not found in {class_name}")
    if not hasattr(target, "end_lineno") or target.end_lineno is None:
        _die("Python AST missing end_lineno; please run with Python>=3.8")
    lines = src.splitlines(True)
    return "".join(lines[target.lineno - 1 : target.end_lineno])


def _method_calls_func(method: ast.FunctionDef, func_name: str) -> bool:
    for n in ast.walk(method):
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == func_name:
            return True
    return False


def _find_deprecated_key_usages(src: str, keys: set[str]) -> list[tuple[int, str]]:
    tree = ast.parse(src)
    bad: list[tuple[int, str]] = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == "get":
            if n.args and isinstance(n.args[0], ast.Constant) and isinstance(n.args[0].value, str):
                k = n.args[0].value
                if k in keys:
                    bad.append((n.lineno or -1, k))
        if isinstance(n, ast.Subscript) and isinstance(n.value, ast.Name):
            sl = n.slice
            if isinstance(sl, ast.Constant) and isinstance(sl.value, str) and sl.value in keys:
                bad.append((n.lineno or -1, sl.value))
    return bad


def _find_bad_spot_fraction_calls(src: str) -> int:
    tree = ast.parse(src)
    bad = 0
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Name) or node.func.id != "_spot_fraction_from_mat":
            continue
        if len(node.args) < 2:
            continue
        second = node.args[1]
        if not isinstance(second, ast.Call):
            continue
        if not isinstance(second.func, ast.Name) or second.func.id != "list":
            continue
        if len(second.args) != 1:
            continue
        arg0 = second.args[0]
        if isinstance(arg0, ast.Attribute) and isinstance(arg0.value, ast.Name):
            if arg0.value.id == "sc_expr" and arg0.attr == "index":
                bad += 1
    return bad


def _get_class_method(src: str, class_name: str, method_name: str) -> ast.FunctionDef:
    tree = ast.parse(src)
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for m in node.body:
                if isinstance(m, ast.FunctionDef) and m.name == method_name:
                    return m
    _die(f"{method_name} not found in {class_name}")


def _extract_meta_dict(method: ast.FunctionDef, *, method_name: str) -> ast.Dict:
    meta_nodes: list[ast.Dict] = []
    for n in ast.walk(method):
        if not isinstance(n, ast.Assign):
            continue
        for t in n.targets:
            if isinstance(t, ast.Name) and t.id == "meta" and isinstance(n.value, ast.Dict):
                meta_nodes.append(n.value)
    if not meta_nodes:
        _die(f"meta dict assignment not found in {method_name}")
    return meta_nodes[-1]


def _dict_str_keys(d: ast.Dict) -> dict[str, ast.AST]:
    out: dict[str, ast.AST] = {}
    for k, v in zip(d.keys, d.values):
        if isinstance(k, ast.Constant) and isinstance(k.value, str):
            out[k.value] = v
    return out


def _unused_private_funcs(src: str) -> list[str]:
    tree = ast.parse(src)
    funcs = [
        n.name
        for n in tree.body
        if isinstance(n, ast.FunctionDef) and n.name.startswith("_") and not n.name.startswith("__")
    ]
    refs = {
        n.id
        for n in ast.walk(tree)
        if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load) and n.id in funcs
    }
    unused = sorted(set(funcs) - refs)
    allow_unused: set[str] = set()
    return [f for f in unused if f not in allow_unused]


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    target = repo_root / "src" / "stages" / "backends" / "cytospace_backend.py"
    src = target.read_text(encoding="utf-8")

    # C1) 禁止重复定义/多入口
    def_counts = {
        "_spot_fraction_from_mat": _count_defs(src, "_spot_fraction_from_mat"),
        "_parse_assigned_locations": _count_defs(src, "_parse_assigned_locations"),
        "_harden_assignment_quota_matching": _count_defs(src, "_harden_assignment_quota_matching"),
        "_validate_and_resolve_config": _count_defs(src, "_validate_and_resolve_config"),
    }
    for k, v in def_counts.items():
        if v != 1:
            _die(f"FAIL: def {k} count != 1 (got {v})")

    run_plus_chunk = _get_run_plus_chunk(src)
    run_baseline_chunk = _get_method_chunk(src, "CytoSPACEBackend", "run_baseline")
    harden_calls = len(re.findall(r"_harden_assignment_quota_matching\(", run_plus_chunk))
    build_outputs_calls = len(re.findall(r"_build_outputs\(", run_plus_chunk))
    if harden_calls not in (1, 2):
        _die(f"FAIL: run_plus harden calls should be 1 or 2 (got {harden_calls})")
    if build_outputs_calls != 1:
        _die(f"FAIL: run_plus _build_outputs calls should be 1 (got {build_outputs_calls})")

    # C2) 禁止旧模式残留（字符串规则）
    if "header=False" in src:
        _die("FAIL: found forbidden 'header=False' in cytospace_backend.py")
    if "else 1.0" in src:
        _die("FAIL: found forbidden fallback 'else 1.0' (ablation change_rate must not fake 1.0)")

    bad_spot_fraction = _find_bad_spot_fraction_calls(src)
    if bad_spot_fraction > 0:
        _die(f"FAIL: found _spot_fraction_from_mat calls with list(sc_expr.index) as cell_types: {bad_spot_fraction}")

    # C2.5) Risk#7: 废弃 key 不允许在代码中被读取（允许出现在 deprecation map 里）
    deprecated_keys = {"capacity_normalize_mode", "capacity_normalize_factor", "capacity_norm_mode"}
    deprecated_usages = _find_deprecated_key_usages(src, deprecated_keys)
    if deprecated_usages:
        _die(f"FAIL: deprecated config keys used in code: {deprecated_usages[:5]}")

    # C2.6) Risk#7: run_baseline/run_plus 必须调用统一配置校验入口
    run_plus_node = _get_class_method(src, "CytoSPACEBackend", "run_plus")
    run_baseline_node = _get_class_method(src, "CytoSPACEBackend", "run_baseline")
    if not _method_calls_func(run_plus_node, "_validate_and_resolve_config"):
        _die("FAIL: run_plus must call _validate_and_resolve_config(...)")
    if not _method_calls_func(run_baseline_node, "_validate_and_resolve_config"):
        _die("FAIL: run_baseline must call _validate_and_resolve_config(...)")

    # C3) 未使用函数检查
    unused = _unused_private_funcs(src)
    if unused:
        _die(f"FAIL: unused private funcs detected: {unused}")

    # C4) Risk#6：meta.mode 必须是 baseline/plus（禁止 plus_svg_type），且必须包含 run_id/variant
    run_plus_node = _get_class_method(src, "CytoSPACEBackend", "run_plus")
    meta_dict = _extract_meta_dict(run_plus_node, method_name="run_plus")
    meta_keys = _dict_str_keys(meta_dict)

    if "mode" not in meta_keys:
        _die("FAIL: meta missing key 'mode' in run_plus")
    mode_val = meta_keys["mode"]
    if isinstance(mode_val, ast.Constant) and mode_val.value == "plus_svg_type":
        _die("FAIL: meta['mode'] must not be 'plus_svg_type' (use 'plus')")

    for required_key in ("run_id", "variant", "out_dir_name"):
        if required_key not in meta_keys:
            _die(f"FAIL: meta missing key '{required_key}' in run_plus")

    # Risk#7: meta 必须包含可审计配置与 capacity 审计
    for required_key in ("config_validation", "capacity_audit", "config_effective_subset"):
        if required_key not in meta_keys:
            _die(f"FAIL: meta missing key '{required_key}' in run_plus")

    base_meta_dict = _extract_meta_dict(run_baseline_node, method_name="run_baseline")
    base_meta_keys = _dict_str_keys(base_meta_dict)
    for required_key in ("config_validation", "capacity_audit", "config_effective_subset"):
        if required_key not in base_meta_keys:
            _die(f"FAIL: meta missing key '{required_key}' in run_baseline")

    if 'run_id.startswith("plus")' not in run_plus_chunk and "run_id.startswith('plus')" not in run_plus_chunk:
        _die("FAIL: run_plus must enforce run_id.startswith('plus')")

    print("OK")
    print("TARGET:", str(target))
    print("COUNTS:", def_counts)
    print("RUN_PLUS_CALLS:", {"harden": harden_calls, "build_outputs": build_outputs_calls})
    print("UNUSED_PRIVATE_FUNCS:", 0)
    print(
        "META_RULES:",
        {
            "mode_forbidden_plus_svg_type": True,
            "has_run_id": True,
            "has_variant": True,
            "has_out_dir_name": True,
            "enforce_run_id_startswith_plus": True,
        },
    )
    print(
        "CONFIG_RULES:",
        {
            "has_validate_and_resolve_config": True,
            "run_plus_calls_validate": True,
            "run_baseline_calls_validate": True,
            "deprecated_key_usages": len(deprecated_usages),
            "meta_has_config_validation": True,
            "meta_has_capacity_audit": True,
            "meta_has_config_effective_subset": True,
        },
    )


if __name__ == "__main__":
    main()
