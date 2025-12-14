"""
Base backend interface for Stage4 mapping orchestrator.
"""
from __future__ import annotations

import abc
from pathlib import Path
from typing import Dict, Any


class MappingBackend(abc.ABC):
    """
    抽象基类：所有映射 backend 需实现 baseline / plus 两个入口。
    输出必须落在指定 out_dir，并写标准三件套 + meta。
    """

    def __init__(self, name: str):
        self.name = name

    @abc.abstractmethod
    def run_baseline(self, stage1_dir: Path, out_dir: Path, config: Dict[str, Any]) -> None:
        """
        baseline：不使用 Stage2/3 插件。
        """
        raise NotImplementedError

    @abc.abstractmethod
    def run_plus(
        self,
        stage1_dir: Path,
        stage2_dir: Path,
        stage3_dir: Path,
        out_dir: Path,
        config: Dict[str, Any],
    ) -> None:
        """
        plus：必须使用 Stage2/3 插件（plugin_genes/final_weight/plugin_type/type_prior）。
        """
        raise NotImplementedError
