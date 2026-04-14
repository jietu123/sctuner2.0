# Current Positioning

本文件是当前项目定位的单一事实来源（SSOT）。

## 当前主线
- 方法主线：`Stage1 -> Stage3 -> Stage4 -> Stage5/6 -> Stage7`
- 核心定位：`type-aware / mismatch-aware` 的映射调谐框架
- 映射后端：当前以 CytoSPACE baseline vs route2 对比为主

## 当前核心贡献
1. 在映射前增加类型错配识别与处理（Stage3）。
2. 不改动 CytoSPACE 核心求解器，仅通过前置调谐改变输入细胞池与类型标签。
3. 形成可复核的评估链路（Stage5/Stage6/Stage7）。

## 非当前主线内容
- `Stage2 (SVG-aware / SVG+HVG)` 属于历史探索，不作为当前主 claim。
- 历史模拟场景、历史可视化批次、历史报告草稿已从仓库中清理。

## 文档使用建议
- 文档入口：`docs/README.md`
- Stage3 设计：`docs/Stage3 插件 V5 升级设计报告：智能去噪与生态位拯救.md`
- 插件接口草案：`docs/插件统一接口规范api.md`
