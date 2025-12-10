# 映射后端 & 插件统一接口规范（v0.1）

> 目标：
>
> - 让 **模块一（SVG+HVG 动态加权）** 和 **模块二（Unknown-aware 类型匹配）** 成为**通用插件**；
> - 阶段四及之后的所有映射 / 评估，只依赖统一接口，不绑死在 CytoSPACE 上；
> - 以后要接 Tangram / Cell2location / Cell2Spatial，只需实现对应 backend，**不用改插件与评估**。

------

## 1. 总体设计思路

- 阶段 2 / 3：输出**标准化的插件结果**（基因权重、类型重写、先验）。
- 阶段 4：通过一个 `MappingBackend` 抽象接口，将插件结果接入具体映射方法：
  - `CytoSPACEBackend`
  - `TangramBackend`
  - `Cell2locationBackend`
  - ...
- 所有 backend **必须输出统一格式的 mapping 结果**：
  - `cell_assignment.csv`
  - `cell_spot_matrix.npz`
  - `spot_type_fraction.csv`
- 阶段 5 / 6 以及“通用评估”只依赖这三类输出，不关心内部是 CytoSPACE 还是 Tangram。

------

## 2. 上游插件统一输出接口（阶段 2 & 3）

> 说明：**这些文件格式必须尽量保持稳定**，后端和评估都靠它们。

### 2.1 阶段二：SVG+HVG 动态加权插件

输出目录示例：`result/stage2_svg/`

#### 2.1.1 `gene_weights.csv`

- 描述：每个基因一行，记录 HVG/SVG 分数和最终权重。

必需列（后续 backend 和评估会用到）：

| 列名                 | 类型  | 说明                                 |
| -------------------- | ----- | ------------------------------------ |
| `gene`               | str   | 基因名                               |
| `hvg_norm`           | float | HVG 标准化分数（0~1）                |
| `svg_norm`           | float | SVG 标准化分数（0~1）                |
| `final_weight`       | float | 综合权重（HVG+SVG 动态加权结果）     |
| `selected_in_plugin` | int   | 是否被选入插件基因集合（0/1）        |
| `is_sc_hvg`          | int   | 是否是 Seurat HVG（0/1，可用于分析） |

后端主要使用：

- `gene`
- `final_weight`
- `selected_in_plugin`

#### 2.1.2 `plugin_genes.txt`

- 描述：插件基因集合 `G_plugin`。
- 格式：一行一个基因名。
- 生成方式：按 `final_weight` 排序取 Top-K。

> **约定：所有 backend 的 plus 模式，默认只在 `plugin_genes` 这些基因上工作**，除非有特殊说明。

------

### 2.2 阶段三：Unknown-aware 类型匹配插件

输出目录示例：`result/stage3_typematch/`

#### 2.2.1 `type_support.csv`（偏分析用，backend 可不用）

- 每个原始 celltype 一行，记录其在 ST 中的支持度。

关键列（评估时用到）：

| 列名                | 类型    | 说明                               |
| ------------------- | ------- | ---------------------------------- |
| `orig_type`         | str     | 原始 celltype 名称                 |
| `n_cells`           | int     | 该类型细胞数                       |
| `support_score`     | float   | 在 ST 中的支持度（0~1 或相关系数） |
| `support_category`  | str     | `strong` / `weak` / `unsupported`  |
| `mapped_st_cluster` | str/int | 最相似的 ST cluster / 区域         |

> 后端一般不直接用，但评估阶段可以用来分析“不匹配类型被怎么处理了”。

#### 2.2.2 `cell_type_relabel.csv`

- 描述：每个细胞的类型重标记结果，**后端在 plus 模式中必须使用**。

关键列：

| 列名          | 类型 | 说明                                          |
| ------------- | ---- | --------------------------------------------- |
| `cell_id`     | str  | 细胞 ID，与 sc 表达矩阵行名一致               |
| `orig_type`   | str  | 原始类型（阶段一的 celltype）                 |
| `plugin_type` | str  | 插件后的类型标签（含 Unknown_sc_only 等）     |
| `status`      | str  | `kept` / `downweighted` / `unknown_merged`... |

> **约定：所有 backend 在 plus 模式下，必须将 `plugin_type` 当成细胞类型标签**；
>  baseline 模式则使用 `orig_type` 或阶段一原始标注。

#### 2.2.3 `type_prior_matrix.csv`

- 描述：spot×plugin_type 的类型先验矩阵。
- 行索引：`spot_id`
- 列名：所有 `plugin_type`（包括 `Unknown_sc_only` 等）
- 单元格：`prior[s, t]`，表示该 spot 对类型 t 的支持度（0~1，行和约为 1）

> **约定：所有 backend 在 plus 模式下，应尽量利用这一矩阵作为“类型先验/约束”的来源**；
>  baseline 模式使用各自原始方法的先验构造方式。

------

## 3. 映射后端抽象接口（伪 API）

> 下面是一个 Python 伪类，只表达“接口长什么样”，不用完全照抄。

```python
class MappingBackend:
    """
    映射后端抽象接口。

    一切具体方法（CytoSPACE、Tangram、Cell2location...）
    都应该继承本类，并实现 run_baseline / run_plus。
    """

    def __init__(self, name: str):
        self.name = name  # 比如 "cytospace", "tangram"

    def run_baseline(
        self,
        stage1_dir: str,
        out_dir: str,
        config: dict | None = None,
    ) -> None:
        """
        运行 baseline 映射（不使用插件）。

        参数
        ----
        stage1_dir : 阶段一输出目录（包含 sc/st 表达 & 元数据）
        out_dir    : 本 backend 的 baseline 输出目录
        config     : (可选) backend 内部需要的超参字典

        要求
        ----
        在 out_dir 下，必须生成统一格式的结果文件：
          - cell_assignment_baseline.csv
          - cell_spot_matrix_baseline.npz (可选但强烈推荐)
          - spot_type_fraction_baseline.csv
        详见第 4 节“统一输出格式”。
        """
        raise NotImplementedError

    def run_plus(
        self,
        stage1_dir: str,
        stage2_dir: str,
        stage3_dir: str,
        out_dir: str,
        config: dict | None = None,
    ) -> None:
        """
        运行 plus 映射（接入模块一+模块二）。

        参数
        ----
        stage1_dir : 阶段一输出目录
        stage2_dir : 阶段二输出目录 (gene_weights, plugin_genes)
        stage3_dir : 阶段三输出目录 (cell_type_relabel, type_prior_matrix)
        out_dir    : 本 backend 的 plus 输出目录
        config     : (可选) backend 内部需要的超参字典

        要求
        ----
        在 out_dir 下，必须生成统一格式的结果文件：
          - cell_assignment_plus.csv
          - cell_spot_matrix_plus.npz
          - spot_type_fraction_plus.csv

        注意
        ----
        run_plus 中必须显式使用：
          - plugin_genes + gene_weights.final_weight（改变表达空间/距离）
          - plugin_type（作为细胞类型）
          - type_prior_matrix（作为 spot×type 先验）
        """
        raise NotImplementedError
```

> 后续你实现 `CytoSPACEBackend(MappingBackend)` 时，就按这个接口来填即可。

------

## 4. 映射结果统一输出格式（所有 backend 必须遵守）

> 这是后面“通用评估”的唯一依赖点，非常关键。

假设某个 backend 为 `cytospace`，则阶段四 orchestrator 可以组织结构如下：

```text
result/
  stage4_mapping/
    cytospace/
      baseline/
        cell_assignment_baseline.csv
        cell_spot_matrix_baseline.npz
        spot_type_fraction_baseline.csv
      plus_svg_type/
        cell_assignment_plus.csv
        cell_spot_matrix_plus.npz
        spot_type_fraction_plus.csv
```

对于其他 backend（如 tangram），目录结构类似：

```text
result/
  stage4_mapping/
    tangram/
      baseline/...
      plus_svg_type/...
```

下面说明三个关键文件的格式。

------

### 4.1 `cell_assignment_*.csv`

- 作用：记录「每个细胞最终被分配到哪个 spot，以及相关信息」。

#### 必需列

| 列名      | 说明                |
| --------- | ------------------- |
| `cell_id` | 细胞 ID             |
| `spot_id` | 被分配到的空间点 ID |

#### 推荐列（方便评估、画图，但不强制）

| 列名           | 说明                                       |
| -------------- | ------------------------------------------ |
| `type`         | 该细胞视角下的类型标签                     |
| `mode`         | `"baseline"` 或 `"plus"`                   |
| `backend`      | `"cytospace"`, `"tangram"` 等              |
| `assign_score` | 相似度 / 概率 / -cost 等，越大/越小越好    |
| `extra_info`   | 字符串或 JSON，留给各 backend 填充（可选） |

> 约定：
>
> - baseline：`type` 一般用原始类型（或方法内部使用的类型）；
> - plus：`type` 一般用 `plugin_type`。

------

### 4.2 `cell_spot_matrix_*.npz`（推荐）

- 作用：保存一个稀疏或 dense 的 **cell×spot 指派矩阵**，供后续精细指标使用。
- 内容约定：
  - 若是 dense：可以是一个 `numpy.ndarray`，shape = (n_cells, n_spots)；
  - 若是稀疏：可使用 `scipy.sparse` 的 CSR/CSC，然后保存为 `.npz`。
- 意义：
  - 行：某个细胞；
  - 列：某个 spot；
  - 值：指派强度（如 0/1，或概率，或 soft assignment）。

> 如果某个 backend 没法轻易给出 soft assignment，
>  至少需要给一个 one-hot 的 0/1 矩阵（单一最佳 spot）。

------

### 4.3 `spot_type_fraction_*.csv`

- 作用：表示「每个 spot 的类型组成」，用于下游评估（比如 type-level 空间分布比较）。
- 行索引：`spot_id`
- 列：
  - 第一列：`spot_id`
  - 之后每列：具体类型名（baseline 用原始类型，plus 用 plugin_type）
  - 单元格：该 spot 内该类型的占比 / 归一化 cell 数（行和可约为 1）

> 评估时，可以：
>
> - 比较 baseline vs plus 的 `spot_type_fraction` 与真值（若有）；
> - 比较这些分布的 KL / JS / L1 距离；
> - 分析 Unknown 类型在空间中的布置等。

------

对，你抓得很准——**我上一版文档只把「通用 backend 接口」说清楚了，还没把「CytoSPACE 专用验证」这条线单独标出来**，这是我没对齐好。我们现在就把“两条对照试验线”一起写清楚：

- 线 A：**只针对 CytoSPACE 的深度验证（阶段 5 / 6 的主线）**
- 线 B：**方法无关的通用验证（阶段 5B / 6B）**

下面这份是**在刚才 mapping 接口基础上的补充版 API 说明**：不改动你已经接受的那部分，只是“再加一层评估接口”。

------

## 5、评估层接口设计：CytoSPACE 专用 vs 通用版

> 前提：
>  映射阶段（阶段 4）已经按照上一个文档输出了统一的结果：
>
> - `cell_assignment_*.csv`
> - `cell_spot_matrix_*.npz`
> - `spot_type_fraction_*.csv`
>
> 下面是**阶段 5 / 6 的评估接口**，明确分成两条线：

- **实验线 A：CytoSPACE 专用评估（5A / 6A）**
- **实验线 B：方法无关的通用评估（5B / 6B）**

------

### 一、实验线 A：CytoSPACE 专用评估接口（5A / 6A）

#### 1.1 目标

- 只针对 `backend = "cytospace"`；
- 深度利用 CytoSPACE 独有的信息（count 约束、类型约束、cost、内部日志等）；
- 在模拟数据（阶段 5A）和真实数据（阶段 6A）上，对比：
  - **Baseline vs Plus（带插件）** 的表现差异；
  - 尽可能做“指标级别的精细说明：插件到底改好了什么”。

#### 1.2 额外输出约定（CytoSPACE backend 自己多导一点）

在我们前面定义的 mapping 输出基础上，**CytoSPACEBackend 可以多给一些“自己特有”的文件**（这些别的 backend 不需要）：

例如在：

```text
result/stage4_mapping/cytospace/baseline/
result/stage4_mapping/cytospace/plus_svg_type/
```

里额外输出：

1. `cyto_cost_matrix.npz`（可选但推荐）
   - cell×spot 的 cost 矩阵（至少可以给出“最终使用的 cost”，哪怕是粗略的）。
2. `cyto_constraints.json`
   - 记录 CytoSPACE 内部用到的约束信息，例如：
     - 每个 spot 的 cell count 约束（expected cells per spot）；
     - spot×type 的期望比例（baseline / plus）。
3. `cyto_objective_stats.json`
   - 比如：
     - 总 cost；
     - 约束违反程度（如果有此类统计）；
     - 优化迭代次数、收敛情况等。

> 这些只有 CytoSPACE backend 提供，**不用强迫其他方法实现**。
>  实验线 A 的评估函数可以专门利用这些信息。

#### 1.3 CytoSPACE 专用评估 API（伪代码）

我们可以这样定义两个函数：

```python
def eval_cytospace_simulation(
    mapping_root: str,
    truth_dir: str,
    out_dir: str,
):
    """
    阶段 5A：CytoSPACE + 模拟数据的深度评估。

    参数
    ----
    mapping_root : 比如 "result/stage4_mapping/cytospace"
                   要求内部有 baseline/ 和 plus_svg_type/ 两个子目录。
    truth_dir    : 模拟真值目录，包含：
                   - cell_true_spot.csv  (cell_id, true_spot_id)
                   - spot_true_type_fraction.csv (可选)
    out_dir      : 评估结果输出目录

    核心逻辑
    --------
    1. 分别读取 baseline 和 plus 的：
         - cell_assignment_*.csv
         - cell_spot_matrix_*.npz (如果有)
         - spot_type_fraction_*.csv
         - cyto_cost_matrix.npz (如果有)
         - cyto_constraints.json, cyto_objective_stats.json (如果有)
    2. 读取真值：
         - cell_true_spot.csv
         - spot_true_type_fraction.csv
    3. 计算一组“偏 CytoSPACE”的深度指标，例如：
         - cell-level mapping accuracy / top-k accuracy
         - spot-level 类型组成 vs 真值的 KL/JS/L1
         - 总 cost & 每细胞平均 cost 的变化 (baseline vs plus)
         - 约束满足度：实际 assigned cells vs expected per spot/type
         - cost 分布对比 (histogram)
    4. 把所有指标写到：
         - metrics_cytospace_sim_baseline_vs_plus.json
         - 并生成若干可视化图 (png/pdf)
    """
    ...
def eval_cytospace_realdata(
    mapping_root: str,
    out_dir: str,
):
    """
    阶段 6A：CytoSPACE + 真实数据（无真值）的深度评估。

    逻辑类似，但不再依赖真值，而是侧重：
      - 邻域一致性 / 空间聚集度 (type-wise spatial autocorrelation)
      - marker gene 表达 vs 映射类型空间分布的一致性
      - cost 分布、约束满足情况 (baseline vs plus)
      - 高 SVG 基因热点 vs 类型热点的共定位程度
    """
    ...
```

> 总结：
>
> - 这两个函数只为 `backend="cytospace"` 服务；
> - 可以大胆使用 CytoSPACE 特有的 cost、约束等信息；
> - 在论文中作为“**CytoSPACE case study**”部分，展示插件的细节提升。

------

### 二、实验线 B：方法无关的通用评估接口（5B / 6B）

这一部分其实就是我上一份文档里重点讲的那块，我们现在**明确标注它是“第二条验证线”**。

#### 2.1 目标

- 面向任意 backend（cytospace, tangram, cell2location...）；
- 完全只依赖统一的 mapping 输出三件套：
  - `cell_assignment_*.csv`
  - `cell_spot_matrix_*.npz`
  - `spot_type_fraction_*.csv`
- 不使用任何“内部私货”（比如 cost、约束实现细节），从而做到**方法无关**。
- 在模拟数据（5B）和真实数据（6B）上，对比 baseline vs plus 的通用指标。

#### 2.2 通用评估 API（伪代码）

```python
def eval_backend_simulation_generic(
    backend_name: str,
    mapping_root: str,
    truth_dir: str,
    out_dir: str,
):
    """
    阶段 5B：通用模拟数据评估（方法无关）。

    参数
    ----
    backend_name : "cytospace", "tangram", "cell2location", ...
    mapping_root : 比如 "result/stage4_mapping"
                   下级目录应为 f"{backend_name}/baseline" 和 ".../plus_svg_type"
    truth_dir    : 模拟真值目录，格式与 eval_cytospace_simulation 类似
    out_dir      : 评估结果输出目录

    核心逻辑
    --------
    1. 从 mapping_root 中，读取该 backend 的 baseline & plus：
         - cell_assignment_baseline.csv
         - cell_spot_matrix_baseline.npz
         - spot_type_fraction_baseline.csv
         - cell_assignment_plus.csv
         - cell_spot_matrix_plus.npz
         - spot_type_fraction_plus.csv
    2. 加载真值：
         - cell_true_spot.csv
         - spot_true_type_fraction.csv
    3. 只用这些“统一结果”，计算 method-agnostic 指标，例如：
         - cell-level accuracy / top-k accuracy (baseline vs plus)
         - spot-level type 分布 vs 真值的距离 (KL/JS/L1)
         - 邻域结构一致性 (graph-based metrics)
         - 类型层面的空间分布 overlap / correlation
    4. 输出统一格式的指标文件：
         - f"metrics_generic_sim_{backend_name}_baseline_vs_plus.json"
    """
    ...
def eval_backend_realdata_generic(
    backend_name: str,
    mapping_root: str,
    out_dir: str,
):
    """
    阶段 6B：通用真实数据评估（方法无关）。

    不依赖真值，且不依赖任何 backend 内部的 cost / 约束，只用 mapping 输出结果。

    可以计算的指标例子：
      - 同类型细胞的空间聚集/分散（Moran's I / Ripley's K 等）
      - marker gene 热点 vs 映射类型热点的空间相关性
      - SVG 基因空间模式 vs 预测 type 分布的相关性
      - baseline vs plus 在这些指标上的明显提升情况

    评估结果写为：
      - f"metrics_generic_real_{backend_name}_baseline_vs_plus.json"
      - 若干通用可视化图
    """
    ...
```

> 这里你可以看到：
>
> - **CytoSPACE 也可以参加“通用评估”这一条线**（它既有 A 线也有 B 线）；
> - 其他 backend（Tangram 等）只参加 B 线（因为它们没有 A 线里那些 CytoSPACE 私货）。

------

### 三、两条验证线是怎么配合的？

可以简单用一个小表总结：

|         | CytoSPACE backend  | 其他 backend (Tangram, ...) |
| ------- | ------------------ | --------------------------- |
| 5A / 6A | ✅ 专用评估（深度） | ❌ 不适用                    |
| 5B / 6B | ✅ 通用评估         | ✅ 通用评估                  |

也就是说：

- 在论文里，你可以有这样的结构：
  1. **Section X：Case study on CytoSPACE**
     - 用 5A / 6A 的结果，深度展示插件对一个强 baseline 的改善细节。
  2. **Section Y：Method-agnostic evaluation**
     - 用 5B / 6B 的结果，横向比较多个 backend，在相同通用指标下展示插件的普适性。

------

### 四、回到你最初那句话，明确对齐一下

> “我们不是说好这 2 种对照试验吗？一种是针对 CytoSPACE 的；一种是通用的。”

现在可以这样理解我们已经设计好的结构：

- 映射层（阶段 4）：
  - 已经通过 `MappingBackend` 抽象 + 统一输出格式，实现了**方法无关**；
  - CytoSPACE 只是实现了一个 `CytoSPACEBackend`。
- 评估层（阶段 5 / 6）：
  - **线 A（5A/6A）**：`eval_cytospace_simulation` / `eval_cytospace_realdata`，只给 CytoSPACE 用；
  - **线 B（5B/6B）**：`eval_backend_simulation_generic` / `eval_backend_realdata_generic`，所有 backend 都可用。

