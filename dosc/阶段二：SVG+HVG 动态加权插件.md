好，我们把之前关于阶段二的所有讨论都整理成一份**完整、成稿级的设计报告**。这份可以直接当成你项目文档里的「阶段二：SVG+HVG 动态加权插件设计说明」。

------

# 阶段二：SVG+HVG 动态加权插件

（Detailed Design Document）

## 1. 阶段二在整体项目中的角色

在整个项目链路中：

- **阶段 0**：选定 BRCA HER2+ 主数据集 + 验证 CytoSPACE 环境；

- **阶段 1**：用 R + Seurat 做统一预处理，得到干净的 sc/ST 表达矩阵及基础基因过滤信息；

- **阶段 2（当前阶段）**：**第一个创新点** ——

  > 在 sc/ST 的公共基因集合上，为每个基因计算「HVG 分数」和「SVG 分数」，再通过动态加权得到一个**最终基因权重 `final_weight(g)`**，并选出一批最有价值的基因 `G_plugin`，供后续 CytoSPACE-plus 映射使用。

与之前“非侵入式 Tangram”的区别是：

> 阶段二的结果**不是只做统计，而是实实在在改变后续映射的输入矩阵**：
>
> - 只用插件选出的基因集合 `G_plugin`；
> - 并根据 `final_weight(g)` 对每个基因加权。

因此，阶段二是我们项目中**第一个真正会改变 mapping 行为的模块化插件**。

------

## 2. 阶段二的输入与输出规范

### 2.1 输入（来自阶段一 + 原始数据）

阶段二工作在已预处理好的数据之上，主要依赖阶段一输出的文件：

来自 **阶段 1 输出目录**（例如 `result/stage1/`）：

- `sc_expression_normalized.csv`
  - 单细胞归一化表达矩阵
  - 行：`cell_id`，列：`gene`
- `st_expression_normalized.csv`
  - 空间归一化表达矩阵
  - 行：`spot_id`，列：`gene`
- `sc_metadata.csv`
  - 单细胞元数据，至少包含：
    - `cell_id`
    - `celltype`
    - `nCount_RNA`
    - `nFeature_RNA`
    - （可扩展）
- `st_coordinates.csv`
  - 空间点位信息，包含：
    - `spot_id`
    - `x`, `y`（或者 pixel 坐标）
    - `in_tissue`（0/1，可选）
- `sc_hvg_genes.txt`
  - Seurat 阶段得到的 HVG 基因列表（字符串列表）
- `st_all_genes.txt`
  - ST 预处理前的全部基因名（SVG 过滤敏感性分析用）
- `st_filtered_genes.txt`
  - ST 预处理后保留的基因名（过滤后基因集合）

来自 **原始数据目录**（可选，用于敏感性分析）：

- `brca_STdata_GEP.txt`
  - 原始 ST count 矩阵（基因 × spot）

> 这些输入由阶段一保证存在，阶段二只需要约定好的目录与文件名。

------

### 2.2 输出（给阶段四 / 评估阶段用）

阶段二的核心产物是：一套**基因权重表 + 插件基因集合**，以及一个**诊断报告**。

输出目录例如 `result/stage2_svg/`，包含：

1. `gene_weights.csv`
    表头约定如下（可扩展）：

   | 列名                 | 含义                                          |
   | -------------------- | --------------------------------------------- |
   | `gene`               | 基因名                                        |
   | `hvg_score`          | HVG 原始分数（例如 log 方差）                 |
   | `svg_score`          | SVG 原始分数（例如 Moran’s I 或类似空间指标） |
   | `hvg_norm`           | HVG 标准化分数（0~1，越大越「高变」）         |
   | `svg_norm`           | SVG 标准化分数（0~1，越大越「空间成片」）     |
   | `final_weight`       | HVG+SVG 动态融合后的最终权重                  |
   | `selected_in_plugin` | 是否被选入插件基因集合 `G_plugin`（0/1）      |
   | `is_sc_hvg`          | 是否在阶段 1 的 Seurat HVG 列表中（0/1）      |

2. `plugin_genes.txt`

   - 一行一个基因名
   - 由 `final_weight` 排序并取 Top-K 基因组成，记作 `G_plugin`

3. `svg_filter_sensitivity.json`（推荐）

   - JSON 结构，记录 **“预处理过滤是否误伤 SVG 候选基因”** 的统计信息，例如：

   ```json
   {
     "topN_svg": 500,
     "num_topN_in_filtered": 470,
     "num_topN_lost": 30,
     "lost_genes_example": ["GeneA", "GeneB", "GeneC"]
   }
   ```

> 对于后续阶段，真正要用的是：
>
> - `plugin_genes.txt`：告诉 CytoSPACE-plus 用哪些基因；
> - `gene_weights.csv`：告诉 CytoSPACE-plus 每个基因该乘多少权重。

`svg_filter_sensitivity.json` 则是偏“监控/诊断”用，为我们后续写论文时提供「SVG vs 过滤」的分析证据。

------

## 3. 基础概念与直观解释

### 3.1 HVG：细胞内在差异

- HVG（高度可变基因）反映的是：**在单细胞层面，这个基因在不同细胞间的表达差异有多大**。
- HVG 分数高 → 这个基因更能区分细胞亚群，是细胞分类的重要特征。

### 3.2 SVG：空间可变基因

- SVG 反映的是：**在空间上，这个基因是否呈现“成片的高表达/低表达区域”**。
- SVG 分数高 → 在某些区域这个基因持续偏高，在另一些区域持续偏低，说明它携带明显的空间结构信息。

可以用一个简单比喻来记：

- HVG：**在「细胞集合」里，它是不是一盏“辨别不同细胞类型的亮灯”**；
- SVG：**在「组织地图」上，它是不是在某一块区域成片亮着的灯**。

阶段二的目标，就是把这两种“灯的价值”融合起来，让我们挑选出**既能区分细胞、又能体现空间结构**的一批基因。

------

## 4. 详细算法流程

阶段二由 Python 实现（例如 `src/stage2_svg_plugin.py`），流程可以拆成 6 个步骤：

### Step 0：加载阶段一结果 + 确定公共基因集合 G

1. 从 `result/stage1/` 读入：

   ```python
   sc_expr = pd.read_csv("sc_expression_normalized.csv", index_col=0)
   st_expr = pd.read_csv("st_expression_normalized.csv", index_col=0)
   sc_meta = pd.read_csv("sc_metadata.csv", index_col=0)
   st_coords = pd.read_csv("st_coordinates.csv", index_col=0)
   sc_hvg = 读取 sc_hvg_genes.txt 为列表
   st_all_genes = 读取 st_all_genes.txt
   st_filtered_genes = 读取 st_filtered_genes.txt
   ```

2. 取 sc/ST 公共基因集合：

   ```python
   G = sorted(set(sc_expr.columns) & set(st_expr.columns))
   ```

> ✅ 之后所有 HVG/SVG/权重计算都在 G 上进行，保证上下游一致。

------

### Step 1：为每个基因计算 HVG 分数（HVG_score）

**直观目标：**
 看这个基因在**单细胞**层面是否“变化很大”。

具体做法（设计）：

1. 基于单细胞表达矩阵 `sc_expr`，对每个基因 g 计算表达方差：

   ```python
   var_sc = sc_expr[G].var(axis=0)  # 按列计算方差
   ```

2. 为了防止数值跨度太大，可以先 `log1p` 平滑，再做排序归一化：

   ```python
   import numpy as np
   
   hvg_raw = np.log1p(var_sc)                         # 原始 HVG 分数
   hvg_norm = hvg_raw.rank(method="average") / len(hvg_raw)  # 映射到 (0,1]
   ```

- `hvg_raw(g)`：表示基因 g 在所有细胞中的“变异度原始值”；
- `hvg_norm(g)`：0~1 的标准化分数，越接近 1 越“高变”。

1. 利用 `sc_hvg_genes.txt` 标记该基因是否为 Seurat HVG：

   ```python
   is_sc_hvg = hvg_raw.index.isin(sc_hvg)
   ```

> HVG 这一块不需要复杂模型，方差+标准化足够用，也能和 Seurat 的 HVG 列表形成互相印证。

------

### Step 2：为每个基因计算 SVG 分数（SVG_score）

**直观目标：**
 看这个基因在**空间**上是不是“在某一片区域持续高表达”，而不是乱跳。

#### 2.1 空间邻居图构建

1. 从 `st_coordinates.csv` 中取出每个 spot 的坐标 `(x, y)`。
2. 用 kNN（例如 k = 6 或 8）构建邻居列表：
   - 对每个 spot 找 k 个最近邻；
   - 得到邻接结构 `neighbors[spot_id] = [邻居 spot 列表]`。
   - 或者构建邻接矩阵 W（稀疏矩阵）。

> 这个邻居定义是 SVG 计算的基础，后期可调 k 作为参数。

#### 2.2 对每个基因计算空间自相关指标

对于每个基因 g：

1. 取它在所有 spots 上的表达向量 `x_g`（长度 = N_spots）。

2. 先判断是否太过稀疏：

   - 若表达非零的 spot 比例 <某个阈值（如 5%），则：
     - 直接设 `svg_score(g) = 0`（认为难以可靠判断空间结构），避免噪声过大。

3. 否则，计算一个**空间自相关指标**，可以是 Moran’s I 或类似形式：

   概念上可以理解为：

   > 如果 “表达高的 spot 的邻居也大多表达高”、“表达低的 spot 的邻居也大多表达低”，
   >  就会得到一个较大的正值；
   >  如果表达值与邻居之间毫无关系，则接近 0。

4. 得到 `svg_raw[g]` 后，同样进行排序归一化：

   ```python
   svg_norm = svg_raw.rank(method="average") / len(svg_raw)
   ```

> 最终：
>
> - `svg_score(g)`：原始空间自相关分值；
> - `svg_norm(g)`：0~1 标准化的“空间成片程度”。

------

### Step 3：HVG / SVG 动态加权融合

有了两个标准化分数：

- `hvg_norm(g)`：单细胞变异度；
- `svg_norm(g)`：空间成片度；

我们使用一个**可调配比**将其融合为一个最终权重：

> ```
> final_weight(g) = α * hvg_norm(g) + β * svg_norm(g)
> ```

其中：

- α：HVG 权重；
- β：SVG 权重；
- 默认 α = β = 0.5（各占一半）。

在实现上加一层保护（让 α+β 归一到 1）：

```python
def combine_scores(hvg_norm, svg_norm, alpha=0.5, beta=0.5):
    s = alpha + beta
    alpha, beta = alpha / s, beta / s
    final = alpha * hvg_norm + beta * svg_norm
    return final
```

> 这就是我们说的**动态加权**：
>
> - 每个基因的 `final_weight` 都是根据它在数据里的实际表现算出来的；
> - 不同数据集、不同组织会自然得到不同的权重分布；
> - 用户还能通过 main.py 参数（`--svg_alpha`、`--svg_beta`）调整“更偏细胞差异还是空间结构”。

------

### Step 4：选择插件基因集合 G_plugin

使用 `final_weight(g)` 对基因进行排名，选出真正要进入 CytoSPACE-plus 的一批基因。

#### 默认策略：Top-K 选择

- 在 main.py 中暴露参数 `--svg_topk`（例如 1000 或 2000）；
- 按 `final_weight` 从大到小排序，取前 K 个基因：

```python
def select_plugin_genes(final_weight, topk=1000):
    final_sorted = final_weight.sort_values(ascending=False)
    return list(final_sorted.index[:topk])
```

然后：

- 将结果写入 `plugin_genes.txt`；
- 在 `gene_weights.csv` 中对这些基因标记 `selected_in_plugin = 1`。

> 这样：
>
> - 即使 SVG 整体不算很强，我们也始终能选出“相对最有用的一批基因”；
> - 不依赖绝对阈值，避免“SVG 全部低→一个都用不上”的极端情况。

------

### Step 5：SVG 过滤敏感性分析（可选但推荐）

这是专门用来回答你之前那句话：

> “会不会在阶段 1 的过滤中把潜在 SVG 滤掉？”

流程设计：

1. 从 `data_dir` 再次读入原始 ST count 矩阵 `brca_STdata_GEP.txt`。
2. 对所有基因做一个简化预处理（例如 CPM + log1p），得到 `st_expr_rawnorm`。
3. 使用同样的空间邻居结构，对“所有基因”计算一个快速版 `svg_raw_all`。
4. 取 `svg_raw_all` 中 SVG 分数最高的 Top-N（比如 N=500），得到强 SVG 候选集合 `S_top`。
5. 用 `st_filtered_genes.txt` 比对：
   - `lost_genes = S_top - st_filtered_genes`
6. 汇总统计信息，写入 `svg_filter_sensitivity.json`：
   - `topN_svg`: N
   - `num_topN_in_filtered`: `len(S_top ∩ st_filtered_genes)`
   - `num_topN_lost`: `len(lost_genes)`
   - `lost_genes_example`: 列出前若干个 lost gene 名字

> 解释：
>
> - 如果 `num_topN_lost` 很少 → 说明前面过滤对 SVG 的打击很小，比较安全；
> - 如果 `num_topN_lost` 很多 → 说明过滤太狠，可以考虑回到阶段 1 把 `filter_mode` 调得更温和（例如从 `medium` 改为 `mild`），再重跑 1+2。

------

## 5. 与 main.py / 阶段四的接口设计

### 5.1 main.py 调用阶段二的参数设计

在主控脚本 `main.py` 中，为阶段二设计以下命令行参数（示例）：

- `--enable_svg_plugin`（bool）：是否启用阶段二（默认 False）
- `--svg_alpha`（float）：HVG 权重 α（默认 0.5）
- `--svg_beta`（float）：SVG 权重 β（默认 0.5）
- `--svg_topk`（int）：Top-K 基因数（默认 1000）
- `--svg_method`（str）：SVG 计算方法（先固定为 `moran`，未来可扩展）
- `--svg_sensitivity`（bool）：是否做过滤敏感性分析（默认 True）

伪调用示意：

```python
from src.stage2_svg_plugin import run_stage2_svg

if args.enable_svg_plugin:
    stage2_dir = os.path.join(args.work_dir, "stage2_svg")
    run_stage2_svg(
        data_dir=args.data_dir,
        stage1_dir=stage1_dir,
        out_dir=stage2_dir,
        alpha=args.svg_alpha,
        beta=args.svg_beta,
        topk=args.svg_topk,
        method=args.svg_method,
        do_filter_sensitivity=args.svg_sensitivity,
    )
```

> 这样，阶段二在项目里就是一个“可开关的插件”：
>
> - 开启时：生成基因权重与插件基因集合，后续映射用 plus 版本；
> - 关闭时：跳过阶段二，后续使用 baseline 版本。

------

### 5.2 阶段四如何使用阶段二结果

在 **阶段四（plus 映射）** 内部，我们约定：

1. 读取阶段一结果：`sc_expr`, `st_expr`；
2. 读取阶段二的输出：
   - `plugin_genes.txt` → `G_plugin`
   - `gene_weights.csv` → 包含 `final_weight` 列
3. 构造用于 CytoSPACE-plus 的加权表达矩阵：

```python
plugin_genes = load_plugin_genes(stage2_dir)
weights = load_gene_weights(stage2_dir)  # DataFrame: index=gene, col='final_weight'

X_sc = sc_expr[plugin_genes].values * weights.loc[plugin_genes, "final_weight"].values
X_st = st_expr[plugin_genes].values * weights.loc[plugin_genes, "final_weight"].values
```

1. 把 `X_sc`、`X_st` 传入 CytoSPACE-plus 作为表达输入，再跑后续的优化与映射。

> 这一步保证了：
>
> - **SVG+HVG 动态加权插件真实改变了映射输入**；
> - 插件关掉时，阶段四用的是原始表达；
> - 插件开启时，用的是选基因+加权表达。

------

## 6. 关键参数与默认建议

| 参数              | 含义                                        | 默认值     |
| ----------------- | ------------------------------------------- | ---------- |
| `filter_mode`     | 阶段一过滤温和程度（间接影响 SVG 候选空间） | `mild`     |
| `svg_alpha`       | HVG 权重                                    | `0.5`      |
| `svg_beta`        | SVG 权重                                    | `0.5`      |
| `svg_topk`        | 插件基因集合大小                            | `1000`     |
| `svg_method`      | SVG 评分方法（目前固定 Moran’s I 类似）     | `moran`    |
| `svg_k`           | 构建邻居图时的 kNN 邻居数（内部参数）       | 例如 `6~8` |
| `svg_sensitivity` | 是否做过滤敏感性分析                        | `True`     |

后续可以在模拟和真实数据上做敏感性实验，验证不同 α/β、topK 的影响，并给出论文中的推荐设置。

------

## 7. 阶段二的潜在风险与应对策略

1. **SVG 信号整体较弱 / 可识别 SVG 较少**
   - 可能原因：
     - 阶段一过滤过于严格；
     - 数据本身空间结构不明显；
   - 对策：
     - 调整 `filter_mode` 为更温和；
     - 使用排名+Top-K，而不是 SVG 绝对阈值；
     - 结合 `svg_filter_sensitivity.json` 评估是否严重误伤。
2. **邻居构建方式影响 SVG 结果**
   - k 太小 → 太局部、易噪；
   - k 太大 → 变成全局平均。
   - 对策：
     - 固定一个经验合理的默认值（例如 k=6 或 8）；
     - 将来可把 k 暴露为配置参数；
     - 对几个高 SVG 基因做空间可视化 sanity check。
3. **HVG 与 SVG 信息高度重叠，权重分布过于尖锐**
   - 某些 marker 基因同时 HVG 和 SVG 极高 → `final_weight` 极大；
   - 有可能让少数基因主导映射。
   - 对策：
     - 调整 α / β 平衡，减弱一侧；
     - 对 `final_weight` 做上限（cap）或平滑处理；
     - 在分析中查看 HVG_norm 与 SVG_norm 的相关性。
4. **技术噪音被误判为 SVG（如行/列效应）**
   - ST 技术噪音可能表现为空间条纹；
   - 对策：
     - 阶段一做基础 QC（UMI 正规化，去除极低质量 spot）；
     - 通过可视化检查明显“条纹”或“边缘伪影”基因；
     - 可以在后续版本中加入“技术模式回归”。
5. **计算开销**
   - 每个基因算一遍空间自相关会有一定时间成本；
   - 对策：
     - 仅在 sc/ST 公共基因集合上计算（减少规模）；
     - 跳过极度稀疏基因；
     - 必要时优化为向量化实现或使用专用库。

------

## 8. 阶段二完成的判断标准

你可以按以下 checklist 判断阶段二是否“完成且可用”：

1. 在 `result/stage2_svg/` 下成功生成：
   - `gene_weights.csv`
   - `plugin_genes.txt`
   - （若启用）`svg_filter_sensitivity.json`
2. 所有基因 g 满足：
   - `0 ≤ hvg_norm(g) ≤ 1`
   - `0 ≤ svg_norm(g) ≤ 1`
   - `final_weight(g)` 为有限实数（非 NaN）
3. `plugin_genes.txt` 中的基因全部属于 sc/ST 公共基因集合 G。
4. 在阶段四中，用 `plugin_genes.txt` + `gene_weights.csv` 成功构造加权表达矩阵，并完成至少一次 CytoSPACE-plus 运行。
5. 使用 `svg_filter_sensitivity.json` 检查时：
   - `num_topN_lost` 不至于过高；
   - 若过高，能通过调整 `filter_mode` 改善。

------

这就是一份**完整、细节齐全的阶段二设计报告**。
 如果你觉得这一版可以接受，我们就可以：

- 下一步开始设计 **阶段三：细胞类型不匹配 / unknown-aware 插件** 的详细方案；
- 或者你先挑几条（比如 SVG 计算方式、α/β 和 topK 的默认值）和我再做一轮微调，我们再锁定最终版本。