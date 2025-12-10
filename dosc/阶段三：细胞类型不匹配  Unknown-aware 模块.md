好，我们就把**阶段三 = 第二个创新点（细胞类型不匹配 / unknown-aware 模块）**系统地设计一遍。
 这版还是偏“策划书风格”，不是代码实现。

------

# 阶段三：细胞类型不匹配 / Unknown-aware 模块

（详细设计说明）

## 一、阶段三在整体项目中的角色

整体骨架回顾一下：

- **阶段 1**：预处理（R + Seurat），得到干净的 sc/ST 表达矩阵 + 元数据；

- **阶段 2**：SVG+HVG 动态加权插件 → 选出一批对映射最有用的基因，并为其加权；

- **✅ 阶段 3（当前）**：

  > 在“细胞类型”这一维度上检查 **scRNA 参考图谱 和 ST 空间样本 是否匹配**，
  >  识别：
  >
  > - 明显“没有出现在空间”的 cell type（sc-only）；
  > - 在 ST 中出现但 sc 参考中没对应的“未知类型 / 新类型”（st-only）；
  >    然后为后续 CytoSPACE-plus 提供一个 **unknown-aware 的类型标签 + 类型先验矩阵**。

一句话概括：

> **阶段三 = 把原始的 cell type 注释，变成一套“空间支持感知 + Unknown-aware 的新类型体系”，并把结果以表格形式喂给阶段四。**

------

## 二、阶段三的输入 / 输出规范

### 2.1 输入：来自阶段一 + 阶段二

来自 **阶段 1** 输出目录（例如 `result/stage1/`）：

- `sc_expression_normalized.csv`
  - 单细胞归一化表达矩阵，行 = `cell_id`，列 = `gene`
- `sc_metadata.csv`
  - 至少包含：
    - `cell_id`
    - `celltype`（原始类型，可包含 `unknown`）
- `st_expression_normalized.csv`
  - 空间归一化表达矩阵，行 = `spot_id`，列 = `gene`
- `st_coordinates.csv`
  - 空间坐标（`spot_id`, `x`, `y`, `in_tissue`）

来自 **阶段 2** 输出目录（例如 `result/stage2_svg/`）：

- `plugin_genes.txt`
  - 阶段二选择的插件基因集合 `G_plugin`
- `gene_weights.csv`
  - 至少包含 `gene` + `final_weight`（可用于加权计算相似度）

> 阶段三可以选择“只在 `G_plugin` 上计算类型匹配”，这样保证和阶段二是一致的特征空间，也减轻计算量。

（可选）来自 **CytoSPACE baseline 运行结果**：

- 如果我们在某个阶段已经跑过一次 baseline CytoSPACE，并得到：
  - `spot × type` 的 cell type fraction 矩阵；
  - 也可以当作额外信息，辅助判定某些类型是否在 ST 中出现。
     设计上先不强依赖，可留作 V2 加强。

------

### 2.2 输出：给阶段四 / 评估用的产物

设阶段三输出目录为 `result/stage3_typematch/`，主要文件：

1. `type_support.csv`

   - 每个原始 cell type 一行，用于量化“这个类型在空间中是否有支持”：

     | column              | 说明                                   |
     | ------------------- | -------------------------------------- |
     | `orig_type`         | 原始 celltype 名称                     |
     | `n_cells`           | 该类型细胞数                           |
     | `support_score`     | 该类型在 ST 中的最大或平均相似度       |
     | `support_category`  | `strong` / `weak` / `unsupported` 等级 |
     | `mapped_st_cluster` | 与之最相似的 ST cluster（后面会解释）  |

2. `cell_type_relabel.csv`

   - 每个单细胞一行，用于更新 celltype 标签（unknown-aware）：

     | column        | 说明                                           |
     | ------------- | ---------------------------------------------- |
     | `cell_id`     | 细胞 ID                                        |
     | `orig_type`   | 原始类型                                       |
     | `plugin_type` | 插件后的类型标签（可能等于原始或重编码）       |
     | `status`      | 如 `kept`, `downweighted`, `unknown_merged` 等 |

   - 这里 `plugin_type` 就是 **阶段三之后，将要交给 CytoSPACE-plus 使用的类型标签**。

3. `type_prior_matrix.csv`

   - 维度：spot × plugin_type

   - 含义：**每个 spot 对各“插件类型”的匹配强度（可看成 soft 先验）**，供阶段四构建 CytoSPACE 的“类型约束”。

   - 比如某行：

     | spot_id | B_cell | T_cell | Myeloid | Unknown |
     | ------- | ------ | ------ | ------- | ------- |
     | spot_01 | 0.60   | 0.20   | 0.10    | 0.10    |

   - 数值可以是归一化后的相似度（soft assignment），不必是严格的 cell 数量。

4. `typematch_config.json`（可选）

   - 记录本次运行时使用的阈值参数、相似度方法等，方便复现实验。

------

## 三、核心思路（先用直话说一遍）

可以用一个直觉来理解阶段三：

> 我们手里有一堆 **“参考 cell type”**（来自 scRNA），
>  还有一幅 **“空间表达地图”**（ST），
>  我们要回答三个问题：

1. **哪些参考类型真的在这张空间地图里出现了？**
    → 这种类型在 ST 某些区域的基因表达，很像它在 sc 中的平均表达 → “strong support”。
2. **哪些类型几乎完全匹配不到 ST？**
    → 在任何 ST 区域都找不到和它相似的表达模式 → “unsupported” → 可以考虑合并进 `Unknown_sc_only`。
3. **ST 中是否有一些“区域表达模式”，没办法匹配到任何已知 cell type？**
    → 这些区域的表达向量和所有已知类型都不太像 → “ST-only / Unknown” → 需要在先验矩阵里留出一个 `Unknown` 类型来承接它们。

**阶段三做的事情：**

- 先在 **插件基因空间（G_plugin）** 上，为每个 sc 类型和每个 ST 区域（或 cluster）建立表达画像；
- 计算“sc-type ↔ ST-cluster/spot”的相似度矩阵；
- 基于这个矩阵给每个 cell type 打一个 **空间支持度分数**：
  - 高 → 保留，标签不变；
  - 中 → 可以保留，但在先验中权重稍微减小；
  - 低 → 认为在 ST 中基本不存在 → 重编码到 `Unknown_sc_only`；
- 同时，根据某些 ST cluster 对所有 type 都低相似度，推断出 **“ST-only 空间模式”**，在 `type_prior_matrix` 中分配给一个 `Unknown` 类型。

------

## 四、详细算法流程

下面用 Step 0~5 分解整个模块逻辑（实现语言默认为 Python）。

------

### Step 0：数据加载 & 特征准备

1. 从阶段一目录读取：

   - `sc_expression_normalized.csv` → `sc_expr`
   - `st_expression_normalized.csv` → `st_expr`
   - `sc_metadata.csv` → `sc_meta`
   - `st_coordinates.csv` → `st_coords`

2. 从阶段二目录读取：

   - `plugin_genes.txt` → `G_plugin`
   - `gene_weights.csv` → `w_gene`（`final_weight`）

3. 对 sc/ST 表达矩阵进行列子集：

   ```python
   genes = sorted(set(G_plugin) & set(sc_expr.columns) & set(st_expr.columns))
   ```

4. 选定对齐后的特征矩阵：

   - `sc_X = sc_expr[genes]`
   - `st_X = st_expr[genes]`

> 这样，阶段三的类型判断是在 **“与映射同一套基因子空间”** 上进行的，更一致。

（可选）可以在后续相似度计算时，对每个基因特征乘上 `final_weight`，形成加权特征空间。

------

### Step 1：构建参考单细胞类型画像

对于每个 sc 类型（`orig_type`）：

1. 在 `sc_meta` 中分组：

   ```python
   groups = sc_meta.groupby("celltype").indices  # 每个 celltype 对应的 cell_id 列表
   ```

2. 对每个类型 t，取该类型所有细胞的表达向量，在 `sc_X` 上求平均（或中位数），得到：

   > `profile_sc[t, g] = mean_{cells in t} sc_X[cell, g]`

3. 得到一个 **type × gene** 的矩阵 `Profile_SC`。

> 这一步相当于：
>  为每个 cell type 画出一个“平均表型”；
>  后面会拿它去和 ST 区域进行匹配。

------

### Step 2：构建 ST 空间“环境画像”（spot 或 cluster）

我们可以有两种做法：

- **简版**：直接以单个 spot 为单位，与 type 做相似度；
- **略稳一点的版**：先用表达 + 坐标对 ST spot 做一个粗的 cluster（比如 K=20~50 个 cluster），为每个 cluster 算平均 profile，减少噪音。

建议设计采用第二种（cluster 版），但保留“按 spot 直接匹配”用于构建先验矩阵。

#### 2.1 ST 聚类（可选但推荐）

1. 对 `st_X` 做一个 PCA + k-means 或 Leiden 聚类，得到 `st_cluster`：

   - 聚类数 K 可以作为参数，例如 `st_cluster_k=30`。
   - 结果保存在向量：`cluster_id[spot_id] = 0..K-1`。

2. 对每个 cluster k，计算：

   > `profile_st_cluster[k, g] = mean_{spots in k} st_X[spot, g]`

   得到矩阵 `Profile_ST_cluster`（cluster × gene）。

#### 2.2 单 spot 画像（直接用 st_X）

- 后续在构建 `type_prior_matrix` 时，会直接使用 `st_X[spot, genes]` 与 `profile_sc[type]` 做相似度，用于生成 spot-level 类型先验。

------

### Step 3：类型匹配矩阵 & 类型支持度评分

这一步要构造两个层面的相似度矩阵：

1. `S_type_cluster[t, k]`：sc 类型 t 与 ST cluster k 之间的相似度；
2. `S_type_spot[t, s]`：sc 类型 t 与单个 spot s 之间的相似度。

#### 3.1 相似度定义（建议）

- 使用 **加权 Pearson 相关 / 余弦相似度**：
  - 特征向量为：`profile_sc[t, genes]` 与 `Profile_ST_cluster[k, genes]`；
  - 若使用加权，则每个基因 g 的元素乘上 `final_weight[g]`。

伪代码（概念）：

```python
def compute_similarity(vec1, vec2, weights=None):
    # vec1, vec2: gene expression vectors
    # weights: optional gene weights
    # 输出：一个 -1~1 的相似度分数
```

#### 3.2 类型支持度：support_score

对 `S_type_cluster[t, k]`：

- 对每个类型 t，取：

  - `max_sim = max_k S_type_cluster[t, k]`
  - `mean_topK = mean of top-M similarities`（比如 M=3）

- 用一个简单的指标，比如：

  > `support_score(t) = mean_topK` 或者 `max_sim`

- 然后根据 `support_score(t)` 把类型分为：

  - `strong`: support_score ≥ θ_strong
  - `weak`: θ_weak ≤ score < θ_strong
  - `unsupported`: score < θ_weak

阈值（可作为参数）：

- `θ_strong`（比如 0.4~0.5）
- `θ_weak`（比如 0.2~0.3）

最终产出 `type_support.csv`，记录每个类型的：

- `orig_type`
- `n_cells`
- `support_score`
- `support_category`（strong / weak / unsupported）
- `mapped_st_cluster`（argmax_k S_type_cluster[t, k]）

------

### Step 4：unknown-aware 的细胞类型重标记（plugin_type）

有了每个原始类型的支持等级，我们要给每个细胞分配一个 **新类型标签 `plugin_type`**，以及一个状态标记 `status`。

#### 分类逻辑示意

对每一个原始类型 t：

- 如果 `support_category == "strong"`：
  - 所有该类型细胞保留原 label：
    - `plugin_type = t`
    - `status = "kept"`
- 如果 `support_category == "weak"`：
  - 类型 t 在 ST 中有一点信号，但不很确定：
    - 可以**保留**这个类型，但在类型先验中给予稍小权重；
    - 单细胞标签：
      - `plugin_type = t`
      - `status = "downweighted"`
- 如果 `support_category == "unsupported"`：
  - 认为这个类型基本不出现在当前 ST 样本中：
    - 把这类细胞归入一个统一的 `Unknown_sc_only`：
      - `plugin_type = "Unknown_sc_only"`
      - `status = "unknown_merged"`

特别地：

- 如果原本 `orig_type` 就是 `unknown` 或 `Unknown`：
  - 可以直接划入 `Unknown_sc_only` 或一个单独 `Unknown_from_annotation`；
  - 具体命名可在实现时统一。

最终得到一个表 `cell_type_relabel.csv`：

| cell_id | orig_type | plugin_type     | status         |
| ------- | --------- | --------------- | -------------- |
| c1      | T_cell    | T_cell          | kept           |
| c2      | rare_type | Unknown_sc_only | unknown_merged |
| c3      | B_cell    | B_cell          | downweighted   |
| ...     | ...       | ...             | ...            |

> 重点是：
>
> - 阶段四不会再直接使用 `orig_type`，而是使用 `plugin_type`；
> - `Unknown_sc_only` 在类型空间中就是一个明确的额外类别。

------

### Step 5：构建 CytoSPACE-plus 的类型先验矩阵（spot × plugin_type）

CytoSPACE 在论文中是利用 spot 的细胞类型比例先验 + UMI 计算约束来构建优化问题的。
 我们在不改内部算法的前提下，可以通过构造一个 **spot × plugin_type 的“类型兼容矩阵”** 来影响它的行为，这个矩阵就是 `type_prior_matrix.csv`。

#### 5.1 spot-level 相似度矩阵 S_type_spot

对每个 spot s 与每个类型 t（这里 t 用 `plugin_type`）：

1. 从阶段一读取：
   - `x_s = st_X[s, genes]`
   - 对应 sc 类型画像：
     - 如果 t 是普通类型：`profile_sc[t, genes]`
     - 如果 t == `Unknown_sc_only`：
       - 可以简单定义为所有被合并到 Unknown_sc_only 的细胞的平均；
       - 或者先不强求 Unknown 有明确 profile，只在人口统计层面存在。
2. 计算相似度 `S_type_spot[t, s]`，同样可以用加权余弦或相关。

#### 5.2 softmax 归一化，得到先验分布

对每个 spot s：

- 对所有 plugin_type t，取 `sim[t] = max(S_type_spot[t, s], 0)`（截断负相似度）；

- 再做一个 softmax 或归一化：

  ```python
  prior[t, s] = sim[t] / sum_t sim[t]
  ```

- 如果 sum_t sim[t] 非常小（说明这个 spot 对所有已知类型都不太像），
   → 则可以把一部分质量强制分配给 `Unknown` 类别（比如至少 0.3），其余按比例给其他类型。

最终得到 `type_prior_matrix`：

- 行：`spot_id`
- 列：所有 `plugin_type`（包括 `Unknown_sc_only`，必要时还可以有 `Unknown_st_only`）
- 值：每个 spot 对各类型的先验兼容程度（0~1，行和近似 1）。

写成 `type_prior_matrix.csv` 输出。

> 阶段四在调用 CytoSPACE-plus 时：
>
> - 使用 `plugin_type` 作为每个 cell 的类型标签；
> - 使用这一矩阵作为 spot × type 的“期望比例”的软约束（可以映射到算法允许的约束接口上）。

------

## 五、与 main.py / 阶段四的接口（概念版）

### 5.1 在 main.py 中的接口

可以在 `main.py` 中增加如下参数：

- `--enable_typematch_plugin`：是否启用阶段三（bool）
- `--typematch_strong_th`：强支持阈值 θ_strong（默认如 0.4）
- `--typematch_weak_th`：弱支持阈值 θ_weak（默认如 0.2）
- `--st_cluster_k`：ST 聚类数 K（默认如 30）

调用示意：

```python
from src.stage3_typematch import run_stage3_typematch

if args.enable_typematch_plugin:
    stage3_dir = os.path.join(args.work_dir, "stage3_typematch")
    run_stage3_typematch(
        stage1_dir=stage1_dir,
        stage2_dir=stage2_dir,
        out_dir=stage3_dir,
        strong_th=args.typematch_strong_th,
        weak_th=args.typematch_weak_th,
        st_cluster_k=args.st_cluster_k,
    )
```

------

### 5.2 阶段四如何使用阶段三的结果

阶段四（映射）需要使用两部分关键信息：

1. **单细胞的 `plugin_type` 标签**
   - 从 `cell_type_relabel.csv` 读取；
   - 用它替代原来的 `celltype` 列，作为 CytoSPACE-plus 的类型定义。
2. **spot × plugin_type 的类型先验矩阵**
   - 从 `type_prior_matrix.csv` 读取；
   - 在构建 CytoSPACE 的优化问题时，用这个矩阵作为“spot 对不同类型的兼容程度 / 期望比例”；
   - 未被支持的类型在这里就自然权重很小甚至为零。

> 这样，**细胞类型不匹配 / unknown-aware 插件是通过“改变类型标签 + 提供新类型先验”两个通道，真实介入到映射过程中的**，而不是停留在分析层。

------

## 六、阶段三完成的判断标准（简要）

可以用以下检查表来判断阶段三是否“跑通 + 对下游有用”：

1. `result/stage3_typematch/` 中存在：
   - `type_support.csv`
   - `cell_type_relabel.csv`
   - `type_prior_matrix.csv`
2. `type_support.csv` 中每个原始类型都有合理的 `support_score` 和 `support_category`（不全是 NaN）。
3. `cell_type_relabel.csv` 中：
   - 每个 `cell_id` 都有对应的 `plugin_type`；
   - `plugin_type` 的取值集合合理（原始类型 + 至少一个 Unknown 类别）。
4. `type_prior_matrix.csv` 中：
   - 行索引覆盖所有 ST spot；
   - 列覆盖所有 `plugin_type`；
   - 每行的值大致在 0~1 内，且和接近 1。
5. 阶段四能够成功读入这些结果，并以 `plugin_type` + `type_prior_matrix` 跑完一次 CytoSPACE-plus。

------

