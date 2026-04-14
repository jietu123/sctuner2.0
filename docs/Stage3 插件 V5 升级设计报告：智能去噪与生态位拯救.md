# Stage3 插件 V5 升级设计报告：智能去噪与生态位拯救

本报告详细规划了Stage3插件的 **V5** 版本升级方案。在经历了V4.1-V4.3的统计学基础建设，并吸取了V4.4（空间热点）因过拟合导致的不稳定教训后，V5版本旨在建立一套**“更鲁棒、更符合生物学直觉”**的弱信号处理机制。

V5的核心哲学是：**双轨制（Dual-Track）**。
1.  **VIP通道**：强信号直接保留。
2.  **拯救通道**：弱信号需通过“信号提纯”或“伴生关联”证明自己，并最终通过“熵值质控”方可放行。

本设计严格遵循插件式开发原则，**完全不修改下游 Stage4 (CytoSPACE) 的代码**，仅通过精细化清洗 scRNA-seq 输入数据来提升最终映射指标。

---

## 总体架构与数据流

Stage3 V5 在**当前主线**中位于 `Stage1` 导出与 `Stage4` 空间映射之间，起到“智能安检门”的作用。仓库历史上曾规划过 `Stage2 (SVG-aware)` 分支，但该分支已退出当前稳定主线；因此 Stage3 V5 的现行输入应理解为来自 `Stage1` 的标准导出，而不是依赖 Stage2。

*   **输入**：
    *   `scRNA_annotations.csv` (原始细胞类型注释)
    *   `scRNA_expression` & `ST_expression` (表达矩阵)
*   **处理流程**：
    *   **Base Check (V4.3逻辑)**：Fisher Z 检验 P < 0.05 ? $\rightarrow$ Yes: **Confident Keep**; No: **Grey Zone**.
    *   **V5.1 (信号提纯)**：针对 Grey Zone，剔除噪音 Marker 重算 P 值 $\rightarrow$ 拟拯救候选 A。
    *   **V5.2 (生态位关联)**：针对 Grey Zone，寻找与 Confident Keep 类型的空间共现 $\rightarrow$ 拟拯救候选 B。
    *   **V5.3 (熵值质控)**：针对 候选 A & B，计算空间分布熵。熵过高（太糊）则驳回。
*   **输出**：
    *   `scRNA_annotations_v5_cleaned.csv`：**经过清洗和重标注的最终名单**。
    *   Stage4 (CytoSPACE) 读取此文件，只会看到“Confident Keep”和“成功拯救”的细胞，不存在的细胞已被移除或合并。

---

## Stage3 插件 V5.1：信号提纯拯救 (Signal Denoising Rescue)

### 模块目标
解决因 **Marker 基因纯度不高**（包含高背景噪音基因）导致真实存在的稀有细胞被误判为缺失的问题。通过迭代优化 Marker 组合，挖掘被噪音掩盖的真实信号。

### 核心算法设计

1.  **识别目标**：锁定在 V4.3 基础检测中 P 值处于 **灰色区间 (0.05 $\le$ P < 0.2)** 的细胞类型集合 $C_{grey}$。
2.  **噪音贡献度分析**：
    对于每个 
    $$
    c \in C_{grey}
    $$
    ，获取其用于计算 Fisher Z 的原始 Marker 列表 
    $$
    M_{orig}（假设 N 个基因）
    $$
    。
    计算每个基因 
    $$
    g \in M_{orig}
    $$
     的**噪音分 (Noise Score)**：
    *   方法：计算基因 $g$ 在空间数据中的变异系数（CV）或与全片平均背景的偏差。偏差越小（越像背景），噪音分越高。
3.  **迭代修剪 (Pruning)**：
    
    * 剔除噪音分最高的 
      $$
      Top k\%
      $$
      个基因（例如 20%），得到 $M_{clean}$。
    *   或者采用“Leave-worst-out”策略：尝试剔除表现最差的基因，直到剩余基因数低于下限（如 5个）。
4.  **重算验证**：
    *   使用 $M_{clean}$ 重新计算 Fisher Z 支持度和 P 值，记为 $P_{new}$。
    *   **判定规则**：如果 $P_{new} < 0.05$（显著性大幅提升），则标记该类型为 **Candidate_A (拟拯救-提纯)**。

### 配置项 (YAML)
```yaml
Stage3:
  V5_denoising:
    enable: true
    p_value_upper_limit: 0.2     # 只有P值小于此值的才尝试拯救
    max_pruning_ratio: 0.3       # 最多允许剔除 30% 的原始 Marker
    min_markers_left: 10         # 剔除后至少保留 10 个基因
```

---

## Stage3 插件 V5.2：生态位关联拯救 (Co-occurrence Niche Rescue)

### 模块目标
解决 **弥散分布且信号微弱** 的细胞（如浸润免疫细胞）因自身自相关性弱而被漏检的问题。利用生物学上的“伴生关系”进行拯救。

### 核心算法设计

1.  **确立锚点 (Anchors)**：
    定义集合 $S_{sure}$ 为 V4.3 判定中 P < 0.01 的**高置信度存在**细胞类型（如肿瘤细胞、上皮细胞等结构性细胞）。
2.  **空间相关性扫描**：
    对于未被 V5.1 拯救的剩余灰色区间细胞 $c$，获取其在空间各点的支持度向量 $V_c$（长度为 Spot 数）。
    计算 $V_c$ 与 $S_{sure}$ 中每一个细胞类型 $s$ 的支持度向量 $V_s$ 的 **Spearman 相关系数** $\rho(c, s)$。
3.  **判定规则**：
    * 找出 $c$ 与所有锚点的最大相关性 
      $$
      \rho_{max} = \max_s(\rho(c, s))
      $$
      
    * 如果 
      $$
      \rho_{max} > \text{Threshold}
      $$
       (例如 0.35)，说明细胞 c 虽然信号弱，但它总是伴随着某个强信号细胞出现（处于同一生态位）。
    *   标记该类型为 **Candidate_B (拟拯救-伴生)**。
    *   *注：同时可记录“伴侣类型”信息到日志，供生物学解释。*

### 配置项 (YAML)
```yaml
Stage3:
  V5_niche_rescue:
    enable: true
    anchor_p_threshold: 0.01     # 只有极显著的类型才能当锚点(大哥)
    correlation_metric: "spearman"
    correlation_threshold: 0.35  # 弱信号必须与大哥有强相关
```

---

## Stage3 插件 V5.3：熵值质控与终极裁决 (Entropy QC & Final Verdict)

### 模块目标
**这是V5版本的核心风控模块。**
解决“过度拯救”导致 Stage4 映射指标下降的问题。通过计算空间分布的熵值，剔除那些虽然勉强通过了 V5.1/V5.2，但在空间上表现为均匀随机分布（纯噪音）的“伪信号”。

### 核心算法设计

**注意：本步骤仅针对 Candidate_A 和 Candidate_B 执行。V4.3 判定的 Confident Keep 类型豁免此检查。**

1.  **概率化映射**：
    对于每个拟拯救类型 $c$，获取其在所有 $N$ 个 Spot 上的 Fisher Z Score 向量 $Z$。
    使用 Softmax 函数将其转换为空间概率分布 $P$：
    $$
     P_i = \frac{e^{Z_i / T}}{\sum_{j=1}^{N} e^{Z_j / T}} 
    $$
    *(其中 $T$ 为温度系数，通常设为 1)*
2.  **计算归一化空间熵 (Normalized Spatial Entropy)**：
    
    *   Shannon 熵：$H(c) = - \sum_{i=1}^{N} P_i \ln(P_i)$
    *   理论最大熵（完全均匀分布）：$H_{max} = \ln(N)$
    *   **归一化熵值**：$E(c) = \frac{H(c)}{H_{max}}$
    *   $E(c)$ 的范围是 0~1。接近 1 表示全片均匀（噪音），接近 0 表示极度聚集。
3.  **终极裁决**：
    *   如果 $E(c) < \text{Entropy\_Threshold}$ (例如 0.9，表示有一定结构)：$\rightarrow$ ✅ **Rescued (Keep)**。
    *   如果 $E(c) \ge \text{Entropy\_Threshold}$ (表示太糊了)：$\rightarrow$ ❌ **Rescue Failed (Drop/Relabel)**。

### 输出对接 (CytoSPACE Interaction)
基于最终裁决结果，Stage3 生成更新后的 `annotations_v5.csv`：
*   **Keep/Rescued 类型**：保留原始 Label，不做变动。
*   **Missing/Failed 类型**：
    *   策略 A (Drop): 从 CSV 中直接删除对应细胞行。
    *   策略 B (Relabel): 将 Label 修改为 V4.3 计算出的相似类型名称（如 `CD8_T_Weak` -> `T_Cells_Global`）。
*   **Stage4 调用**：CytoSPACE 读取此 CSV，此时它看到的 scRNA-seq 数据已经没有了噪音类型，从而被迫只映射真实存在的细胞，PCC/RMSE 指标自然提升。

### 配置项 (YAML)
```yaml
Stage3:
  V5_entropy_qc:
    enable: true
    entropy_threshold: 0.90      # 超过0.9视为高熵噪音，予以剔除
  final_action:
    method: "relabel"            # 或 "drop"
    target_similarity: 0.85      # Relabel的相似度门槛
```

---

## 验收标准 (Acceptance Criteria)

### 场景一：稀有但真实的细胞 (如 少量浸润T细胞)
*   **预期行为**：
    *   V4.3 可能判为缺失 (P~0.1)。
    *   V5.1 可能通过提纯 Marker 使 P 变显著；或 V5.2 发现它与肿瘤细胞共现。
    *   V5.3 检查其空间熵，虽然弥散但有聚类倾向 (熵 < 0.9)。
    *   **结果**：成功拯救，Stage4 映射该细胞到肿瘤区域。

### 场景二：纯噪音/不存在的细胞 (如 神经元出现在肝脏)
*   **预期行为**：
    *   V4.3 判为缺失。
    *   V5.1 提纯后 P 依然不显著；V5.2 找不到任何细胞与之共现。
    *   即使勉强进入 V5.3，由于其在肝脏切片上是随机背景噪音，分布极其均匀，熵值接近 1.0。
    *   **结果**：拯救失败，执行 Drop。Stage4 不会映射该细胞，避免了“满天星”式的错误预测。

### 场景三：指标验证
*   对比开启 V5 插件前后，Stage4 (CytoSPACE) 输出的 **PCC (Pearson Correlation Coefficient)** 和 **RMSE**。
*   **合格标准**：V5 版本的 PCC 应高于 V4.3 版本，且不仅保留了主要细胞类型，还正确映射了至少一种稀有细胞类型。
