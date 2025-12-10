![](images/CytoSPACE_logo.jpeg)

# 

单细胞和空间转录组的稳健且快速比对
单细胞和空间转录组的稳健且快速比对

**CytoSPACE** 是一种新颖的计算工具，用于将单细胞转录组分配到原位空间转录组（ST）数据中。我们的方法通过最短增强路径优化程序，最小化基于相关性的代价函数来解决单细胞/点点分配问题。

![](images/CytoSPACE_overview.png)

我们方法的关键创新包括：
我们方法的关键创新包括：

*   与传统通过点点计算细胞类型分解的方法不同，CytoSPACE生成的组织标本既具备高基因覆盖率，又具有空间分辨的scRNA-seq数据，适合后续分析。
    与传统通过点点计算细胞类型分解的方法不同，CytoSPACE 生成的组织标本既具备高基因覆盖率，又具有空间分辨的 scRNA-seq 数据，适合后续分析。
*   CytoSPACE对噪声具有高度鲁棒性，并返回全球最优的细胞到斑点分配。（详见[论文](https://www.nature.com/articles/s41587-023-01697-9)。）
*   与其他通常在预选标记基因或共享嵌入空间上作的方法（后者可能抹去真实的生物变异）不同，CytoSPACE使用完整的转录组而无需批次校正，有助于保持对细微细胞状态的敏感性。
    与其他通常在预选标记基因或共享嵌入空间上作的方法（后者可能抹去真实的生物变异）不同，CytoSPACE 使用完整的转录组而无需批次校正，有助于保持对细微细胞状态的敏感性。
*   CytoSPACE执行快速且简单。即使用个人笔记本的单CPU也能在几分钟内运行，无需超参数调优或基因/特征选择。
    CytoSPACE 执行快速且简单。即使用个人笔记本的单 CPU 也能在几分钟内运行，无需超参数调优或基因/特征选择。

CytoSPACE 通过 [cytospace.stanford.edu](https://cytospace.stanford.edu/) 的网页界面提供，用户可以在无需下载源代码的情况下以默认设置运行 CytoSPACE。

## 安装说明（5-10分钟）
安装说明（5-10分钟）

展开部分

1.  如果还没有[Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)，可以安装。
    
2.  克隆这个仓库：
    克隆这个仓库：
    

```bash
  git clone https://github.com/digitalcytometry/cytospace
```

3.  导航至目录：`cytospace`
    导航至目录：`cytospace`

```bash
  cd cytospace
```

4.  （5-10分钟）创建一个带有所需依赖的 conda 环境：
    （5-10 分钟）创建一个带有所需依赖的 conda 环境：

```bash
  conda env create -f environment.yml
```

5.  激活你刚刚创建的环境：
    激活你刚刚创建的环境：

```bash
  conda activate cytospace_v1.1.0
```

`cytospace_v1.1.0`是最新版本文件中的环境名称。如果你有不同的版本，可以在文件的第一行检查，在命令里用你的版本名替换。`environment.yml``cytospace_v1.1.0``conda activate`
`cytospace_v1.1.0` 是最新版本文件中的环境名称。如果你有不同的版本，可以在文件的第一行检查，在命令里用你的版本名替换。`environment.yml``cytospace_v1.1.0``conda activate`

注意：如果你运行在 M1/M2 Mac上，解决 conda 环境时可能会因为安装 conda 依赖时遇到错误。要解决这个问题，请尝试以下方法：
注意：如果你运行在 M1/M2 Mac 上，解决 conda 环境时可能会因为安装 conda 依赖时遇到错误。要解决这个问题，请尝试以下方法：

```bash
    conda create -n cytospace_v1.1.0
    conda activate cytospace_v1.1.0
    conda config --env --set subdir osx-64
    conda env update --file environment.yml
```

6.  （~30秒）通过执行以下方式安装 CytoSPACE ：
    （~30 秒）通过执行以下方式安装 CytoSPACE ：

```bash
  pip install .
```

7.  （推荐，~1分钟）通过执行以下方式安装包：`lapjv`
    （推荐，~1 分钟）通过执行：`lapjv` 安装套件

```bash
  pip install lapjv==1.3.14
```

我们强烈建议您安装此软件包，它在 CytoSPACE 中快速实现了默认的核心优化算法。
我们强烈建议您安装此软件包，它在 CytoSPACE 中快速实现了默认的核心优化算法。

**第7步可能会**因系统不同而出错，因为该软件包需要CPU支持AVX2指令，而并非所有系统都支持。要判断你的系统是否支持该软件包，通常最简单的方法是尝试安装它。如果安装顺利，你的系统就能支持它！如果你遇到错误，很可能你的系统不支持它，你可以在调用CytoSPACE时指定（）参数，使用我们提供的其他选项。详情请参见[**“运行CytoSPACE**](#running-cytospace) - **选择求解器**”。`--solver-method``-sm`
第 7 步可能会根据你的系统出现错误，因为该软件包需要 CPU 支持 AVX2 指令，而并非所有系统都支持。要判断你的系统是否支持该软件包，通常最简单的方法是按照上述方法尝试安装。如果安装顺利，你的系统会支持的！如果遇到错误，很可能是你的系统不支持，你可以在调用 CytoSPACE 时指定 （） 参数，使用我们提供的其他选项。详情请参见运行 CytoSPACE - 选择求解器。`--求解器-方法` `-SM`

请注意，如果软件包成功安装，但在运行 CytoSPACE 时遇到“非法指令”错误，你可能可以通过以下命令构建该软件包：
请注意，如果软件包成功安装，但在运行 CytoSPACE 时遇到“非法指令”错误，你可能可以通过以下命令构建该软件包：

```bash
   pip3 install git+https://github.com/src-d/lapjv
```

更多信息请参见[lapjv文档页面](https://pypi.org/project/lapjv/)。

## 输入文件

展开部分

默认情况下，CytoSPACE 需要 4 个文件作为输入。所有文件应以制表符分隔的表格输入格式（.txt）提供，且不得加双引号。每个输入文件的更多格式细节如下所述。在本节末尾，我们还提供了使用脚本从Seurat对象生成输入文件的指导。
默认情况下，CytoSPACE 需要 4 个文件作为输入。所有文件应以制表符分隔的表格输入格式（.txt）提供，且不得加双引号。每个输入文件的更多格式细节如下所述。在本节末尾，我们还提供了使用脚本从 Seurat 对象生成输入文件的指导。

1.  **一个scRNA-seq基因表达文件：**

*   矩阵必须由基因（行）组成，按细胞（列）排列。
    矩阵必须由基因（行）组成，按细胞（列）排列。
*   第一行必须包含单个细胞ID，第一列必须包含基因名称。
    第一行必须包含单个细胞 ID，第一列必须包含基因名称。
*   第一列（基因名）必须有标题。
    第一列（基因名）必须有标题。
*   基因表达数据应以非归一化计数表示。
    基因表达数据应以非归一化计数表示。
*   所有重复基因名称的实例将在运行时被剔除。
    所有重复基因名称的实例将在运行时被剔除。

![](images/scRNAfile.png)

2.  **一个单元类型标签文件：**

*   对应scRNA-seq基因表达矩阵中单个细胞ID的细胞类型标记。
    对应 scRNA-seq 基因表达矩阵中单个细胞 ID 的细胞类型标记。
*   单元格类型的标签字符串不应包含特殊字符。
    单元格类型的标签字符串不应包含特殊字符。
*   表格应包含两列，第1列包含对应scRNA-seq矩阵列的单个单元ID，第2列包含对应的单元类型标签。
    表格应包含两列，第 1 列包含对应 scRNA-seq 矩阵列的单个单元 ID，第 2 列包含对应的单元类型标签。
*   列必须有头部。
    列必须有头部。

![](images/celllabelfile.png)

3.  **空间转录组学（ST）基因表达文件：**

*   矩阵必须由基因（行）组成，并按ST点（列）。
    矩阵必须由基因（行）组成，并按 ST 点（列）。
*   第一行必须包含ST点ID，第一列必须包含基因名称。
    第一行必须包含 ST 点 ID，第一列必须包含基因名称。
*   第一列（基因名）必须有标题。
    第一列（基因名）必须有标题。
*   基因表达数据应以非归一化计数表示。
    基因表达数据应以非归一化计数表示。
*   所有重复基因名称的实例将在运行时被剔除。
    所有重复基因名称的实例将在运行时被剔除。

![](images/STdatafile.png)

4.  **空间转录组坐标文件：**

*   一个由3列组成的表格，第一列包含对应ST基因表达矩阵列的ST点ID，第2列和第3列分别包含每个ST点的行和列索引。
    一个由 3 列组成的表格，第一列包含对应 ST 基因表达矩阵列的 ST 点 ID，第 2 列和第 3 列分别包含每个 ST 点的行和列索引。
*   列必须有头部。
    列必须有头部。

![](images/STcoordfile.png)

### 来自太空游侠的输出
来自太空游侠的输出

如果用户从Space Ranger输出开始，可以将ST输入文件作为单一tar.gz，**替代**（3）基因表达和（4）坐标。如果指定了Space Ranger输出，CytoSPACE会自动尝试解压提供的tarball，并加载对应的ST表达和坐标数据。如果基因表达文件中有重复的基因名称，它们将被设置为唯一，正如Seurat所做的那样。

沥青球应仅包含以下内容：
沥青球应仅包含以下内容：

*   一个包含ST基因表达的单个H5文件（扩展名.h5）
    一个包含 ST 基因表达的单个 H5 文件（扩展名.h5）
*   一个包含图像数据的单一子目录
    一个包含图像数据的单一子目录

将上述项目置于名为 的目录中，可以通过以下命令生成一个 tarball：`spaceranger_input`
将上述项目放在名为 的目录中，可以通过以下命令生成一个 tarball： `spaceranger_input`

```bash
  tar -cvzf sr_input.tar.gz spaceranger_input
```

或者更广义地说：
或者更广义地说：

```bash
  tar -cvzf [name_of_tarball] [name_of_directory]
```

左侧展示了一个未拉链沥青球的文件树示例。如果从公开的10倍Visium数据下载，用户可以下载右侧所示的文件。
左侧展示了一个未拉链沥青球的文件树示例。如果从公开的 10 倍 Visium 数据下载，用户可以下载右侧所示的文件。

![](images/VisiumTar.png) ![](images/Visium.png)

从 CytoSPACE v1.1.0 开始，该工具支持 Space Ranger v2.0.0+ 输出。这通过检测具有更新命名规范的文件（即文件本身），然后将这些文件的副本保存到与 scipy 读取 Visium 数据集功能的命名规范兼容的目录中实现。`tissue_positions.csv`
从 CytoSPACE v1.1.0 开始，该工具支持 Space Ranger v2.0.0+ 输出。这通过检测具有更新命名规范的文件（即文件本身），然后将这些文件的副本保存到与 scipy 读取 Visium 数据集功能的命名规范兼容的目录中实现。`tissue_positions.csv`

### 使用稀疏矩阵进行基因表达
使用稀疏矩阵进行基因表达

从 CytoSPACE v1.0.5 起，用户还可以以稀疏矩阵形式提供 scRNA 或 ST 基因表达，而非用制表符或逗号分隔的文件。这将需要非常特定的文件名集合，以避免解析问题，我们建议用户使用下面的 R 辅助脚本生成这些输入。请参阅下面的小节，了解更多关于这些辅助脚本的信息。`Prepare_input_files`
这将需要非常特定的文件名集合，以避免解析问题，我们建议用户使用下面的 R 辅助脚本生成这些输入。 请参阅下面的小节，了解更多关于这些辅助脚本的信息。`Prepare_input_files`

如果以稀疏矩阵形式输入，你需要在同一目录下三个文件来表示一个表达式矩阵：、 、 和 。文件会以稀疏矩阵格式列出数值，也可以是压缩格式。括号内的文件名（）可能有所不同。和 将是矩阵对应的基因名和细胞（或点/条形码）名称，位于与矩阵文件相同的目录中。`[expression].mtx``genes.tsv``cells.tsv``.mtx``[expression].mtx.gz``expression``genes.tsv``cells.tsv`
如果以稀疏矩阵形式输入，你需要在同一目录下三个文件来表示一个表达式矩阵：、 、 和 。文件会以稀疏矩阵格式列出数值，也可以是压缩格式。括号内的文件名（）可能有所不同。和 将是矩阵对应的基因名和细胞（或点/条形码）名称，位于与矩阵文件相同的目录中。`[表达].mtx` `基因.tsv` `细胞.``tsv.mtx``[表达].mtx.gz` 表达`基因.tsv``cells.tsv`

从 CytoSPACE v1.1.0 开始，该工具还支持基于 Space Ranger 和 Cell Ranger 常见稀疏矩阵输出的替代命名规范（和）。这些文件也支持扩展名和压缩格式。文件每行应有一个条目（无头部）。如果文件中有多个列，CytoSPACE会将第一列作为基因/细胞名。请注意，对于常见的Cell Ranger和Space Ranger稀疏矩阵输入，基因文件的第一列通常包含Ensembl ID，而基因符号在第二列，因此如果你的其他输入模态包含基因名称的基因符号，你需要更改列的顺序。`features.tsv``barcodes.tsv``.csv``.gz`
从 CytoSPACE v1.1.0 开始，该工具还支持基于 Space Ranger 和 Cell Ranger 常见稀疏矩阵输出的替代命名规范（和）。这些文件也支持扩展名和压缩格式。文件每行应有一个条目（无头部）。如果文件中有多个列，CytoSPACE 会将第一列作为基因/细胞名。请注意，对于常见的 Cell Ranger 和 Space Ranger 稀疏矩阵输入，基因文件的第一列通常包含 Ensembl ID，而基因符号在第二列，因此如果你的其他输入模态包含基因名称的基因符号，你需要更改列的顺序。`features.tsv` `条码.tsv``.csv``.gz`

文件应作为或的参数提供，此时CytoSPACE会自动尝试从同一目录中查找对应的基因和细胞文件。`[expression].mtx``--st-path``--scRNA-path`
文件应作为或的参数提供，此时 CytoSPACE 会自动尝试从同一目录中查找对应的基因和细胞文件。`[表达式].mtx``--st-path``--scRNA-path`

请确保基因文件包含与提供ST数据中相同命名规范的基因名（例如，Ensembl基因ID、HGNC基因符号等）。
请确保基因文件包含与提供 ST 数据中相同命名规范的基因名（例如，Ensembl 基因 ID、HGNC 基因符号等）。

**从Seurat对象准备输入文件**

如果你有Seurat对象形式的数据，可以通过我们在脚本中提供的辅助函数生成适合CytoSPACE输入格式的文件。要使用这些辅助函数，首先通过以下方式导入`R``generate_cytospace_from_seurat_object.R``cytospace/Prepare_input_files``generate_cytospace_from_seurat_object.R`
如果你有 Seurat 对象形式的数据，可以通过我们在脚本中提供的辅助函数生成适合 CytoSPACE 输入格式的文件。要使用这些辅助函数，首先通过包含 `R` `generate_cytospace_from_seurat_object.R` `cytospace/Prepare_input_files` `generate_cytospace_from_seurat_object.R` 来导入它们

```bash
  source('/path/to/generate_cytospace_from_seurat_object.R')
```

用你的R字体。
用你的 R 字体。

### 来自scRNA-seq的Seurat对象
来自 scRNA-seq 的 Seurat 对象

对于从 scRNA Seurat 对象生成 CytoSPACE 输入，我们提供函数，可以称为`generate_cytospace_from_scRNA_seurat_object`
对于从 scRNA Seurat 对象生成 CytoSPACE 输入，我们提供函数，可以称为 `generate_cytospace_from_scRNA_seurat_object`

```bash
  generate_cytospace_from_scRNA_seurat_object(scRNA_Seurat_Object, dir_out='.', fout_prefix='', write_sparse=FALSE, rna_assay='RNA')
```

在你的R脚本里。 （必填） ： 输入Seurat对象
（可选，默认为工作目录） ： 输出目录的路径以存储结果
（必填） ：输入 Seurat 对象（可选，默认为无） ： 一个前缀用于添加到输出文件名中，否则生成为 和
（可选，默认为FALSE） ： 是否保存表达式数据为稀疏矩阵格式
（可选，默认为工作目录）：存储输出目录的路径（可选，默认为 ） ： 从哪个检测法取计数矩阵
请注意，必须设置为包含单元类型。
（可选，默认为无）：一个前缀，用于添加到输出文件名中，否则这些文件名通常生成为 和`scRNA_Seurat_Object``dir_out``fout_prefix``scRNA_data.txt``cell_type_labels.txt``write_sparse``rna_assay``RNA``Idents(scRNA_Seurat_Object)`
请注意，必须设置为包含单元类型。`scRNA_Seurat_Object``dir_out``fout_prefix``scRNA_data.txt``cell_type_labels.txt``write_sparse``rna_assay``RNA` `鉴定（scRNA_Seurat_Object）`

### From Spatial Seurat object
来自空间修拉对象

For producing CytoSPACE inputs from ST Seurat objects, we provide the function which may be called as`generate_cytospace_from_ST_seurat_object`
为了从 ST Seurat 对象生成 CytoSPACE 输入，我们提供了函数，可以称为 `generate_cytospace_from_ST_seurat_object`

```bash
  generate_cytospace_from_ST_seurat_object(ST_Seurat_Object, dir_out='.', fout_prefix='', write_sparse=FALSE, slice='slice1')
```

within your R script.
在你的 R 脚本里。
(required) : input Seurat object
（必填） ：输入 Seurat 对象
(optional, default is working directory) : the path to the output directory to store the results
（可选，默认为工作目录）：存储输出目录的路径
(optional, default is none) : a prefix to add to output file names, which otherwise are generated as and
（可选，默认为无）：一个前缀，用于添加到输出文件名中，否则这些文件名通常生成为 和
(optional, default is FALSE) : whether to save the expression data in sparse matrix format
（可选，默认为 FALSE）：是否保存表达式数据为稀疏矩阵格式
(optional, default is ) provides the name of your slice as stored in your Seurat object.`ST_Seurat_Object``dir_out``fout_prefix``ST_data.txt``Coordinates.txt``write_sparse``slice``slice1`
`切片` （可选，默认为`切片 1`）提供你在 Seurat 对象中存储的切片名称。

## 运行CytoSPACE
运行 CytoSPACE

展开部分
展开部分

通过 激活 conda 环境后，可以通过命令行从任意文件夹调用 CytoSPACE。关于如何运行 CytoSPACE 的更多示例可见下方“[**运行 CytoSPACE 的示例数据集**](#example-dataset-for-running-cytospace)”部分。`cytospace``conda activate cytospace``cytospace`
通过 `conda activate cytospace 激活 cytospace` 后，可以通过命令行从任何使用 `cytospace` 的文件夹调用 CytoSPACE。 关于如何运行 CytoSPACE 的更多示例可见下方“运行 CytoSPACE 的示例数据集”部分。

默认设置的典型CytoSPACE运行会是这样的：
默认设置的典型 CytoSPACE 运行会是这样的：
默认设置的典型 CytoSPACE 运行会是这样的：

```bash
cytospace \
   --scRNA-path /path/to/scRNA_geneexpression \
   --cell-type-path /path/to/scRNA_celllabels \
   --st-path /path/to/ST_geneexpression \
   --coordinates-path /path/to/ST_coordinates
```

或者参数名称更简洁：
或者参数名称更简洁：
或者参数名称更简洁：

```bash
cytospace \
   -sp /path/to/scRNA_geneexpression \
   -ctp /path/to/scRNA_celllabels \
   -stp /path/to/ST_geneexpression \
   -cp /path/to/ST_coordinates
```

或者，如果从Space Ranger输出开始，命令可能如下：
或者，如果从 Space Ranger 输出开始，命令可能如下：
或者，如果从 Space Ranger 输出开始，命令可能如下：

```bash
 cytospace \
    --scRNA-path /path/to/scRNA_geneexpression \
    --cell-type-path /path/to/scRNA_celllabels \
    --spaceranger-path /path/to/spaceranger_output.tar.gz
```

```bash
 cytospace -sp /path/to/scRNA_geneexpression \
    -ctp /path/to/scRNA_celllabels \
    -srp /path/to/spaceranger_output.tar.gz
```

有关完整使用详情及更多选项，请参见下方[**“扩展使用详情**](#extended-usage-details)”部分。
有关完整使用详情及更多选项，请参见下方“扩展使用详情”部分。
有关完整使用详情及更多选项，请参见下方“扩展使用详情”部分。

### 选择求解器

CytoSPACE 提供三种求解器选项。简而言之，如果您的系统支持 AVX2（即在安装步骤第 7 时成功安装），建议使用默认选项。使用默认求解器无需选项。使用时，将参数传递给调用。完整求解器细节请参见下方 [**CytoSPACE 求解器选项**](#cytospace-solver-options)部分。`lapjv``pip install lapjv==1.3.14``lap_CSPR``lapjv``lap_CSPR``-sm lap_CSPR``cytospace`
CytoSPACE 提供三种求解器选项。简而言之，如果你的系统支持 AVX2（即你在安装步骤 7 时成功安装了），我们建议使用默认选项。使用默认求解器无需任何选项。用的是把论点交给你的叫牌。有关完整求解器细节，请参见下方 CytoSPACE 求解器选项部分。`lapjv``pip install lapjv==1.3.14``lap_CSPR``lapjv` `lap_CSPR-sm lap_CSPR``cytospace`

**CytoSPACE的其他运行方式
CytoSPACE 的其他运行方式**

*   你可以导入 Python 里的方法或函数，然后修改或创建自己的 管道。例如：`CytoSPACE`
    你可以导入 Python 里的方法或函数，然后修改或创建自己的 管道。例如：`CytoSPACE`

```python
  from cytospace import cytospace

  for mean_cell_numbers in [5, 10, 20]:
      cytospace.main_cytospace(..., mean_cell_numbers=mean_cell_numbers)
```

## CytoSPACE 输出
CytoSPACE 输出

展开部分
展开部分

CytoSPACE默认会生成六个输出文件。
CytoSPACE 默认会生成六个输出文件。
CytoSPACE 默认会生成六个输出文件。

1.  `cell_type_assignments_by_spot.pdf`
    ST样本中细胞类型分配的热图。除了显示每个点被映射到单元格总数的图外，这些图展示了单元类型分配的空间分布。颜色条表示每个点推断出的相应单元格数量。
    ST 样本中细胞类型分配的热图。除了显示每个点被映射到单元格总数的图外，这些图展示了单元类型分配的空间分布。颜色条表示每个点推断出的相应单元格数量。
    ST 样本中细胞类型分配的热图。除了显示每个点被映射到单元格总数的图外，这些图展示了单元类型分配的空间分布。颜色条表示每个点推断出的相应单元格数量。
2.  `cell_type_assignments_by_spot_jitter.pdf`
    一个单一散点图，显示所有按位置分配的单元格。每个单元格根据其单元类型进行着色。
    一个单一散点图，显示所有按位置分配的单元格。每个单元格根据其单元类型进行着色。
    一个单一散点图，显示所有按位置分配的单元格。每个单元格根据其单元类型进行着色。
3.  `assigned_locations.csv`
    该文件将提供每个单单元格映射到ST点的指定位置。由于根据输入scRNA-seq集的大小，某些细胞可能被映射到多个位置，因此每个细胞会被分配新的细胞ID（）并在第一列中标注。第二列包含原始单元ID （）;第三列包含对应的单元类型（）;第四列包含分配的点位编号（）;第五和第六列分别包含和索引，或如初始坐标文件中提供的xy坐标，适用于相应的点。`UniqueCID``OriginalCID``CellType``SpotID``row``column``X``Y`
    该文件将提供每个单单元格映射到 ST 点的指定位置。由于部分细胞可能根据输入 scRNA-seq 集的大小被映射到多个位置，因此每个细胞会被分配新的细胞 ID（`UniqueCID`），并在第一列中给出。第二列包含原始单元代码 ID（`OriginalCID`）;第三列包含对应的单元类型（`CellType`）;第四列包含分配的点位编号（`SpotID`）;第五列和第六列分别包含对应点的`行`和`列`索引，或如初始坐标文件中提供的 X 和 `Y` 等 xy 坐标 。
4.  `assigned_expression`，一个包含 、 、 和 的目录（CytoSPACE v1.0.6+）。这表示生成后分配的基因表达，行为基因，列为细胞（ ），从原始输入的scRNA矩阵中恢复出来。表达式数据可以通过Seurat或SciPy等函数读取，用于后续分析。与 一样，请注意，可能存在被分配到多个位置的单元格，因此在表达矩阵中多次出现（在不同的 s 下）。
    为了兼容Seurat的默认参数，将基因名列在第二列，第一列用s填充。
    `assigned_expression`，一个包含 、 和 的目录（CytoSPACE v1.0.6+）`barcodes.tsv``genes.tsv``matrix.mtx``UniqueCID``assigned_locations.csv``Read10X()``io.mmread()``assigned_locations.csv``UniqueCID``Read10X()``genes.tsv``NA`
    `assigned_expression`，一个包含 `barcodes.tsv`、`genes.tsv` 和 `matrix.mtx` 的目录（CytoSPACE v1.0.6+）
5.  `cell_type_assignments_by_spot.csv`
    该文件提供了每个位置每种单元格类型的原始单元数，以及分配给该单元格的总单元数。`SpotID`
    该文件提供了每个位置每种单元格类型的原始单元数，以及分配给该单元格的总单元数。`SpotID`
    该文件通过 `SpotID` 提供了每个点的原始单元格数，以及分配给该点的总单元格数。
6.  `fractional_abundances_by_spot.csv`
    该文件给出了分配给每个点的单元类型分数丰度。`SpotID`
    该文件给出了分配给每个点的单元类型分数丰度。`SpotID`
    该文件显示了 `SpotID` 分配给每个点的单元类型比例丰度。
7.  `unassigned_locations.csv`
    该文件包含算法未分配单元格的点（位置）列表。列包括点ID（行名）、点位的行（）和列（）索引，以及每个点（），中所有点的单元格总数（），即该文件中所有点的总数为0。`row``col``Number of cells`
    该文件包含算法未分配单元格的点（位置）列表。列包括点 ID（行名）、点位的行（）和列（）索引，以及每个点（），中所有点的单元格总数（），即该文件中所有点的总数为 0。 `行``列单元`单元`数`
    该文件包含算法未分配单元格的点（位置）列表。列包括点 ID（行名）、点的行（ `行` ）和列（ `列` ）索引，以及每个点的总单元格数（单元`数` ），该文件中所有点数均为 0。
8.  `log.txt`
    该文件包含 CytoSPACE 运行参数和运行时间的日志。
    该文件包含 CytoSPACE 运行参数和运行时间的日志。

## 常见问题解答（FAQ）
常见问题解答（FAQ）

展开部分
展开部分

1.  我的ST数据集来自10x Visium以外的平台。我应该指定哪些额外参数？
    我的 ST 数据集来自 10x Visium 以外的平台。我应该指定哪些额外参数？请参阅以下相应章节中的后续说明：
    请参阅下方相应章节中的进一步说明：
    和
    。
    
2.  我的scRNA-seq数据集格式不是UMI计数矩阵。
    我的 scRNA-seq 数据集格式不是 UMI 计数矩阵。你可能需要对默认的CytoSPACE工作流程做以下两项调整。
    你可能需要对默认的 CytoSPACE 工作流程做以下两项调整。
    （1） CytoSPACE内部估计细胞类型分数的算法通常以UMI计数矩阵为参考。作为替代方案，你可以选择提供一个文件。更多信息请参见——
    。[](#advanced-options) （2） CytoSPACE（v1.0.4+）在分配前将scRNA-seq数据集降采样至每个细胞一定数量（默认为1500个），以便分配不依赖于每个细胞的总转录本数。你可以通过在CytoSPACE调用后附加标志来关闭此功能，并选择将之前下采样的scRNA-seq数据作为输入。
    （1） CytoSPACE 内部估计细胞类型分数的算法通常以 UMI 计数矩阵为参考。作为替代方案，你可以选择提供一个文件。更多信息请参见高级选项——用户提供的每种单元格类型的分数组合。`--cell-type-fraction-estimation-path``--downsample-off`
    （1） CytoSPACE 内部估计细胞类型分数的算法通常以 UMI 计数矩阵为参考。作为替代方案，你可以选择提供 `--cell-type-fraction-estimation-path` 一个文件。更多信息请参见高级选项——用户提供的每种单元格类型的分数组合。
    
3.  运行时，我会收到 / / / 错误。虽然这些错误可能由多种原因引起，但很可能CytoSPACE运行需要比可用更多的内存。
    我们提供了一种子采样程序，将ST数据集划分为更小的子集，并逐个子集评估，从而减少了内存需求。更多信息请参见 -
    。
    运行时，我会收到 `“被杀` ”/ `终止` /分`段错误` / `concurrent.futures.process.BrokenProcessPool` 错误。[](#advanced-options)如果你在单单元ST数据集中遇到这种错误，建议降低参数。
    虽然这些错误可能由多种原因引起，但很可能 CytoSPACE 运行需要比可用更多的内存。`Killed``Terminated``Segmentation fault``concurrent.futures.process.BrokenProcessPool``-noss`
    虽然这些错误可能由多种原因引起，但很可能 CytoSPACE 运行需要比可用更多的内存。
    
4.  我的输入数据非常稀少。有没有办法提供稀疏矩阵作为输入？
    我的输入数据非常稀少。有没有办法提供稀疏矩阵作为输入？自 CytoSPACE v1.1.0 起，CytoSPACE 支持稀疏矩阵作为输入，适用于单单元和空间数据集。这既适用于主细胞空间函数，也适用于get\_cellfracs\_seuratv3 R脚本用于估算每种细胞类型分数成分。虽然输入仍需遵循固定的文件命名逻辑，但我们扩大了兼容性，涵盖了最常见的命名惯例和不同版本 Cell Ranger 和 Space Ranger 稀疏矩阵输出的扩展。详情请参阅
    ——[基因表达使用稀疏矩阵](#input-files)部分。
    自 CytoSPACE v1.1.0 起，CytoSPACE 支持稀疏矩阵作为输入，适用于单单元和空间数据集。这既适用于主细胞空间函数，也适用于用于估算每种细胞类型分数组成的 get\_cellfracs\_seuratv3R 脚本。虽然输入仍需遵循固定的文件命名逻辑，但我们扩大了兼容性，涵盖了最常见的命名惯例和不同版本 Cell Ranger 和 Space Ranger 稀疏矩阵输出的扩展。详情请参阅输入文件——基因表达使用稀疏矩阵部分。
    自 CytoSPACE v1.1.0 起，CytoSPACE 支持稀疏矩阵作为输入，适用于单单元和空间数据集。这既适用于主细胞空间函数，也适用于用于估算每种细胞类型分数组成的 get\_cellfracs\_seuratv3R 脚本。虽然输入仍需遵循固定的文件命名逻辑，但我们扩大了兼容性，涵盖了最常见的命名惯例和不同版本 Cell Ranger 和 Space Ranger 稀疏矩阵输出的扩展。详情请参阅输入文件——基因表达使用稀疏矩阵部分。
    

5.  提供Space Ranger（v2.0.0+）的输出会导致错误。
    提供 Space Ranger（v2.0.0+）的输出会导致错误。自 CytoSPACE v1.1.0 起，CytoSPACE 支持 Space Ranger（v2.0.0+）的输出。详情请参阅“来自[太空游侠](#input-files)
    ”部分。
    自 CytoSPACE v1.1.0 起，CytoSPACE 支持 Space Ranger（v2.0.0+）的输出。详情请参阅“来自太空游侠的输入文件”部分。
    自 CytoSPACE v1.1.0 起，CytoSPACE 支持 Space Ranger（v2.0.0+）的输出。详情请参阅“来自太空游侠的输入文件”部分。

## 运行 CytoSPACE 的示例数据集
运行 CytoSPACE 的示例数据集
运行 CytoSPACE 的示例数据集

展开部分
展开部分

为了让用户测试 CytoSPACE，我们提供了示例运行的文件：
为了让用户测试 CytoSPACE，我们提供了示例运行的文件：

*   Wu等人（[Nature Genetics，2021](https://www.nature.com/articles/s41588-021-00911-1)）制作的HER2+乳腺癌scRNA-seq图谱，以及Visium平台（[10x Genomics](https://www.10xgenomics.com/resources/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0)）谱写的HER2+ 乳腺癌FFPE标本。默认参数的选择是基于Visium样本，适用于此处。
    Wu 等人（Nature Genetics，2021）制作的 HER2+乳腺癌 scRNA-seq 图谱，以及 Visium 平台（10x Genomics）谱写的 HER2+ 乳腺癌 FFPE 标本。默认参数的选择是基于 Visium 样本，适用于此处。

### 下载示例数据集
下载示例数据集

包含示例数据集的压缩文件可从以下链接下载：
包含示例数据集的压缩文件可从以下链接下载：

*   [乳腺癌
    乳腺癌](https://drive.google.com/file/d/1G8gK4MxCmRG4JZi588wloMsP8iZlQf_z/view?usp=share_link)

### 运行示例分析的命令：
运行示例分析的命令：

一旦示例文件下载并解压，以下命令可以在解压目录内执行：
一旦示例文件下载并解压，以下命令可以在解压目录内执行：

```bash
  cytospace -sp brca_scRNA_GEP.txt -ctp brca_scRNA_celllabels.txt -stp brca_STdata_GEP.txt -cp brca_STdata_coordinates.txt -o cytospace_results_brca -sm lap_CSPR
```

请注意，这里我们使用求解器来实现兼容性。如果你的系统支持 AVX2 内在元件，你可以运行相同的命令，但没有最终参数，直接用求解器。**CytoSPACE运行大约需要5分钟。**`lap_CSPR``lapjv`
请注意，这里我们使用求解器来实现兼容性。如果你的系统支持 AVX2 内在元件，你可以运行相同的命令，但没有最终参数，直接用求解器。CytoSPACE 运行大约需要 5 分钟。`lap_CSPR``lapjv`
请注意，这里我们使用 `lap_CSPR` 求解器来实现兼容性。如果你的系统支持 AVX2 内在元件，你可以运行相同的命令，但没有最后一个参数，改用 `lapjv` 求解器。CytoSPACE 运行大约需要 5 分钟。

### 例如，CytoSPACE输出文件，例如乳腺癌数据
例如，CytoSPACE 输出文件，例如乳腺癌数据

CytoSPACE运行的主要输出是名为的文件，该文件提供单个细胞被分配到的ST点。`assigned_locations.csv`
CytoSPACE 运行的主要输出是名为的文件，该文件提供单个细胞被分配到的 ST 点。`assigned_locations.csv`
CytoSPACE 运行的主要输出是名为`assigned_locations.csv`的文件，该文件提供单个细胞被分配到的 ST 点。

![](images/assigned_locations.png)

CytoSPACE的结果以热图形式可视化，显示每种细胞类型在ST斑点上的单细胞分布。颜色条表示每个点推断出的相应单元格数量。以下是为示例BRCA数据制作的热力图。`cell_type_assignments_by_spot.pdf`
CytoSPACE 的结果以热图形式可视化，显示每种细胞类型在 ST 斑点上的单细胞分布。颜色条表示每个点推断出的相应单元格数量。以下是为示例 BRCA 数据制作的热力图。 `cell_type_assignments_by_spot.pdf`
CytoSPACE 的结果以热图形式可视化，显示 `cell_type_assignments_by_spot.pdf` 每种细胞类型在 ST 斑点上的单细胞分布。颜色条表示每个点推断出的相应单元格数量。以下是为示例 BRCA 数据制作的热力图。

![](images/BRCA_cell_type_assignments_by_spot.png)

作为对比，请参考该ST样本的病理学家注释，由10x提供：

![](images/Visium_FFPE_Human_Breast_Cancer_Pathologist_Annotations.png)

CytoSPACE还提供了一个散点图，显示所有类型的细胞在其位置附近同时存在，保存为。每个单元格根据其单元类型进行着色。以下是为示例BRCA数据绘制的图表。`cell_type_assignments_by_spot_jitter.pdf`
CytoSPACE 还提供了一个散点图，显示所有类型的细胞在其位置附近同时存在，保存为 `cell_type_assignments_by_spot_jitter.pdf` 。每个单元格根据其单元类型进行着色。以下是为示例 BRCA 数据绘制的图表。

![](images/BRCA_cell_type_assignments_by_spot_jitter.png)

文件中提供了每个点按单元类型及总数的单元数。每种单元类型的分数丰度返回在文件中。记录 CytoSPACE 输入和运行时间的日志文件输出于 文件 。`cell_type_assignments_by_spot.csv``fractional_abundances_by_spot.csv``log.txt`
文件中提供了 `cell_type_assignments_by_spot.csv` 每个点按单元类型及总数的单元数。每种单元类型的分数丰度返回在文件 `fractional_abundances_by_spot.csv` 中。文件 `log.txt` 中输出记录 CytoSPACE 输入和运行时间的日志文件。

可通过以下链接下载预期CytoSPACE输出的压缩文件（使用CytoSPACE v1.0.0求解器）：`lap_CSPR`
可通过以下链接下载预期 CytoSPACE 输出的压缩包（使用 `lap_CSPR` 求解器支持 CytoSPACE v1.0.0）：

*   [乳腺癌治疗结果
    乳腺癌治疗结果](https://drive.google.com/file/d/1ZMA0XEl_pjC12mb8bZL8zI9yzYd2djbq/view?usp=share_link)

**模拟数据集
模拟数据集**

除了上述示例数据集外，我们生成的用于评估CytoSPACE在不同条件下鲁棒性的模拟数据集也可在下方下载。 这些数据是利用Rodriques等人（[Science，2019](https://www.science.org/doi/10.1126/science.aaw1219)）中小脑和海马切片的带注释Slide-seq数据集生成的。每个模拟数据集包含使用不同特定点分辨率（每个点5、15和30个细胞）生成的数据子目录，以及包含定义百分比基因扰动的参考单细胞数据集子目录。更多信息请参阅[论文](https://www.nature.com/articles/s41587-023-01697-9)的方法部分。`scRNA`
除了上述示例数据集外，我们生成的用于评估 CytoSPACE 在不同条件下鲁棒性的模拟数据集也可在下方下载。这些数据是利用 Rodriques 等人（Science，2019）中小脑和海马切片的带注释 Slide-seq 数据集生成的。每个模拟数据集包含使用不同斑点分辨率（每个斑点 5、15 和 30 个细胞）生成的数据子目录，以及一个 `scRNA` 子目录，包含在特定基因百分比中扰动的参考单细胞数据集。更多信息请参阅论文的方法部分。

1.  [小脑](https://drive.google.com/file/d/1qfz2T8u3HRG4qdZc9qafcO4aCvjA91Rb/view?usp=share_link)
2.  [Hippocampus](https://drive.google.com/file/d/1Jyd14n-ISc5lF65pnJWLhCCgSkpjtbsr/view?usp=share_link)

## 在遗留ST数据上运行CytoSPACE
在遗留 ST 数据上运行 CytoSPACE

展开部分
展开部分

默认情况下，CytoSPACE参数已针对标准10倍Visium空间切片进行了优化。由传统ST平台生成的数据集可以用类似命令运行，但我们建议调整以下参数：
默认情况下，CytoSPACE 参数已针对标准 10 倍 Visium 空间切片进行了优化。由传统 ST 平台生成的数据集可以用类似命令运行，但我们建议调整以下参数：

1.  `--mean_cell_numbers`，或，应设为。传统ST平台的点块大小更大，因此我们建议每个点平均映射20个单元。`-mcn``20`
    `——mean_cell_numbers`，或应设置为。传统 ST 平台的点块大小更大，因此我们建议每个点平均映射 20 个单元。`-MCN``20`
    `——mean_cell_numbers`，或者说 `-mcn`，应该设置为 `20`。传统 ST 平台的点块大小更大，因此我们建议每个点平均映射 20 个单元。
2.  `--geometry`，或应设置为。这样图函数就能将每个点塑造成正方形而非六边形。`-g``square`
    `--几何体，` 或应设置为 。这样图函数就能将每个点塑造成正方形而非六边形。`-g` `平方`
    `--几何`形状，或 `-g` 应设置为`正方`形。这样图函数就能将每个点塑造成正方形而非六边形。

与上述乳腺癌示例数据集类似，我们提供以下示例数据集：
与上述乳腺癌示例数据集类似，我们提供以下示例数据集：

*   Tirosh 等人制作的黑色素瘤 scRNA-seq 图谱（[Science，2016](https://www.science.org/doi/10.1126/science.aad0501?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)），以及由传统 ST 平台分析的黑色素瘤标本（Thrane 等，[Cancer Research，2018](https://aacrjournals.org/cancerres/article/78/20/5970/631815/Spatially-Resolved-Transcriptomics-Enables)）。
    Tirosh 等人制作的黑色素瘤 scRNA-seq 图谱（Science，2016），以及由传统 ST 平台分析的黑色素瘤标本（Thrane 等，Cancer Research，2018）。

包含该数据集的压缩文件可[在此](https://drive.google.com/file/d/1hwK_sh355chdmW50yrPJq7_W8j6HuRHh/view?usp=share_link)下载。
包含该数据集的压缩文件可在此下载。

使用以下命令运行 CytoSPACE 即可生成[此处](https://drive.google.com/file/d/1bX4SqrYzIXov_A5ivlJ8U0qD8_lXmmBf/view?usp=share_link)显示的结果。输出格式将与上述乳腺癌数据集相同。请注意，这里我们指定参数，而非使用 CytoSPACE 内部估算细胞分数算法（参见[**高级选项**](#advanced-options) - **用户提供的每种细胞类型分数**组成），因为用作参考的 scRNA-seq 图谱是通过 Smart-seq2 生成的。如果使用 CytoSPACE v1.0.4+，还应额外指定该标志。`-ctfep``--downsample-off`
使用以下命令运行 CytoSPACE 即可生成此处显示的结果。输出格式将与上述乳腺癌数据集相同。请注意，这里我们指定参数，而非使用 CytoSPACE 内部估算细胞分数算法（参见高级选项 - 用户提供的每种细胞类型分数组成），因为用作参考的 scRNA-seq 图谱是通过 Smart-seq2 生成的。如果使用 CytoSPACE v1.0.4+，还应额外指定该标志。`-ctfep``--下采样`
使用以下命令运行 CytoSPACE 即可生成此处显示的结果。输出格式将与上述乳腺癌数据集相同。请注意，这里我们指定了 `-ctfep` 参数，而非使用 CytoSPACE 内部估算细胞分数算法（参见高级选项——用户提供的每种细胞类型分数组成），因为作为参考的 scRNA-seq 图谱是使用 Smart-seq2 生成的。如果使用 CytoSPACE v1.0.4+，还应指定 `--downsample-off` 标志。

```bash
  cytospace -sp melanoma_scRNA_GEP.txt -ctp melanoma_scRNA_celllabels.txt -stp melanoma_STdata_slide1_GEP.txt -cp melanoma_STdata_slide1_coordinates.txt -ctfep melanoma_cell_fraction_estimates.txt -o cytospace_results_melanoma -mcn 20 -g square -sm lap_CSPR
```

## 在单单元ST数据上运行CytoSPACE
在单单元 ST 数据上运行 CytoSPACE

展开部分
展开部分

虽然CytoSPACE设计用于大多数Visium类型数据，其中大多数部位含有多个细胞的RNA，但CytoSPACE也可用于单细胞分辨率空间数据，如[Vizgen的MERSCOPE平台](https://vizgen.com/resources/meet-the-merscope-platform/)。我们预计该扩展有助于减少噪声并扩大ST数据中每个细胞的转录组覆盖范围，从而能够识别比单细胞分辨率ST平台更多样化的基因间空间依赖性变化。
虽然 CytoSPACE 设计用于大多数 Visium 类型数据，其中大多数部位含有多个细胞的 RNA，但 CytoSPACE 也可用于单细胞分辨率空间数据，如 Vizgen 的 MERSCOPE 平台。我们预计该扩展有助于减少噪声并扩大 ST 数据中每个细胞的转录组覆盖范围，从而能够识别比单细胞分辨率 ST 平台更多样化的基因间空间依赖性变化。

在单单元分辨率模式下，CytoSPACE将ST数据划分为更小的子集，并利用多个CPU核心将参考scRNA-seq数据的下采样版本分配给这些区域。
在单单元分辨率模式下，CytoSPACE 将 ST 数据划分为更小的子集，并利用多个 CPU 核心将参考 scRNA-seq 数据的下采样版本分配给这些区域。

要以单单元模式运行 CytoSPACE，请在您的命令中添加以下参数：
要以单单元模式运行 CytoSPACE，请在您的命令中添加以下参数：

1.  `--single-cell` (`-sc`），它告诉CytoSPACE以单单元模式运行。
    `--单细胞` （`-sc`），它告诉 CytoSPACE 以单细胞模式运行。
    `--单细胞` （`-sc`），它告诉 CytoSPACE 以单细胞模式运行。
2.  ST数据集的单元类型。请注意，对于单细胞模式，CytoSPACE不支持细胞类型分数的内部估计，用户需指定以下两种选项之一。
    ST 数据集的单元类型。请注意，对于单细胞模式，CytoSPACE 不支持细胞类型分数的内部估计，用户需指定以下两种选项之一。
    （1）（）
    **如果有单个位置的细胞类型，我们建议使用此选项。**该文件将列出每个斑点的细胞类型标签，格式与 下指定的 scRNA-seq 细胞类型标签相同。所有存在于中的细胞类型也必须存在于。
    （2）（）
    但如果用户无法访问每个具体点的单元类型，可以选择此选项。请参阅“[**高级选项**](#advanced-options)——**用户提供每种单元格类型的分数组合**”部分，了解该文件的格式化方式。
    如果有单个斑点的细胞类型，我们建议使用此选项。该文件将列出每个斑点的细胞类型标签，格式与 下指定的 scRNA-seq 细胞类型标签相同。所有存在于中的细胞类型也必须存在于。`--st-cell-type-path``-stctp``--cell-type-path``--st-cell-type-path``--cell-type-path``--cell-type-fraction-estimation-path``-ctfep`
    然而，如果用户无法访问每个具体位置的单元类型，他们可以选择使用这个选项。请参阅“高级选项——用户提供每种单元格类型的分数组合”部分，了解该文件的格式化方式。`--st-cell-类型-路径` `-stctp``--cell-type-path``--st-cell-type-path`\--cell `类型路径--细胞类型路径` `--cell-type-fraction-estimation-path` -`ctfep`
    如果有单个斑点的细胞类型，我们建议使用此选项。该文件将列出每个斑点的细胞类型标签，格式与 `--cell-type-path` 下指定的 scRNA-seq 细胞类型标签相同。-`-st-cell-type-path` 中存在的所有细胞类型也必须存在于 `--cell-type-path` 中。
3.  `--number-of-processors` (`-nop`），表示可使用的核心数量。
    `--处理器数量` （`-nop`），表示可使用的核心数。
    `--处理器数量` （`-nop`），表示可使用的核心数。
4.  `--number-of-selected-spots` (`-noss`），表示每个子集中的ST点数。我们通常推荐。`-noss 10000`
    `--被选点`数（`-noss`），表示每个子集中的 ST 点数。我们通常推荐。`-10000 号`
    `--被选点`数（`-noss`），表示每个子集中的 ST 点数。我们一般推荐 `-noss 10000`。

要运行带有单细胞分辨率空间数据的CytoSPACE：
要运行带有单细胞分辨率空间数据的 CytoSPACE：

```bash
cytospace --single-cell \
   --scRNA-path /path/to/scRNA_geneexpression \
   --cell-type-path /path/to/scRNA_celllabels \
   --st-path /path/to/ST_geneexpression \
   --coordinates-path /path/to/ST_coordinates \
   --st-cell-type-path /path/to/ST_celllabels \
   --number-of-processors NUMBER_OF_PROCESSORS \
   --number-of-selected-spots NUMBER_OF_SELECTED_SPOTS
```

或者参数名称更简洁：
或者参数名称更简洁：

```bash
cytospace -sc \
   -sp /path/to/scRNA_geneexpression \
   -ctp /path/to/scRNA_celllabels \
   -stp /path/to/ST_geneexpression \
   -cp /path/to/ST_coordinates \
   -stctp /path/to/ST_celllabels \
   -nop NUMBER_OF_PROCESSORS \
   -noss NUMBER_OF_SELECTED_SPOTS
```

一个包含单单元输入示例的压缩文件可从[Google Drive下载](https://drive.google.com/file/d/1odOcIfY3oqvLCNdXHLRaSmTraRxqnHLp/view?usp=share_link)。
一个包含单单元输入示例的压缩文件可从 Google Drive 下载。

要用这个示例数据集运行 CytoSPACE，请从解压输入位置执行以下命令，并且你的 CytoSPACE conda 环境处于激活状态：
要用这个示例数据集运行 CytoSPACE，请从解压输入位置执行以下命令，并且你的 CytoSPACE conda 环境处于激活状态：

```bash
cytospace \
   -sp HumanColonCancerPatient2_scRNA_expressions_cytospace.tsv \
   -ctp HumanColonCancerPatient2_scRNA_annotations_cytospace.tsv \
   -stp HumanColonCancerPatient2_ST_expressions_cytospace.tsv \
   -cp HumanColonCancerPatient2_ST_coordinates_cytospace.tsv \
   -stctp HumanColonCancerPatient2_ST_celltypes_cytospace.tsv \
   -o cytospace_results_crc \
   -sm lap_CSPR \
   -sc -noss 10000 -nop 2
```

Running CytoSPACE in the mode will output the assignments , the plot , and the log file . The plot generated will be a scatterplot of the cells colored by cell type, as shown below for the example dataset. The full results for the example dataset using the above command is available for download [here](https://drive.google.com/file/d/1LTTDVGAuQ4QYkyCX6WtyBNXcnZe9fxKG/view?usp=share_link).`--single-cell``assigned_locations.csv``cell_type_assignments_by_spot_single_cell.pdf``log.txt`
在 模式下运行 CytoSPACE 会输出赋值、图和日志文件。生成的图将是按细胞类型着色的散点图，示例数据集如下所示。使用上述命令完成的示例数据集完整结果可在此下载。`——单细胞` `assigned_locations.csv` `cell_type_assignments_by_spot_single_cell.pdf` `log.txt`
以`单单元`模式运行 CytoSPACE 将输出分配 `assigned_locations.csv`、图 `cell_type_assignments_by_spot_single_cell.pdf` 和日志文件 `log.txt`。生成的图将是按细胞类型着色的散点图，示例数据集如下所示。使用上述命令完成的示例数据集完整结果可在此下载。

![](images/CRC_cell_type_assignments_by_spot_single_cell.png)

## 高级选项
高级选项

虽然大多数用例推荐默认选项，但我们也提供额外的高级选项。
虽然大多数用例推荐默认选项，但我们也提供额外的高级选项。

**用户提供的每种单元类型的分数组成（-ctfep）
用户提供的每种单元类型的分数组成（-ctfep）**

为了考虑scRNA-seq和ST数据在每种细胞类型中细胞数量的差异，CytoSPACE要求ST组织中每种细胞类型的分数组成。默认情况下，CytoSPACE 会通过内部调用输入文件来生成这些信息。该脚本使用 ，作为 CytoSPACE 环境的一部分安装。我们强烈建议使用上方用于细胞类型分数估计。`get_cellfracs_seuratv3.R``Seurat v3``Seurat v3``Seurat v4`
为了考虑 scRNA-seq 和 ST 数据在每种细胞类型中细胞数量的差异，CytoSPACE 要求 ST 组织中每种细胞类型的分数组成。默认情况下，CytoSPACE 会通过内部调用输入文件来生成这些信息。该脚本使用 ，作为 CytoSPACE 环境的一部分安装。我们强烈建议使用上方用于细胞类型分数估计。`get_cellfracs_seuratv3.R``Seurat v3``Seurat v3``Seurat v4`
为了考虑 scRNA-seq 和 ST 数据在每种细胞类型中细胞数量的差异，CytoSPACE 要求 ST 组织中每种细胞类型的分数组成。默认情况下，CytoSPACE 会通过内部调用`get_cellfracs_seuratv3.R`输入文件来生成这些信息。该脚本使用 `Seurat v3`，作为 CytoSPACE 环境的一部分安装。我们强烈建议在细胞类型分数估计中使用`修拉 v3` 而非`修拉 v4`。

虽然我们提供的脚本使用[了空间修拉](https://satijalab.org/seurat/articles/spatial_vignette.html)，但还有多种其他方法可用，如[cell2location](https://www.sanger.ac.uk/tool/cell2location/)、[SPOTlight](https://github.com/MarcElosua/SPOTlight)或[CIBERSORTx](https://cibersortx.stanford.edu/)。
虽然我们提供的脚本使用了空间修拉，但还有多种其他方法可用，如 cell2location、SPOTlight 或 CIBERSORTx。

用户可以选择跳过 CytoSPACE 的内部算法，而是提供自己的 ST 数据集中估计的单元类型组成文件，该文件由 （） 标志指定。特别建议，如果参考scRNA-seq数据集来自非基于UMI计数的技术，如Smart-seq，应提供单独文件。`--cell-type-fraction-estimation-path``-ctfep``-ctfep`
用户可以选择跳过 CytoSPACE 的内部算法，而是提供自己的 ST 数据集中估计的单元类型组成文件，该文件由 （） 标志指定。特别建议，如果参考 scRNA-seq 数据集来自非基于 UMI 计数的技术，如 Smart-seq，应提供单独文件。 `--cell-type-fraction-estimation-path` `-CTFEP-CTFEP`
用户可以选择跳过 CytoSPACE 的内部算法，而是提供自己的文件，用于用 `--cell-type-fraction-estimation-path` （`-ctfep`） 标志指定 ST 数据集的估计单元类型组成。特别建议，如果参考 scRNA-seq 数据集来自非基于 UMI 计数的技术，如 Smart-seq，应单独提供`一个-ctfep` 文件。

所提供文件必须是一个由两行组成的表格，行名为一行，第一行包含单元格类型，第二行包含每种单元格类型的单元分数，比例介于0到1之间。**请确保第一行的单元类型标签与单元格类型标签文件中的标签一致，且单元格类型分数之和为1。两行的行名称必须同时存在。**
所提供的文件必须是一个由两行组成的表格，行名为，第一行包含单元类型，第二行包含每种单元格类型的单元分数，以0到1之间的比例表示。请确保第一行的单元类型标签与单元类型标签文件中的标签一致，且单元格类型分数之和为1。两行的行名必须同时存在。

![](images/cell_type_fractions_file.png)

**用户提供每个点的单元格数估计（-ncpsp）
用户提供每个点的单元格数估计（-ncpsp）**

用户无需使用CytoSPACE内部机制估算每个点的细胞数，而是可以在一个带有标题的两列文件中提供自己的估计（例如图像分割），第一列包含点ID，第二列显示每个点预测的单元格数：
用户无需使用 CytoSPACE 内部机制估算每个点的细胞数，而是可以在一个带有标题的两列文件中提供自己的估计（例如图像分割），第一列包含点 ID，第二列显示每个点预测的单元格数：

![](images/n_cells_per_spot.PNG)

要用这个选项运行 CytoSPACE，只需传递标志或后跟文件位置。`-ncpsp``--n-cells-per-spot-path`
要用这个选项运行 CytoSPACE，只需传递标志或后跟文件位置。`-ncpsp``--n-单元-每个点路径`
要运行带有此选项的 CytoSPACE，只需传递 `-ncpsp` 或 `--n-cells-per spot-path` 标志，再输入文件位置。

**并行化的点子采样（-sss）
并行化的点子采样（-sss）**

运行 CytoSPACE 所需的内存和运行时间可能因位置数量而异。为了让CytoSPACE在不同条件下运行，我们提供了将ST样本中估计细胞数划分为更小块的选项，然后通过多个CPU核心分配类似的下采样参考scRNA-seq数据。请注意，该选项仅适用于非单单元ST数据集（10倍Visium和遗留ST）;在单单元ST数据上运行CytoSPACE的用户应修改参数以实现同样效果。`-noss`
运行 CytoSPACE 所需的内存和运行时间可能因位置数量而异。为了让 CytoSPACE 在不同条件下运行，我们提供了将 ST 样本中估计细胞数划分为更小块的选项，然后通过多个 CPU 核心分配类似的下采样参考 scRNA-seq 数据。请注意，该选项仅适用于非单单元 ST 数据集（10 倍 Visium 和遗留 ST）;在单单元 ST 数据上运行 CytoSPACE 的用户应修改参数以实现同样效果。`-诺斯`

用户可以通过指定（）标志、每个分区所需的子采样单元数（或）以及所需核心数（或）来使用此选项。请注意，正确的组合高度依赖用户的系统配置，可能需要尝试不同的组合以确定哪种能成功运行。理想情况下，使用最高值且不超出内存会是最有效的。我们注意到，在8GB内存环境下，BRCA数据集的效果很好。`--sampling-sub-spots``-sss``--number-of-selected-sub-spots``-nosss``--number-of-processors``-nop``-nosss``-nop``-nosss``-sss -nosss 3000 -nop 2`

以下示例命令将对示例乳腺癌数据集运行 CytoSPACE，使用两个核心一次分配 5000 个细胞 scRNA-seq 数据：

```bash
  cytospace \
    -sp brca_scRNA_GEP.txt \
    -ctp brca_scRNA_celllabels.txt \
    -stp brca_STdata_GEP.txt \
    -cp brca_STdata_coordinates.txt \
    -o cytospace_results_brca \
    -sm lap_CSPR \
    -sss -nosss 5000 -nop 2
```
**备用距离度量（-dm）**

默认情况下，CytoSPACE使用皮尔逊相关法比较细胞和斑点转录组。用户可以选择使用斯皮尔曼相关性或欧几里得距离，分别通过函数调用 或 传递。`-dm Spearman_correlation``-dm Euclidean`

**设置一个新的随机种子（-se）**

虽然CytoSPACE算法大多是确定性的，但对待映射单元格的初始采样步骤是随机进行的。为了提供另一种随机种子，产生不同的随机样本，用户可以先传递，然后通过函数调用所需的整数种子。CytoSPACE的默认随机种子是1。`-se`

**采样（-sam）的替代处理**

CytoSPACE首先创建一个与ST数据预期匹配的细胞池。默认情况下，这通过重新采样单个细胞来实现，以获得组织中估计的整体细胞类型分数和总细胞数。我们建议在所有真实数据分析中都采用此默认设置运行CytoSPACE。不过，我们还提供了一种额外的选项，通过从每种细胞类型内的基因计数分布中抽样生成新的“占位符”细胞，并利用这一方法确保映射细胞的唯一性，以便对模拟数据进行基准测试。要用这种替代模式运行 CytoSPACE，用户可以通过函数调用传递。在占位符模式下运行时，目录不会被生成;相反，新生成的“占位符”细胞的基因表达将作为输出的一部分保存在 。`-sam place_holders``assigned_expression``new_scRNA.csv`

**方法扩展：映射质量**

虽然CytoSPACE作为线性赋值问题的表述保证了成本函数的最优解，但没有基于概率的概率框架来估计映射不确定性。一种可能性是在映射后确定某一单元格类型是否属于某个特定点——即该点是否至少包含一个相同单元格类型。值得注意的是，这并未区分同类型细胞之间的拟合质量。因此，该协议虽然不完整，但提供了某种映射质量的辅助脚本，通过支持向量机生成并训练由输入的scRNA-seq数据生成的伪块数来实现这一点。该脚本 输入路径为 ST 数据集计数矩阵文件、scRNA-seq 计数矩阵文件和 CytoSPACE 输出文件 ，返回附加输出文件，置信度为 。完成 CytoSPACE 运行后执行该脚本的命令如下：`uncertainty_quantification.R``assigned_locations.csv``assigned_locationswConfidenceScores.csv`

```bash
Rscript uncertainty_quantification.R /path/to/ST_geneexpression /path/to/scRNA_geneexpression /path/to/assigned_locations.csv
```

在解读置信度分数时，我们建议将阈值设为0.1，分数越高，表示对某点至少包含同一细胞类型细胞的信心越高。

请注意，这需要与提供文件中包含的conda环境依赖分开。该脚本应在安装以下 R 包的独立环境中运行：（必须是 v4;测试过 v4.0.1）、（测试过 v1.14.0）和（测试过 v1.7.8）。`uncertainty_quantification.R``environment.yml``cytospace``Seurat``data.table``e1071`

## 扩展使用详情

展开部分

```
usage: cytospace [-h] -sp SCRNA_PATH -ctp CELL_TYPE_PATH [-stp ST_PATH] [-cp COORDINATES_PATH] [-srp SPACERANGER_PATH]
                 [-stctp ST_CELL_TYPE_PATH] [-ctfep CELL_TYPE_FRACTION_ESTIMATION_PATH] [-ncpsp N_CELLS_PER_SPOT_PATH]
                 [-o OUTPUT_FOLDER] [-op OUTPUT_PREFIX] [-mcn MEAN_CELL_NUMBERS] [--downsample-off]
                 [-smtpc SCRNA_MAX_TRANSCRIPTS_PER_CELL] [-sc] [-noss NUMBER_OF_SELECTED_SPOTS] [-sss]
                 [-nosss NUMBER_OF_SELECTED_SUB_SPOTS] [-nop NUMBER_OF_PROCESSORS] [-sm {lapjv,lapjv_compat,lap_CSPR}]
                 [-dm {Pearson_correlation,Spearman_correlation,Euclidean}] [-sam {duplicates,place_holders}]
                 [-se SEED] [-p] [-g GEOMETRY] [-nc NUM_COLUMN] [-mp MAX_NUM_CELLS_PLOT]

CytoSPACE is a computational strategy for assigning single-cell transcriptomes to in situ spatial transcriptomics (ST)
data. Our method solves single cell/spot assignment by minimizing a correlation-based cost function through a linear
programming-based optimization routine.

optional arguments:
  -h, --help            show this help message and exit
  -stp ST_PATH, --st-path ST_PATH
                        Path to spatial transcriptomics data (expressions)
  -cp COORDINATES_PATH, --coordinates-path COORDINATES_PATH
                        Path to transcriptomics data (coordinates)
  -srp SPACERANGER_PATH, --spaceranger-path SPACERANGER_PATH
                        Path to SpaceRanger tar.gz data file
  -stctp ST_CELL_TYPE_PATH, --st-cell-type-path ST_CELL_TYPE_PATH
                        Path to ST cell type file (recommended for single-cell ST)
  -ctfep CELL_TYPE_FRACTION_ESTIMATION_PATH, --cell-type-fraction-estimation-path CELL_TYPE_FRACTION_ESTIMATION_PATH
                        Path to ST cell type fraction file (recommended for bulk ST)
  -ncpsp N_CELLS_PER_SPOT_PATH, --n-cells-per-spot-path N_CELLS_PER_SPOT_PATH
                        Path to number of cells per ST spot file
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Relative path to the output folder
  -op OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Prefix of results stored in the 'output_folder'
  -mcn MEAN_CELL_NUMBERS, --mean-cell-numbers MEAN_CELL_NUMBERS
                        Mean number of cells per spot, default 5 (appropriate for Visium). If analyzing legacy spatial
                        transcriptomics data, set to 20
  --downsample-off      Turn off downsampling for scRNA-seq data
  -smtpc SCRNA_MAX_TRANSCRIPTS_PER_CELL, --scRNA_max_transcripts_per_cell SCRNA_MAX_TRANSCRIPTS_PER_CELL
                        Number of transcripts per cell to downsample scRNA-seq dataset to. This allows for assignments
                        that are not dependent on the overall expression level
  -sc, --single-cell    Use single-cell spatial approach if specified
  -noss NUMBER_OF_SELECTED_SPOTS, --number-of-selected-spots NUMBER_OF_SELECTED_SPOTS
                        Number of selected spots from ST data used in each iteration
  -sss, --sampling-sub-spots
                        Sample subspots to limit the number of mapped cells if specified
  -nosss NUMBER_OF_SELECTED_SUB_SPOTS, --number-of-selected-sub-spots NUMBER_OF_SELECTED_SUB_SPOTS
                        Number of selected subspots from ST data to limit the number of mapped cells
  -nop NUMBER_OF_PROCESSORS, --number-of-processors NUMBER_OF_PROCESSORS
                        Number of processors used for the analysis
  -sm {lapjv,lapjv_compat,lap_CSPR}, --solver-method {lapjv,lapjv_compat,lap_CSPR}
                        Which solver to use for the linear assignment problem, default 'lapjv'
  -dm {Pearson_correlation,Spearman_correlation,Euclidean}, --distance-metric {Pearson_correlation,Spearman_correlation,Euclidean}
                        Which distance metric to use for the cost matrix, default 'Pearson_correlation'
  -sam {duplicates,place_holders}, --sampling-method {duplicates,place_holders}
                        Which underlying method to use for dealing with duplicated cells, default 'duplicates'
  -se SEED, --seed SEED
                        Set seed for random generators, default 1
  -p, --plot-off        Turn create plots on/off
  -g GEOMETRY, --geometry GEOMETRY
                        ST geometry, either 'honeycomb' or 'square' accepted
  -nc NUM_COLUMN, --num-column NUM_COLUMN
                        Number of columns in figure
  -mp MAX_NUM_CELLS_PLOT, --max-num-cells-plot MAX_NUM_CELLS_PLOT
                        Maximum number of cells to plot in single-cell visualization

Required arguments:
  -sp SCRNA_PATH, --scRNA-path SCRNA_PATH
                        Path to scRNA-Seq data
  -ctp CELL_TYPE_PATH, --cell-type-path CELL_TYPE_PATH
                        Path to cell type labels
```

You can see this list of variables and default values for running CytoSPACE from the commmand line as well at any time by calling along with the or flag, i.e., .`cytospace``-h``--help``cytospace -h` 

## CytoSPACE 求解器选项

展开部分

1.  `lapjv` **（大多数系统推荐）**默认情况下，CytoSPACE 调用 package 中的求解器。该求解器是Jonker-Volgenant最短增强路径分配算法的快速实现，并在给定我们[论文](https://www.nature.com/articles/s41587-023-01697-9)定义的目标函数时返回一个全局最优解。然而，如上所述，该软件包并非所有系统都支持，因为它通过使用 AVX2 指令实现了加速。该求解器默认被选中，可以通过传递参数或 来显式指定。`lapjv``lapjv``--solver-method lapjv``-sm lapjv``cytospace`
2.  `lap_CSPR` **（推荐用于不支持 `lapjv` 的系统）**第二个求解器选项是包里的方法。该求解器采用了与第一和第三个选项不同的方法，是一种称为成本尺度推送重标法的赋值算法。该算法将赋值成本近似于整数，因此会失去一些数值精度。因此，虽然在给出论文中定义的目标函数后，经过**近似**后返回的全局最优解，但结果与前两种方法相似但通常不完全相同。该求解器运行时间与第一个选项相似，是不支持该软件包的系统的一个好选择。该求解器可通过传递参数或 来选择。`linear_assignment``ortools``lapjv``--solver-method lap_CSPR``-sm lap_CSPR``cytospace`
3.  `lapjv_compat`第三个求解器选项实现了包中的求解器。与第一个选项一样，该求解器也实现了Jonker-Volgenant最短增广路径赋值算法，在论文定义的目标函数下返回相同的全局最优解。此外，该求解器广泛支持，应适用于所有标准作系统。然而，运行时间是第一个求解器选项（包中的求解器）的3-4倍，因此我们仅推荐不支持第一个选项的系统使用。该求解器可通过传递参数或 来选择。`lapjv``lap``lapjv``lapjv``lapjv``--solver-method lapjv_compat``-sm lapjv_compat``cytospace`

## 本地安装更新

展开部分

要在更新本 GitHub 仓库后更新本地安装的 CytoSPACE，请进入你的目录并执行以下命令：`cytospace`

```bash
  git pull
  conda env update --name cytospace --file environment.yml
  conda activate cytospace
  pip install .
```

如果你已经对 CytoSPACE 源代码的本地版本进行了更新，你应该执行

```bash
  pip install .
```

然后就跑了。

## 作者

CytoSPACE由[Newman实验室](https://anlab.stanford.edu/)开发
CytoSPACE 由 [Newman 实验室](https://anlab.stanford.edu/)开发

*   米拉德·R·瓦希德（miladrv）
*   艾琳·L·布朗（艾林布朗）
*   克洛伊·B·斯汀（cbsteen）
*   张武冰（Wubing Zhang）
*   贤素全（hsjeon-k）
*   亚伦·M·纽曼（aaronmnewman）

## 联系

如有任何疑问，请联系CytoSPACE团队 [cytospaceteam@gmail.com](mailto:cytospaceteam@gmail.com)。

## 许可证

请参阅[许可证](LICENSE)文件。

## 引文

如果您使用CytoSPACE，请注明：

*单细胞和空间转录组的高分辨率比对（CytoSPACE*）（Nature Biotechnology 2023） Milad R. Vahid\*， Erin L. Brown\*， Chloé B. Steen\*， Wubing Zhang， Hyun so Jeon， Minji Kang， Andrew J. Gentles， Aaron M. Newman.
*单细胞和空间转录组的高分辨率比对（CytoSPACE）（*Nature Biotechnology 2023） Milad R. Vahid\*， Erin L. Brown\*， Chloé B. Steen\*， Wubing Zhang， Hyun so Jeon， Minji Kang， Andrew J. Gentles， Aaron M. Newman.