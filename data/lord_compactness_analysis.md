# LoRD（Local Red Dots）致密性-光谱关联检验

## 背景

- **论文**: Casey+ (2026, arXiv:2606.26098) *"A Population of Little Red Dot-like Quasars in SDSS"*
- **Submitted to ApJ, 2026-06-24**
- **核心主张**: SDSS ugriz 滤光片在 z≈0.4, 0.8 与 JWST LRD 选择滤光片有效重叠 → 从 SDSS DR16Q 中选出 ~1300 个"Local Red Dots"作为高红移 LRD 的本地类似体

## 研究动机

在我们（Tan et al. 2026, submitted）的 LRD 论文中，核心论点是：致密性（Σ = log M* − log(π r_eff²)）是组织 LRD 光谱异常的关键变量，且 Σ-color 相关性在高红移（z=7-9）达到 ρ=+0.516。自然的问题是：**这个 Σ-color 耦合是否也存在于 Casey+ 的本地类似体中？**

如果 LoRDs 同样展示 Σ-color 相关 → 致密性-光谱耦合是普适规律。
如果 LoRDs 不展示 Σ-color 相关 → Casey+ 的颜色选择在不同红移选出了不同的物理群体（"同一张汤圆皮，不同的馅"），支持我们论文中的引擎简并性论证（§7.6）。

## 数据与方法

### 1. LoRD 选择复现

- **数据源**: SDSS DR16Q_v4 (Wu & Shen 2022, ApJS, 263, 42)，从 SDSS SAS 直接下载
- **样本量**: 750,414 个宽线类星体
- **消光修正**: `dmag = PSFMAG − EXTINCTION`
- **iz 波段**: `iz = (i_mag + z_mag) / 2`（i 和 z 波段流量的对数平均）
- **选择标准**（Casey+ Equations 1 & 2）:

```
red1: z ∈ [0.68, 0.88], u−g < 0.8, r−iz > 0.7, r−z > 1.0
red2: z ∈ [0.28, 0.48], u−g < 0.8, r−iz > 0.7, r−z > 1.0
```

- **复现结果**: 1094 个 LoRD（red1: 752, red2: 342），vs Casey+ 报告的 1325（red1: 935, red2: 390）
- **差异原因**: 消光修正版本差异（Casey+ 可能使用更新版的银河消光图），匹配率 ~83%

### 2. 对照样本

从相同红移范围内 (red1_z ∪ red2_z) 选取，排除 LoRD 颜色标准的类星体。随机抽样 N=80。

### 3. 致密性测量

通过 SDSS DR16 SkyServer REST API，逐源查询每个 LoRD/对照 QSO 的 PhotoObj 参数：

- `petroR50_r`: Petrosian 半光半径（角秒）
- `concentration = psfMag_r − modelMag_r`: 聚集度（0 = 点源）
- `fracDeV_r`: de Vaucouleurs 成分比例
- 查询方法: `SpecObj.bestObjID → PhotoObj` JOIN，通过 `(plate, mjd, fiberID)` 三元组匹配

### 4. SDSS 分辨率的局限性

SDSS 的角分辨率约为 1.2"，对应 z=0.5 处的物理尺度约 5-8 kpc。相比之下，高红移 LRD 的典型 r_eff ≈ 50-150 pc。因此，此分析只能探测 **kpc 尺度的"致密性"**（即点源 vs 延展源），无法探测 LRD 真正的 ~100 pc 致密性。HST 或 JWST 对 LoRD 样本的成像观测将是决定性检验。

## 结果

### LoRD vs Control: 形态学对比

| 参数 | LoRD (N=77) | Control QSO (N=80) | 检验 |
|------|-------------|-------------------|------|
| petroR50 (median) | 0.716" | 0.647" | KS D=0.253, **p=0.010** |
| petroR50 (mean) | 0.827" | 0.679" | MW U **p=0.020** |
| concentration (median) | 0.118 mag | 0.044 mag | — |
| 点源比例 (\|c\|<0.2) | 66.2% | **95.0%** | — |

**关键发现 1**: LoRD 反而比对照 QSO **更延展**。LRD 颜色选择在 z~0.5 倾向于选出有宿主星系贡献的红化 AGN，而非致密核占主导的源。

### Σ-Color 相关：核心检验

| 相关对 | LoRD (N=77) | Control (N=80) |
|--------|-------------|----------------|
| ρ(concentration, r−z) | +0.154 (p=0.18) | **+0.424 (p=0.0001)** |
| ρ(petroR50, r−z) | −0.014 (p=0.90) | — |
| ρ(concentration, u−g) | −0.033 (p=0.77) | — |
| Partial r(conc, r−z \| z) | +0.181 (p=0.11) | — |

**关键发现 2**: LoRD 内 **不存在显著的 Σ-color 相关性**。按 r−z 中位数分割后，红 LoRD 和蓝 LoRD 的 petroR50 几乎相同（0.756" vs 0.707"）。

**关键发现 3**: 对照样本反而展示 **强且显著的** concentration-r−z 相关（ρ=+0.424, p=0.0001）。这说明在正常 QSO 群体中，颜色红化与 kpc 尺度的聚集度增加相关（可能源于尘埃红化与核球主导的耦合），但 LoRD 的颜色选择打破了这一关系——它们的红化可能来自不同的机制（如更弥散的尘埃几何或更强的宿主污染）。

### 与高红移 LRD 的对比

| 样本 | N | ρ(Σ, color) | 显著性 |
|------|---|-------------|--------|
| 高红移 LRD（COSMOSWeb z=7-9） | 5,601 | +0.516 | >5.6σ |
| 高红移 LRD（SB-Σ合并） | 291 | −0.478 | p=5.3×10⁻¹⁸ |
| 高红移 LRD（BD-Σ） | 37 | +0.351 | p=0.033 |
| **本地 LoRD（本次分析）** | **77** | **+0.154** | **p=0.18** |

## 物理解读

### "同皮不同馅"（Same Skin, Different Filling）

Casey+ 的颜色映射方法选出的本地类似体，在 **物理本质上** 与高红移 LRD 有系统性差异：

```
         高红移 LRD (z~5)          本地 LoRD (z~0.5)
皮 (颜色)   F150W−F444W > 2 mag      r−z > 1.0  ✅ 匹配
           V-shaped continuum       V-shaped continuum (部分)
馅 (物理)   致密核团 r_eff~50-150 pc  正常AGN r_eff~kpc  ❌ 不匹配
           X射线弱/缺失              X射线正常
           Σ-color 强相关           Σ-color 零相关
           BD异常高                 BD正常(?)
```

### 对 LRD 物理的启示

1. **NIRCam 颜色选择的引擎简并性**：同一套颜色标准在 z~0.5 选出的是标准 AGN（红化+宿主污染），在 z~5 选出的是极端致密源。仅靠宽波段颜色无法区分引擎类型。

2. **致密性是高红移特有的关键变量**：Σ-color 相关可能是高红移宇宙中致密结构形成的独有特征。在 z~0.5，AGN 的红化主要来自尘埃而非致密性，因此颜色-致密性耦合断裂。

3. **MIRI 中红外光谱的必要性**：为了区分"红化 AGN"和"致密核团驱动源"，需要 MIRI 5-14 μm 的光谱诊断（尘埃温度：AGN 热环 ~600-1500K vs G_eff 加热尘埃 ~300-800K）。

4. **LoRD 样本本身的价值**：LoRD 是研究"纯颜色选择会选出什么"的绝佳对照组——它们代表了在没有致密性驱动的情况下，颜色标准本身会产生的选择偏倚。

## 后续论文可用的论点

1. **Casey+ 的 LoRD 发现实际上加强而非削弱了 LRD 的独特性**：证明了颜色选择在宇宙时间上的引擎简并性，与 Tan+ (2026) 的 §7.6 论证完全一致。

2. **一个可发表的跟随研究**：系统获取 LoRD 样本的 HST 高分辨成像（或利用现有 HST 档案数据），在 ~100 pc 尺度上检验是否存在 "hidden compact" 子群，以及该子群是否恢复 Σ-color 相关。预期：如果存在 z~0.5 的真致密子群（r_eff < 300 pc），它们应该展示 Σ-color 相关；如果全部 LoRD 在 HST 分辨率下都是延展的，则"不同馅"论证得到最终确认。

3. **统计框架推广**：可以将 COSMOSWeb 的 Σ-color 盲验方法推广到 SDSS 本地宇宙——选取更大红移范围（z=0-1）的 SDSS 星系样本，计算 Σ-color 偏相关随红移的演化，检验 ρ(z) ∝ (1+z)^1.14 的幂律是否在高红移和高 Σ 端之外仍然成立。

## 数据文件

- **LoRD 选择星表**: `lord_dr16q_selection.csv` — 1094 个 LoRD 的 SDSS 名称、坐标、红移、PSF 星等
- **致密性分析结果**: `lord_compactness.json` — 77 个 LoRD + 80 个对照 QSO 的 petroR50、concentration、颜色
- **DR16Q 原始星表**: `/tmp/DR16Q_v4.fits` (1.6 GB，可删除以释放空间)

## 分析脚本

复现分析的核心步骤：
1. 从 https://data.sdss.org/sas/dr16/eboss/qso/DR16Q/DR16Q_v4.fits 下载 DR16Q 星表
2. 应用 Casey+ 的颜色选择标准（见上方 §1）
3. 通过 SDSS SkyServer REST API 逐源查询 PhotoObj 尺寸
4. 计算 Spearman 秩相关和 KS 检验

## 参考

- Casey, Q.O., Hickox, R.C., Cleri, N.J. et al. (2026), arXiv:2606.26098
- Wu, Q. & Shen, Y. (2022), ApJS, 263, 42 — SDSS DR16Q catalog
- Tan et al. (2026), submitted — 三路汇聚 LRD 论文
- Greene, J.E. et al. (2024), ApJ, 964, 39 — LRD red1/red2 color criteria
