# Withers+ 2026 (arXiv:2606.06585) — 致密REG 交叉验证备忘录

**日期**: 2026-06-09  
**论文**: Sunna Withers, Adam Muzzin 等 19 人, "A Population of Red Galaxies with Very Strong Emission Lines at z>5: 'Classic' LRDs, Dusty Star-Forming Galaxies, and a Missing Population of LRDs" (已投稿 ApJ)  
**与 ApJL AAS77670 的关系**: 同期投稿（6月4日 vs 6月5日），不构成审稿风险，但为未来交叉验证提供独立样本

---

## 核心发现

Withers+ 用 NIRCam 中波段（非宽带颜色）在 CANUCS/Technicolor/JUMPS 巡天（58.2 arcmin²）中发现 26 个 REG（红色发射线星系），分为三类：

| 子类 | 数量 | 性质 |
|------|------|------|
| Classic LRDs | 8 | 同时满足经典宽带选择标准 |
| Extended REGs | 9 | 可分辨（Reff 1-3 kpc），≈ DSFG |
| **Compact REGs** | **9** | 不可分辨（点源），但经典 LRD 标准全未选中 |

**遗漏原因**：连续谱极暗（F444W>28 mag）+ [OIII]+Hβ 极强（EW 最高 ~2300 Å）→ V 形 SED 出不来 → 宽带颜色选不中。

论文还有一条普遍规律：ELG 连续谱颜色与发射线强度反相关（Spearman ρ=-0.34, p~10⁻²⁵），REG 是这个关系上红偏 2σ 的离群体。

---

## 与 ApJL AAS77670 的 F444W 表面亮度对比

| | ApJL 交叉匹配 38 源 | Withers+ 致密 REG 9 源 |
|---|---|---|
| r_eff | 52–197 pc | **< 509–1020 pc**（PSF 上限） |
| r_eff_arcsec | 0.010–0.033″ | **< 0.081″**（未分辨） |
| F444W mag | 30.9–35.0 | 27.5–28.7 |
| log M* | 8.8–10.6 | 8.0–8.4 |
| **SB_F444W** | 22.1–28.2（**中位 24.5**） | **> 23.3–24.5（下限！中位 24.0）** |

### 关键结论

1. **Withers+ 致密REG 的 F444W 表面亮度下限（中位 24.0）已比 ApJL 样本中位数（24.5）亮 0.5 mag**。真实值更亮。

2. 它们在 SB 维度上与 ApJL 样本高端重叠，是"亮的漏网之鱼"——不是因为暗才漏，是因为**发射线污染把颜色洗平了**。

3. **恒星质量低一个量级**（log M* ≈ 8 vs 9-10），暗示中心源主导了 F444W 光度但总恒星质量不高。

4. **缺 BD 数据**——没有巴尔末减缩就无法验证 Σ-BD 相关性。这是未来的可检验预测：如果拿到光谱，它们应落在相关线的高 Σ、高 BD 端。

---

## 战略意义

- **对 ApJL 审稿无影响**（同期投稿，审稿人不会引用未发表论文）
- **对未来 ApJ Main Journal（260 源 + 1620 对照组）是弹药**：可专门辟一节讨论样本完备性，引用 Withers+ 作为独立验证——"中波段完备样本预测相关性更强"
- **可检验预测**：如果致密REG 是 LRD，其 BD 应落在 Σ-BD 相关线的高端延伸

---

## 文件清单

`11_Withers_MissingLRDs_CrossCheck/`

| 文件 | 内容 |
|------|------|
| `withers_compact_regs.csv` | 9 致密REG 完整参数（含 SB 下限） |
| `apjl_crossmatch_with_sb.csv` | ApJL 38 交叉源 + F444W SB |
| `kokorev_full_sb.csv` | Kokorev 全 260 源 + F444W SB |
| `comparison_stats.json` | 对比统计摘要 |
| `plot_data.json` | 画图数据（供后续可视化） |
| `withers_analysis.py` | 完整分析脚本（可复现） |
