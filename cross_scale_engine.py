"""
跨尺度级联分析引擎 - Cross-Scale Cascade Analysis Engine
提供分子 → 细胞/组织 → 群体三层级联分析框架
包含: 参数传递、对比可视化、多疾病对比、基因追踪
"""

import numpy as np
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field
import math


# ==================== 数据结构 ====================

@dataclass
class ScaleResult:
    """单尺度分析结果"""
    scale: str                          # "molecular" | "cellular" | "population"
    scale_cn: str                       # 中文标签
    model_id: str                       # 使用的模型ID
    model_name: str                     # 模型中文名
    summary: Dict[str, Any] = field(default_factory=dict)
    detail: Any = None
    figure: Optional[go.Figure] = None

@dataclass
class CascadeReport:
    """跨尺度级联报告"""
    disease: str
    results: Dict[str, ScaleResult] = field(default_factory=dict)
    cross_scale_insights: List[str] = field(default_factory=list)
    architecture_html: str = ""
    param_transfer: Dict[str, Any] = field(default_factory=dict)  # 参数传递记录


# ==================== 引擎主类 ====================

class CrossScaleEngine:
    """
    跨尺度级联分析引擎

    核心功能:
        1. 跨尺度参数传递  — 分子层指标自动推导细胞层/群体层参数
        2. 跨尺度对比可视化 — 三层并排网络 + 雷达图
        3. 多疾病对比分析   — 多疾病在各尺度的指标对比
        4. 基因追踪         — 追踪单个基因跨三层的角色变化
    """

    def __init__(self, db=None, social_sim=None, tcga_sim=None):
        self.db = db
        self.social_sim = social_sim
        self.tcga_sim = tcga_sim

    # ==================================================================
    # Layer 1: 分子尺度
    # ==================================================================

    def analyze_molecular(self, disease_name: str,
                          network_type: str = "interaction") -> ScaleResult:
        """分子尺度分析：基因互作/调控网络 + 网络拓扑"""
        result = ScaleResult(
            scale="molecular", scale_cn="分子尺度",
            model_id="gene_interaction" if network_type == "interaction" else "gene_regulation",
            model_name="基因互作网络模型" if network_type == "interaction" else "基因调控网络模型",
        )

        if self.db is None:
            result.summary = {"error": "数据库未初始化"}
            return result

        disease = self.db.diseases.get(disease_name)
        if disease is None:
            result.summary = {"error": f"未找到疾病: {disease_name}"}
            return result

        G = nx.Graph() if network_type == "interaction" else nx.DiGraph()
        for gene in disease.genes:
            G.add_node(gene)

        relations = disease.interactions if network_type == "interaction" else disease.regulations
        for rel in relations:
            if network_type == "interaction":
                G.add_edge(rel.gene1, rel.gene2)
            else:
                G.add_edge(rel.regulator, rel.target)

        n_nodes = G.number_of_nodes()
        n_edges = G.number_of_edges()
        density = nx.density(G) if n_nodes > 1 else 0

        if network_type == "interaction" and n_nodes > 0:
            clustering = nx.average_clustering(G) if n_nodes > 2 else 0
            degrees = dict(G.degree())
            avg_degree = sum(degrees.values()) / max(n_nodes, 1)
            hub_genes = sorted(degrees, key=degrees.get, reverse=True)[:10]
        else:
            clustering = 0
            avg_degree = 0
            hub_genes = []
            degrees = dict(G.degree()) if n_nodes > 0 else {}

        result.summary = {
            "nodes": n_nodes, "edges": n_edges,
            "density": round(density, 4),
            "clustering": round(clustering, 4),
            "avg_degree": round(avg_degree, 2),
            "hub_genes": hub_genes,
        }
        result.detail = {"graph": G, "degrees": degrees}
        return result

    # ==================================================================
    # Layer 2: 细胞/组织尺度
    # ==================================================================

    def analyze_cellular(self, data_type: str = "gene",
                         max_features: int = 100,
                         stratify_by: str = "none",
                         seed_genes: Optional[List[str]] = None,
                         progress_callback=None) -> ScaleResult:
        """
        细胞/组织尺度分析：MRNetB 网络推断

        新增 seed_genes 参数：从分子层传递的Hub基因，用于定向推断
        """
        result = ScaleResult(
            scale="cellular", scale_cn="细胞/组织尺度",
            model_id="mrnetb", model_name="MRNetB网络推断模型",
        )

        if self.tcga_sim is None:
            result.summary = {"error": "TCGA仿真器未初始化"}
            return result

        if self.tcga_sim.gene_data is None:
            try:
                self.tcga_sim.load_data()
            except Exception as e:
                result.summary = {"error": f"数据加载失败: {e}"}
                return result

        try:
            # 如果有种子基因，优先使用定向推断
            if seed_genes and len(seed_genes) >= 2 and data_type == "gene":
                network_df = self.tcga_sim.build_network_for_genes(
                    seed_genes, data_type=data_type,
                    max_features=max_features,
                    progress_callback=progress_callback,
                )
                detail = {"seed_directed": network_df}
                seed_info = f"使用 {len(seed_genes)} 个Hub基因定向推断"
            elif stratify_by == "none":
                data = self.tcga_sim.gene_data if data_type == "gene" else self.tcga_sim.mirna_data
                if data is None:
                    result.summary = {"error": f"无{data_type}数据"}
                    return result
                network_df = self.tcga_sim.build_network_mrnetb(
                    data, data_type=data_type,
                    max_features=max_features,
                    progress_callback=progress_callback,
                )
                detail = {"all": network_df}
                seed_info = "全量推断"
            elif stratify_by == "age":
                detail = self.tcga_sim.analyze_by_age(
                    data_type=data_type, max_features=max_features,
                    progress_callback=progress_callback,
                )
                seed_info = "按年龄分层"
            elif stratify_by == "gender":
                detail = self.tcga_sim.analyze_by_gender(
                    data_type=data_type, max_features=max_features,
                    progress_callback=progress_callback,
                )
                seed_info = "按性别分层"
            elif stratify_by == "stage":
                detail = self.tcga_sim.analyze_by_stage(
                    data_type=data_type, max_features=max_features,
                    progress_callback=progress_callback,
                )
                seed_info = "按阶段分层"
            else:
                detail = {}
                seed_info = ""

            # 构建细胞层NetworkX图用于拓扑分析
            cell_G = nx.Graph()
            total_edges = 0
            for k, df in detail.items():
                if isinstance(df, pd.DataFrame) and len(df) > 0:
                    total_edges += len(df)
                    col1 = "gene1" if "gene1" in df.columns else "mirna1"
                    col2 = "gene2" if "gene2" in df.columns else "mirna2"
                    for _, row in df.iterrows():
                        cell_G.add_edge(row[col1], row[col2], weight=row.get("weight", 1))

            cell_density = nx.density(cell_G) if cell_G.number_of_nodes() > 1 else 0
            cell_clustering = nx.average_clustering(cell_G) if cell_G.number_of_nodes() > 2 else 0

            groups = list(detail.keys())
            result.summary = {
                "groups": groups, "total_edges": total_edges,
                "cell_nodes": cell_G.number_of_nodes(),
                "cell_density": round(cell_density, 4),
                "cell_clustering": round(cell_clustering, 4),
                "seed_info": seed_info,
                "edges_per_group": {
                    k: len(v) for k, v in detail.items() if isinstance(v, pd.DataFrame)
                },
            }
            result.detail = {"networks": detail, "graph": cell_G}

        except Exception as e:
            result.summary = {"error": str(e)}

        return result

    # ==================================================================
    # Layer 3: 群体尺度
    # ==================================================================

    def analyze_population(self, N: int = 100, c: int = 5,
                           beta: float = 0.05, gamma: float = 0.2,
                           ini: float = 0.01, max_step: int = 100) -> ScaleResult:
        """群体尺度分析：SIS 传播动力学"""
        result = ScaleResult(
            scale="population", scale_cn="群体尺度",
            model_id="sis_epidemic", model_name="SIS传播动力学模型",
        )

        if self.social_sim is None:
            result.summary = {"error": "社交网络仿真器未初始化"}
            return result

        try:
            G, communities = self.social_sim.build_community_network(N=N, c=c)
            infected_density, node_state = self.social_sim.SIS_simulation(
                G, beta=beta, gamma=gamma, ini=ini, max_step=max_step,
            )
            metrics = self.social_sim.get_network_metrics(G, communities)

            result.summary = {
                "nodes": G.number_of_nodes(), "edges": G.number_of_edges(),
                "communities": len(communities),
                "beta": round(beta, 4), "gamma": round(gamma, 4),
                "final_density": round(float(infected_density[-1]), 4),
                "peak_density": round(float(np.max(infected_density)), 4),
                "mean_density": round(float(np.mean(infected_density)), 4),
                **metrics,
            }
            result.detail = {
                "graph": G, "communities": communities,
                "infected_density": infected_density, "node_state": node_state,
            }
        except Exception as e:
            result.summary = {"error": str(e)}

        return result

    # ==================================================================
    # 功能1: 跨尺度参数传递 — 完整级联
    # ==================================================================

    @staticmethod
    def _derive_population_params(mol_summary: dict, cell_summary: dict) -> dict:
        """
        从分子层 + 细胞层指标 推导 群体层SIS参数

        方法:
            采用 Hill函数 剂量-响应模型将分子/细胞层拓扑"信号"映射为群体层"行为率"。

            Hill函数: B(s) = B_min + (B_max - B_min) × s^n / (s^n + K^n)

            信号 s 的定义:
                s = (mol_density + cell_clustering) / 2
                即分子层网络密度与细胞层聚类系数的均值，反映跨尺度网络连通性。

            参数:
                β: 感染率 — Hill(s, B_min=0.01, B_max=0.20, K=0.3, n=2)
                c: 社区数 — 由Hub基因数量映射，类比疾病模块数
                γ: 恢复率 — 0.2 (固定基线)
        """
        mol_density = mol_summary.get("density", 0)
        cell_clustering = cell_summary.get("cell_clustering", 0)
        hub_genes = mol_summary.get("hub_genes", [])

        # 跨尺度信号: 分子密度 + 细胞聚类的均值
        s = (mol_density + cell_clustering) / 2.0

        # Hill函数参数
        B_min, B_max = 0.01, 0.20  # β的范围
        K = 0.3   # 半最大信号浓度
        n = 2     # Hill系数（超敏感性）

        # β = B_min + (B_max - B_min) × s^n / (s^n + K^n)
        hill_val = (s ** n) / (s ** n + K ** n) if (s ** n + K ** n) > 0 else 0
        beta = B_min + (B_max - B_min) * hill_val
        beta = round(min(max(beta, 0.01), 0.20), 4)

        # 社区数: 类比疾病模块划分
        c = max(2, min(8, len(hub_genes) // 2))
        N = 100
        gamma = 0.2

        return {
            "beta": beta, "gamma": gamma, "N": N, "c": c,
            "signal_s": round(s, 4),
            "formula_beta": (
                f"Hill函数: β = B_min + (B_max−B_min) × s²/(s²+K²)\n"
                f"  s = (mol_density + cell_clustering)/2 = ({mol_density:.4f}+{cell_clustering:.4f})/2 = {s:.4f}\n"
                f"  β = 0.01 + 0.19 × {s:.4f}²/({s:.4f}²+0.3²) = {beta:.4f}"
            ),
            "formula_c": f"c = max(2, min(8, {len(hub_genes)}//2)) = {c}  [类比疾病模块数]",
            "references": [],
        }

    def run_gene_cascade(self, disease_name: str,
                         network_type: str = "interaction",
                         data_type: str = "gene",
                         max_features: int = 100,
                         progress_callback=None) -> CascadeReport:
        """
        基因网络两层级联：分子层 → 细胞层（不含群体层）

        流程:
            1. 分子层 → Hub基因 + 网络拓扑
            2. 细胞层 → 用Hub基因做定向MRNetB推断
        """
        report = CascadeReport(disease=disease_name)

        # Layer 1: 分子
        if progress_callback:
            progress_callback(0.0, "🧬 分子尺度分析...")
        mol = self.analyze_molecular(disease_name, network_type)
        report.results["molecular"] = mol

        # 参数传递: Hub基因 → 种子节点
        hub_genes = mol.summary.get("hub_genes", []) if "error" not in mol.summary else []
        all_disease_genes = []
        if self.db and disease_name in self.db.diseases:
            all_disease_genes = self.db.diseases[disease_name].genes
        seed_genes = hub_genes + [g for g in all_disease_genes if g not in hub_genes]
        seed_genes = seed_genes[:max_features]

        # Layer 2: 细胞
        if progress_callback:
            progress_callback(0.4, f"🔬 细胞尺度 (种子: {len(seed_genes)} 基因)...")
        cell = self.analyze_cellular(
            data_type=data_type,
            max_features=max_features,
            seed_genes=seed_genes if seed_genes else None,
        )
        report.results["cellular"] = cell

        # 生成洞察
        insights = []
        if "error" not in mol.summary:
            s = mol.summary
            insights.append(
                f"🧬 分子层: {s.get('nodes',0)} 个疾病基因，{s.get('edges',0)} 条互作边，"
                f"密度 {s.get('density',0):.4f}，Hub基因: {', '.join(s.get('hub_genes',[])[:5]) or '无'}。"
            )
        if "error" not in cell.summary:
            s = cell.summary
            insights.append(
                f"🔬 细胞层: 从TCGA-COAD推断出 {s.get('cell_nodes',0)} 个节点、"
                f"{s.get('total_edges',0)} 条功能关联（{s.get('seed_info','')})。"
            )
        if hub_genes and "error" not in cell.summary:
            insights.append(
                f"📐 跨尺度传递: {len(seed_genes)} 个Hub基因作为种子节点定向推断细胞层网络"
                f"（网络邻近性原理）。"
            )
        report.cross_scale_insights = insights or ["未产生洞察。"]

        # 架构HTML
        report.architecture_html = self._gene_cascade_html(report, len(seed_genes))

        if progress_callback:
            progress_callback(1.0, "✅ 完成")
        return report

    def _gene_cascade_html(self, report: CascadeReport, seed_count: int) -> str:
        """基因网络两层级联架构图HTML"""
        mol = report.results.get("molecular")
        cell = report.results.get("cellular")

        def _m(label, val):
            return f"<span style='margin-right:12px;'><b>{label}:</b> {val}</span>"

        mol_metrics = ""
        if mol and "error" not in mol.summary:
            s = mol.summary
            hubs = ", ".join(s.get("hub_genes", [])[:5])
            mol_metrics = (
                f"{_m('节点', s.get('nodes','-'))} {_m('边', s.get('edges','-'))} "
                f"{_m('密度', s.get('density','-'))} {_m('聚类', s.get('clustering','-'))}"
                f"<br/>{_m('Hub基因', hubs or '-')}"
            )
        cell_metrics = ""
        if cell and "error" not in cell.summary:
            s = cell.summary
            cell_metrics = (
                f"{_m('节点', s.get('cell_nodes','-'))} {_m('边', s.get('total_edges','-'))} "
                f"{_m('密度', s.get('cell_density','-'))} {_m('聚类', s.get('cell_clustering','-'))}"
            )

        return f"""
        <div style="font-family:system-ui;padding:16px;display:flex;flex-direction:column;align-items:center;gap:0;">
          <h3>基因网络级联分析 — {report.disease}</h3>
          <div style="background:linear-gradient(135deg,#667eea,#764ba2);color:#fff;
                      border-radius:12px;padding:14px 20px;width:100%;max-width:600px;">
            <b>🧬 Layer 1 — 分子尺度</b><br/><span style="font-size:0.85em;">基因互作/调控网络</span>
            <div style="font-size:0.82em;margin-top:6px;">{mol_metrics}</div>
          </div>
          <div style="text-align:center;margin:3px 0;">⬇️ <span style="font-size:0.75em;color:#555;">传递 {seed_count} 个Hub基因 → 种子节点</span></div>
          <div style="background:linear-gradient(135deg,#f093fb,#f5576c);color:#fff;
                      border-radius:12px;padding:14px 20px;width:100%;max-width:600px;">
            <b>🔬 Layer 2 — 细胞/组织尺度</b><br/><span style="font-size:0.85em;">MRNetB网络推断（TCGA-COAD）</span>
            <div style="font-size:0.82em;margin-top:6px;">{cell_metrics}</div>
          </div>
        </div>
        """

    def run_full_cascade(self, disease_name: str,
                         network_type: str = "interaction",
                         data_type: str = "gene",
                         max_features: int = 100,
                         progress_callback=None) -> CascadeReport:
        """
        完整三层级联分析 — 带跨尺度参数传递

        流程:
            1. 分子层 → 提取 hub_genes + density + clustering
            2. 细胞层 → 用 hub_genes 做定向MRNetB推断 → 提取 cell_clustering
            3. 群体层 → 用公式推导 β, c → 运行SIS仿真
        """
        report = CascadeReport(disease=disease_name)

        # --- Layer 1: 分子 ---
        if progress_callback:
            progress_callback(0.0, "🧬 Layer 1 / 分子尺度分析...")
        mol = self.analyze_molecular(disease_name, network_type)
        report.results["molecular"] = mol

        # --- 参数传递: 分子 → 细胞 ---
        hub_genes = mol.summary.get("hub_genes", []) if "error" not in mol.summary else []
        all_disease_genes = []
        if self.db and disease_name in self.db.diseases:
            all_disease_genes = self.db.diseases[disease_name].genes
        # 种子基因 = Hub基因(优先) + 其余疾病基因, 取前 max_features 个
        seed_genes = hub_genes + [g for g in all_disease_genes if g not in hub_genes]
        seed_genes = seed_genes[:max_features]

        # --- Layer 2: 细胞 ---
        if progress_callback:
            progress_callback(0.3, f"🔬 Layer 2 / 细胞尺度 (种子基因: {len(seed_genes)})...")
        cell = self.analyze_cellular(
            data_type=data_type,
            max_features=max_features,
            seed_genes=seed_genes if seed_genes else None,
        )
        report.results["cellular"] = cell

        # --- 参数传递: 分子+细胞 → 群体 ---
        pop_params = self._derive_population_params(mol.summary, cell.summary)
        report.param_transfer = {
            "seed_genes": seed_genes[:10],
            "seed_count": len(seed_genes),
            **pop_params,
        }

        # --- Layer 3: 群体 ---
        if progress_callback:
            progress_callback(0.7, f"👥 Layer 3 / 群体尺度 (β={pop_params['beta']}, c={pop_params['c']})...")
        pop = self.analyze_population(
            N=pop_params["N"], c=pop_params["c"],
            beta=pop_params["beta"], gamma=pop_params["gamma"],
        )
        report.results["population"] = pop

        # --- 洞察 ---
        report.cross_scale_insights = self._generate_full_insights(report)
        report.architecture_html = self.create_cascade_html(report)

        if progress_callback:
            progress_callback(1.0, "✅ 三层级联分析完成")

        return report

    def _generate_full_insights(self, report: CascadeReport) -> List[str]:
        """三层级联洞察（含文献依据）"""
        insights = []
        mol = report.results.get("molecular", ScaleResult("", "", "", ""))
        cell = report.results.get("cellular", ScaleResult("", "", "", ""))
        pop = report.results.get("population", ScaleResult("", "", "", ""))
        pt = report.param_transfer

        if "error" not in mol.summary:
            n = mol.summary.get("nodes", 0)
            d = mol.summary.get("density", 0)
            hubs = mol.summary.get("hub_genes", [])
            insights.append(
                f"🧬 分子层: {n} 个疾病基因，网络密度 {d:.4f}，"
                f"Hub基因: {', '.join(hubs[:5]) if hubs else '无'}。"
            )

        if "error" not in cell.summary:
            te = cell.summary.get("total_edges", 0)
            si = cell.summary.get("seed_info", "")
            cn = cell.summary.get("cell_nodes", 0)
            insights.append(
                f"🔬 细胞层: {si}，从 TCGA-COAD 推断出 {cn} 个节点、{te} 条功能关联。"
                f" Hub基因作为种子节点定向推断，"
                f"依据网络邻近性原理。"
            )

        if pt:
            insights.append(
                f"📐 参数传递 (Hill函数, PhysiCell框架): "
                f"跨尺度信号 s={pt.get('signal_s', '?')}，"
                f"{pt.get('formula_beta', '').split(chr(10))[-1].strip()}。"
            )

        if "error" not in pop.summary:
            peak = pop.summary.get("peak_density", 0)
            final = pop.summary.get("final_density", 0)
            beta = pop.summary.get("beta", 0)
            insights.append(
                f"👥 群体层: β={beta:.4f} → 峰值感染 {peak:.2%}，最终感染 {final:.2%}。"
            )
            # 基本再生数 R0 近似
            gamma = pt.get("gamma", 0.2)
            avg_deg = pop.summary.get("avg_degree", 0)
            if avg_deg > 0 and gamma > 0:
                R0_approx = beta * avg_deg / gamma
                insights.append(
                    f"📊 基本再生数近似: R₀ ≈ β⟨k⟩/γ = {beta:.4f}×{avg_deg:.1f}/{gamma:.2f} = {R0_approx:.2f}。"
                    f" {'R₀>1 → 疾病可在网络中持续传播。' if R0_approx > 1 else 'R₀<1 → 传播趋于消亡。'}"
                )

        refs = pt.get("references", [])
        if refs:
            insights.append(
                "📚 方法论参考: " + "; ".join(refs) + "。"
            )

        if not insights:
            insights.append("未产生跨尺度洞察；请确认各层分析均成功运行。")

        return insights

    def create_cascade_html(self, report: CascadeReport) -> str:
        """三层级联架构图 + 参数传递公式"""
        mol = report.results.get("molecular")
        cell = report.results.get("cellular")
        pop = report.results.get("population")
        pt = report.param_transfer

        def _m(label, val):
            return f"<span style='margin-right:12px;'><b>{label}:</b> {val}</span>"

        def _layer(color1, color2, icon, title, subtitle, metrics_html, text_color="#fff"):
            return f"""
            <div style="background:linear-gradient(135deg,{color1},{color2});color:{text_color};
                        border-radius:12px;padding:14px 20px;width:100%;max-width:600px;">
              <b>{icon} {title}</b><br/>
              <span style="font-size:0.85em;">{subtitle}</span>
              <div style="font-size:0.82em;margin-top:6px;">{metrics_html}</div>
            </div>"""

        def _arrow(formula=""):
            label = f'<span style="font-size:0.75em;color:#555;">{formula}</span>' if formula else ""
            return f'<div style="text-align:center;margin:2px 0;">⬇️ {label}</div>'

        # Layer 1
        mol_metrics = ""
        if mol and "error" not in mol.summary:
            s = mol.summary
            hubs = ", ".join(s.get("hub_genes", [])[:5])
            mol_metrics = (
                f"{_m('节点', s.get('nodes','-'))} {_m('边', s.get('edges','-'))} "
                f"{_m('密度', s.get('density','-'))} {_m('聚类', s.get('clustering','-'))}"
                f"<br/>{_m('Hub基因', hubs or '-')}"
            )
        l1 = _layer("#667eea", "#764ba2", "🧬", "Layer 1 — 分子尺度",
                     "基因互作/调控网络拓扑分析", mol_metrics)

        # Arrow 1→2
        a1 = _arrow(f"传递 Hub基因 → 种子节点 ({pt.get('seed_count', '?')} 个)")

        # Layer 2
        cell_metrics = ""
        if cell and "error" not in cell.summary:
            s = cell.summary
            cell_metrics = (
                f"{_m('节点', s.get('cell_nodes','-'))} {_m('边', s.get('total_edges','-'))} "
                f"{_m('密度', s.get('cell_density','-'))} {_m('聚类', s.get('cell_clustering','-'))}"
                f"<br/>{_m('模式', s.get('seed_info','-'))}"
            )
        l2 = _layer("#f093fb", "#f5576c", "🔬", "Layer 2 — 细胞/组织尺度",
                     "MRNetB 网络推断（TCGA-COAD 表达数据）", cell_metrics)

        # Arrow 2→3
        a2 = _arrow(pt.get("formula_beta", ""))

        # Layer 3
        pop_metrics = ""
        if pop and "error" not in pop.summary:
            s = pop.summary
            pop_metrics = (
                f"{_m('β', s.get('beta','-'))} {_m('γ', s.get('gamma','-'))} "
                f"{_m('社区', s.get('communities','-'))}"
                f"<br/>{_m('峰值感染', s.get('peak_density','-'))} "
                f"{_m('最终感染', s.get('final_density','-'))}"
            )
        l3 = _layer("#43e97b", "#38f9d7", "👥", "Layer 3 — 群体尺度",
                     "SIS 传播动力学仿真", pop_metrics, text_color="#333")

        html = f"""
        <div style="font-family:system-ui,sans-serif;padding:16px;
                    display:flex;flex-direction:column;align-items:center;gap:0;">
          <h3 style="margin-bottom:10px;">跨尺度级联分析 — {report.disease}</h3>
          {l1}{a1}{l2}{a2}{l3}
        </div>
        """
        return html

    # ==================================================================
    # 功能2: 跨尺度对比可视化
    # ==================================================================

    def create_radar_chart(self, report: CascadeReport) -> go.Figure:
        """
        跨尺度雷达图：对比三层的标准化拓扑指标
        """
        categories = ["网络密度", "聚类系数", "平均度(归一化)", "节点数(归一化)", "边数(归一化)"]
        mol = report.results.get("molecular")
        cell = report.results.get("cellular")
        pop = report.results.get("population")

        def _norm(val, max_val):
            return min(val / max_val, 1.0) if max_val > 0 else 0

        # 收集原始值
        mol_vals = [0, 0, 0, 0, 0]
        cell_vals = [0, 0, 0, 0, 0]
        pop_vals = [0, 0, 0, 0, 0]

        max_nodes = 1
        max_edges = 1
        max_degree = 1

        if mol and "error" not in mol.summary:
            s = mol.summary
            mol_vals = [s.get("density", 0), s.get("clustering", 0),
                        s.get("avg_degree", 0), s.get("nodes", 0), s.get("edges", 0)]
            max_nodes = max(max_nodes, s.get("nodes", 0))
            max_edges = max(max_edges, s.get("edges", 0))
            max_degree = max(max_degree, s.get("avg_degree", 0))

        if cell and "error" not in cell.summary:
            s = cell.summary
            cell_vals = [s.get("cell_density", 0), s.get("cell_clustering", 0),
                         0, s.get("cell_nodes", 0), s.get("total_edges", 0)]
            max_nodes = max(max_nodes, s.get("cell_nodes", 0))
            max_edges = max(max_edges, s.get("total_edges", 0))

        if pop and "error" not in pop.summary:
            s = pop.summary
            pop_vals = [s.get("density", 0), s.get("clustering", 0),
                        s.get("avg_degree", 0), s.get("nodes", 0), s.get("edges", 0)]
            max_nodes = max(max_nodes, s.get("nodes", 0))
            max_edges = max(max_edges, s.get("edges", 0))
            max_degree = max(max_degree, s.get("avg_degree", 0))

        # 归一化
        def normalize(vals):
            return [
                vals[0],  # density 已在 0-1
                vals[1],  # clustering 已在 0-1
                _norm(vals[2], max_degree),
                _norm(vals[3], max_nodes),
                _norm(vals[4], max_edges),
            ]

        mol_n = normalize(mol_vals)
        cell_n = normalize(cell_vals)
        pop_n = normalize(pop_vals)

        fig = go.Figure()
        fig.add_trace(go.Scatterpolar(r=mol_n + [mol_n[0]], theta=categories + [categories[0]],
                                       fill='toself', name='🧬 分子层',
                                       line=dict(color='#764ba2')))
        fig.add_trace(go.Scatterpolar(r=cell_n + [cell_n[0]], theta=categories + [categories[0]],
                                       fill='toself', name='🔬 细胞层',
                                       line=dict(color='#f5576c')))
        fig.add_trace(go.Scatterpolar(r=pop_n + [pop_n[0]], theta=categories + [categories[0]],
                                       fill='toself', name='👥 群体层',
                                       line=dict(color='#43e97b')))
        fig.update_layout(
            polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
            title=dict(text="跨尺度拓扑指标对比", x=0.5),
            height=500, showlegend=True,
        )
        return fig

    def create_two_network_plots(self, report: CascadeReport) -> go.Figure:
        """两层并排网络缩略图 (分子层 + 细胞层)"""
        fig = make_subplots(rows=1, cols=2,
                            subplot_titles=["🧬 分子层网络", "🔬 细胞层网络"],
                            horizontal_spacing=0.08)

        def _add_network(G, col, color, max_nodes=50):
            if G is None or G.number_of_nodes() == 0:
                fig.add_trace(go.Scatter(x=[0], y=[0], mode='text',
                                          text=["无数据"], textfont=dict(size=14)),
                              row=1, col=col)
                return
            if G.number_of_nodes() > max_nodes:
                top_nodes = sorted(G.degree(), key=lambda x: x[1], reverse=True)[:max_nodes]
                G = G.subgraph([n[0] for n in top_nodes]).copy()
            pos = nx.spring_layout(G, seed=42, k=0.8, iterations=30)
            ex, ey = [], []
            for e in G.edges():
                x0, y0 = pos[e[0]]; x1, y1 = pos[e[1]]
                ex.extend([x0, x1, None]); ey.extend([y0, y1, None])
            fig.add_trace(go.Scatter(x=ex, y=ey, mode='lines',
                                      line=dict(width=0.5, color='#ccc'),
                                      hoverinfo='none', showlegend=False), row=1, col=col)
            nl = list(G.nodes())
            nx_arr = [pos[n][0] for n in nl]; ny_arr = [pos[n][1] for n in nl]
            degs = [G.degree(n) for n in nl]; md = max(degs) if degs else 1
            sizes = [5 + 15 * d / md for d in degs]
            fig.add_trace(go.Scatter(x=nx_arr, y=ny_arr, mode='markers',
                                      marker=dict(size=sizes, color=color, opacity=0.8,
                                                  line=dict(width=0.5, color='white')),
                                      text=[str(n) for n in nl], hoverinfo='text',
                                      showlegend=False), row=1, col=col)

        mol = report.results.get("molecular")
        mol_G = mol.detail.get("graph") if mol and isinstance(mol.detail, dict) else None
        _add_network(mol_G, 1, '#764ba2')
        cell = report.results.get("cellular")
        cell_G = cell.detail.get("graph") if cell and isinstance(cell.detail, dict) else None
        _add_network(cell_G, 2, '#f5576c')
        for i in range(1, 3):
            fig.update_xaxes(showgrid=False, showticklabels=False, zeroline=False, row=1, col=i)
            fig.update_yaxes(showgrid=False, showticklabels=False, zeroline=False, row=1, col=i)
        fig.update_layout(height=450, title=dict(text="基因网络两层结构对比", x=0.5),
                          plot_bgcolor='white')
        return fig

    def create_gene_radar_chart(self, report: CascadeReport) -> go.Figure:
        """基因网络两层雷达图（分子层 + 细胞层）"""
        categories = ["网络密度", "聚类系数", "节点数(归一化)", "边数(归一化)"]
        mol = report.results.get("molecular")
        cell = report.results.get("cellular")
        def _norm(v, mx): return min(v / mx, 1.0) if mx > 0 else 0
        mol_vals = [0, 0, 0, 0]; cell_vals = [0, 0, 0, 0]
        max_n = 1; max_e = 1
        if mol and "error" not in mol.summary:
            s = mol.summary
            mol_vals = [s.get("density", 0), s.get("clustering", 0), s.get("nodes", 0), s.get("edges", 0)]
            max_n = max(max_n, s.get("nodes", 0)); max_e = max(max_e, s.get("edges", 0))
        if cell and "error" not in cell.summary:
            s = cell.summary
            cell_vals = [s.get("cell_density", 0), s.get("cell_clustering", 0), s.get("cell_nodes", 0), s.get("total_edges", 0)]
            max_n = max(max_n, s.get("cell_nodes", 0)); max_e = max(max_e, s.get("total_edges", 0))
        mol_n = [mol_vals[0], mol_vals[1], _norm(mol_vals[2], max_n), _norm(mol_vals[3], max_e)]
        cell_n = [cell_vals[0], cell_vals[1], _norm(cell_vals[2], max_n), _norm(cell_vals[3], max_e)]
        fig = go.Figure()
        fig.add_trace(go.Scatterpolar(r=mol_n + [mol_n[0]], theta=categories + [categories[0]],
                                       fill='toself', name='🧬 分子层', line=dict(color='#764ba2')))
        fig.add_trace(go.Scatterpolar(r=cell_n + [cell_n[0]], theta=categories + [categories[0]],
                                       fill='toself', name='🔬 细胞层', line=dict(color='#f5576c')))
        fig.update_layout(polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
                          title=dict(text="基因网络两层拓扑对比", x=0.5), height=450, showlegend=True)
        return fig

    def create_three_network_plots(self, report: CascadeReport) -> go.Figure:
        """三层并排网络缩略图 (1×3 subplots)"""
        fig = make_subplots(rows=1, cols=3,
                            subplot_titles=["🧬 分子层网络", "🔬 细胞层网络", "👥 群体层网络"],
                            horizontal_spacing=0.05)

        def _add_network(G, col, color, max_nodes=50):
            if G is None or G.number_of_nodes() == 0:
                fig.add_trace(go.Scatter(x=[0], y=[0], mode='text',
                                          text=["无数据"], textfont=dict(size=14)),
                              row=1, col=col)
                return

            # 取前max_nodes个节点（按度数排序）
            if G.number_of_nodes() > max_nodes:
                top_nodes = sorted(G.degree(), key=lambda x: x[1], reverse=True)[:max_nodes]
                top_node_ids = [n[0] for n in top_nodes]
                G = G.subgraph(top_node_ids).copy()

            pos = nx.spring_layout(G, seed=42, k=0.8, iterations=30)

            # 边
            ex, ey = [], []
            for e in G.edges():
                x0, y0 = pos[e[0]]
                x1, y1 = pos[e[1]]
                ex.extend([x0, x1, None])
                ey.extend([y0, y1, None])
            fig.add_trace(go.Scatter(x=ex, y=ey, mode='lines',
                                      line=dict(width=0.5, color='#ccc'),
                                      hoverinfo='none', showlegend=False),
                          row=1, col=col)

            # 节点
            nx_list = list(G.nodes())
            node_x = [pos[n][0] for n in nx_list]
            node_y = [pos[n][1] for n in nx_list]
            degrees = [G.degree(n) for n in nx_list]
            max_deg = max(degrees) if degrees else 1
            sizes = [5 + 15 * d / max_deg for d in degrees]
            fig.add_trace(go.Scatter(x=node_x, y=node_y, mode='markers',
                                      marker=dict(size=sizes, color=color, opacity=0.8,
                                                  line=dict(width=0.5, color='white')),
                                      text=[str(n) for n in nx_list],
                                      hoverinfo='text', showlegend=False),
                          row=1, col=col)

        # 分子层
        mol = report.results.get("molecular")
        mol_G = mol.detail.get("graph") if mol and isinstance(mol.detail, dict) else None
        _add_network(mol_G, 1, '#764ba2')

        # 细胞层
        cell = report.results.get("cellular")
        cell_G = cell.detail.get("graph") if cell and isinstance(cell.detail, dict) else None
        _add_network(cell_G, 2, '#f5576c')

        # 群体层
        pop = report.results.get("population")
        pop_G = pop.detail.get("graph") if pop and isinstance(pop.detail, dict) else None
        _add_network(pop_G, 3, '#43e97b')

        for i in range(1, 4):
            fig.update_xaxes(showgrid=False, showticklabels=False, zeroline=False, row=1, col=i)
            fig.update_yaxes(showgrid=False, showticklabels=False, zeroline=False, row=1, col=i)

        fig.update_layout(height=450, title=dict(text="三层网络结构对比", x=0.5),
                          plot_bgcolor='white')
        return fig

    # ==================================================================
    # 功能3: 多疾病对比分析
    # ==================================================================

    def compare_diseases(self, disease_names: List[str],
                         network_type: str = "interaction") -> Tuple[pd.DataFrame, go.Figure]:
        """
        多疾病分子层对比

        返回:
            (对比表DataFrame, 分组柱状图Figure)
        """
        rows = []
        for dn in disease_names:
            result = self.analyze_molecular(dn, network_type)
            s = result.summary
            if "error" in s:
                rows.append({"疾病": dn, "状态": f"❌ {s['error']}"})
            else:
                rows.append({
                    "疾病": dn,
                    "基因数": s.get("nodes", 0),
                    "边数": s.get("edges", 0),
                    "网络密度": s.get("density", 0),
                    "聚类系数": s.get("clustering", 0),
                    "平均度": s.get("avg_degree", 0),
                    "Top Hub": ", ".join(s.get("hub_genes", [])[:3]),
                })
        df = pd.DataFrame(rows)

        # 柱状图
        metrics = ["基因数", "边数", "网络密度", "聚类系数", "平均度"]
        existing_metrics = [m for m in metrics if m in df.columns]

        fig = make_subplots(rows=1, cols=len(existing_metrics),
                            subplot_titles=existing_metrics,
                            horizontal_spacing=0.08)

        colors = ['#667eea', '#f093fb', '#43e97b', '#f5576c', '#ffa07a',
                  '#4ecdc4', '#bb8fce', '#f7dc6f']
        for col_idx, metric in enumerate(existing_metrics, 1):
            for i, row in df.iterrows():
                fig.add_trace(go.Bar(
                    x=[row["疾病"]], y=[row.get(metric, 0)],
                    name=row["疾病"] if col_idx == 1 else None,
                    marker_color=colors[i % len(colors)],
                    showlegend=(col_idx == 1),
                    legendgroup=row["疾病"],
                ), row=1, col=col_idx)

        fig.update_layout(height=400, title=dict(text="多疾病分子网络对比", x=0.5),
                          barmode='group')
        return df, fig

    # ==================================================================
    # 功能4: 基因追踪
    # ==================================================================

    def trace_gene(self, gene_name: str, disease_name: str,
                   network_type: str = "interaction") -> str:
        """
        追踪单个基因在三层尺度中的角色

        返回:
            基因档案卡 HTML
        """
        card_sections = []

        # === 分子层 ===
        mol_info = {"found": False}
        if self.db:
            disease = self.db.diseases.get(disease_name)
            if disease:
                G = nx.Graph() if network_type == "interaction" else nx.DiGraph()
                for g in disease.genes:
                    G.add_node(g)
                relations = disease.interactions if network_type == "interaction" else disease.regulations
                for rel in relations:
                    if network_type == "interaction":
                        G.add_edge(rel.gene1, rel.gene2)
                    else:
                        G.add_edge(rel.regulator, rel.target)

                if gene_name in G:
                    degree = G.degree(gene_name)
                    all_degrees = dict(G.degree())
                    rank = sorted(all_degrees.values(), reverse=True).index(degree) + 1
                    neighbors = list(G.neighbors(gene_name))[:10]

                    # 查找所在通路
                    gene_obj = self.db.genes.get(gene_name)
                    pathways = gene_obj.pathways if gene_obj else []

                    mol_info = {
                        "found": True, "degree": degree, "rank": rank,
                        "total": G.number_of_nodes(),
                        "is_hub": rank <= 10,
                        "neighbors": neighbors,
                        "pathways": pathways[:5],
                    }

        mol_html = self._gene_card_section(
            "🧬 分子尺度", "#667eea",
            f"""
            <b>度数:</b> {mol_info.get('degree', '-')} (排名 #{mol_info.get('rank', '-')}/{mol_info.get('total', '-')})<br/>
            <b>Hub基因:</b> {'✅ 是' if mol_info.get('is_hub') else '❌ 否'}<br/>
            <b>邻居基因:</b> {', '.join(mol_info.get('neighbors', []))}<br/>
            <b>所在通路:</b> {', '.join(mol_info.get('pathways', [])) or '无'}
            """ if mol_info["found"] else f"<i>{gene_name} 不在 {disease_name} 的基因网络中</i>"
        )
        card_sections.append(mol_html)

        # === 细胞层 ===
        cell_html_content = f"<i>TCGA数据未加载</i>"
        if self.tcga_sim:
            if self.tcga_sim.gene_data is None:
                try:
                    self.tcga_sim.load_data()
                except Exception:
                    pass

            if self.tcga_sim.gene_data is not None and gene_name in self.tcga_sim.gene_data.index:
                expr = self.tcga_sim.gene_data.loc[gene_name]
                mean_expr = float(expr.mean())
                std_expr = float(expr.std())
                max_expr = float(expr.max())

                # 表达排名
                all_means = self.tcga_sim.gene_data.mean(axis=1).sort_values(ascending=False)
                expr_rank = list(all_means.index).index(gene_name) + 1 if gene_name in all_means.index else "-"
                total_genes = len(all_means)

                cell_html_content = (
                    f"<b>平均表达量:</b> {mean_expr:.2f} (排名 #{expr_rank}/{total_genes})<br/>"
                    f"<b>标准差:</b> {std_expr:.2f}<br/>"
                    f"<b>最大值:</b> {max_expr:.2f}<br/>"
                    f"<b>表达变异度:</b> {'高' if std_expr / (mean_expr + 1e-6) > 1 else '中' if std_expr / (mean_expr + 1e-6) > 0.5 else '低'}"
                )
            elif self.tcga_sim.gene_data is not None:
                cell_html_content = f"<i>{gene_name} 不在 TCGA-COAD 基因表达数据中</i>"

        cell_html = self._gene_card_section("🔬 细胞/组织尺度", "#f5576c", cell_html_content)
        card_sections.append(cell_html)

        # 组装卡片
        return f"""
        <div style="font-family:system-ui,sans-serif;padding:16px;max-width:700px;margin:auto;">
          <h3 style="text-align:center;">🔍 基因追踪: <code>{gene_name}</code> @ {disease_name}</h3>
          <div style="display:flex;flex-direction:column;gap:8px;margin-top:12px;">
            {''.join(card_sections)}
          </div>
        </div>
        """

    @staticmethod
    def _gene_card_section(title: str, color: str, content: str,
                           text_color: str = "#fff") -> str:
        return f"""
        <div style="background:linear-gradient(135deg,{color},
                    {color}cc);color:{text_color};
                    border-radius:10px;padding:14px 18px;">
          <b>{title}</b>
          <div style="font-size:0.88em;margin-top:6px;line-height:1.6;">
            {content}
          </div>
        </div>
        """

    # ==================================================================
    # 通用工具
    # ==================================================================

    @staticmethod
    def _insights_html(insights: List[str]) -> str:
        items = "".join(f"<li style='margin-bottom:6px;'>{i}</li>" for i in insights)
        return f"""
        <div style="margin-top:16px;background:#f8f9fa;border-left:4px solid #667eea;
                    padding:12px 16px;border-radius:0 8px 8px 0;">
          <b>🔗 跨尺度洞察</b>
          <ul style="margin:8px 0 0 0;padding-left:18px;font-size:0.9em;">{items}</ul>
        </div>
        """

    def cascade_summary_df(self, report: CascadeReport) -> pd.DataFrame:
        """将级联报告转换为汇总 DataFrame"""
        rows = []
        for scale_key, sr in report.results.items():
            row = {
                "尺度": sr.scale_cn, "模型": sr.model_name,
                "状态": "❌ " + sr.summary.get("error", "") if "error" in sr.summary else "✅ 成功",
            }
            for k, v in sr.summary.items():
                if k in ("error", "edges_per_group"):
                    continue
                if isinstance(v, (int, float)):
                    row[k] = v
                elif isinstance(v, list) and len(v) <= 5:
                    row[k] = ", ".join(str(x) for x in v)
                elif isinstance(v, str):
                    row[k] = v
            rows.append(row)
        return pd.DataFrame(rows)
