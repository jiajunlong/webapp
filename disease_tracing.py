"""
疾病多尺度溯源分析模块 - Disease Multi-Scale Tracing

实现从疾病 → 信号通路 → 基因模块 → 基因表达特征的逐层下钻溯源。

方法论依据:
    [1] Weighted Gene Networks Derived from Multi-Omics Reveal Core Cancer Genes
        — 加权网络构建 + 变异系数(CV)加权通路影响力
    [2] Node Attribute Entropy in Complex Network Analysis
        — NAE评分: 综合度中心性、介数中心性、表达方差、通路特异性
    [3] Decoding Colon Cancer Heterogeneity Through Integrated miRNA-Gene Network
        — miRNA-基因负相关调控关系映射
"""

import pandas as pd
import numpy as np
import networkx as nx
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats
from typing import Dict, List, Tuple, Optional
import logging
import os

logger = logging.getLogger(__name__)


class DiseaseTracer:
    """疾病多尺度溯源分析器"""

    def __init__(self):
        self.gene_disease = None
        self.pathway_data = None
        self.expr_data = None
        self.mirna_data = None
        self.clinical_data = None
        self._disease_cache = {}

    def load_data(self):
        """加载所有数据"""
        # gene_disease.tsv
        gd_file = "data/gene_disease.tsv"
        if os.path.exists(gd_file):
            self.gene_disease = pd.read_csv(gd_file, sep='\t')
            logger.info(f"✓ gene_disease: {self.gene_disease.shape}")

        # pathway
        pw_file = "data/pathway(基因名映射版).tsv"
        if os.path.exists(pw_file):
            self.pathway_data = pd.read_csv(pw_file, sep='\t')
            logger.info(f"✓ pathways: {len(self.pathway_data)}")

        # expression (gene symbols)
        for f in ["TCGA-COAD/filtered_hiseq_data.csv", "data/TCGA-COAD/filtered_hiseq_data.csv"]:
            if os.path.exists(f):
                self.expr_data = pd.read_csv(f, index_col=0)
                if not str(self.expr_data.index[0]).startswith('ENSG'):
                    logger.info(f"✓ expression: {self.expr_data.shape}")
                    break
                self.expr_data = None

        # miRNA
        for f in ["TCGA-COAD/filtered_miRNA_with_names.csv", "data/TCGA-COAD/filtered_miRNA_with_names.csv"]:
            if os.path.exists(f):
                self.mirna_data = pd.read_csv(f, index_col=0)
                logger.info(f"✓ miRNA: {self.mirna_data.shape}")
                break

        # clinical
        clin_tsv = "TCGA-COAD/clinical.tsv"
        if os.path.exists(clin_tsv) and self.expr_data is not None:
            raw = pd.read_csv(clin_tsv, sep='\t')
            if 'case_submitter_id' in raw.columns:
                c_idx = raw.set_index('case_submitter_id')
                rows = {}
                for col in self.expr_data.columns:
                    p = col[:12]
                    if p in c_idx.index:
                        r = c_idx.loc[p]
                        rows[col] = r.iloc[0] if isinstance(r, pd.DataFrame) else r
                self.clinical_data = pd.DataFrame(rows).T
                if 'age_at_index' in self.clinical_data.columns:
                    self.clinical_data['age_at_index'] = pd.to_numeric(
                        self.clinical_data['age_at_index'], errors='coerce')
                    self.clinical_data['Age_Group'] = pd.cut(
                        self.clinical_data['age_at_index'],
                        bins=[0, 50, 70, 200], labels=['young', 'middle', 'old'])
                logger.info(f"✓ clinical: {self.clinical_data.shape}")

    # ==========================================================
    # Layer 1: 疾病概览
    # ==========================================================

    def get_disease_list(self) -> List[str]:
        """获取有足够基因的疾病列表"""
        if self.gene_disease is None:
            return []
        counts = self.gene_disease.groupby('disease_name')['gene_symbol'].nunique()
        return sorted([d for d, c in counts.items() if c >= 5], key=lambda x: -counts[x])

    def get_disease_overview(self, disease_name: str) -> dict:
        """Layer 1: 疾病概览信息"""
        if self.gene_disease is None:
            return {"error": "数据未加载"}

        rows = self.gene_disease[self.gene_disease['disease_name'] == disease_name]
        if rows.empty:
            return {"error": f"未找到疾病: {disease_name}"}

        # 提取基因
        genes = list(set(
            str(g).split(',')[0].strip()
            for g in rows['gene_symbol'].dropna()
            if str(g).strip() and str(g).strip() != 'NA'
        ))

        # 提取通路
        all_pathways = set()
        for pw_str in rows['gene_pathway'].dropna():
            for pw in str(pw_str).split(';'):
                pw = pw.strip()
                if pw:
                    all_pathways.add(pw)

        # 疾病分类
        category = rows['disease_category'].dropna().iloc[0] if 'disease_category' in rows.columns and rows['disease_category'].notna().any() else "未知"

        # 表达数据中可用的基因
        avail_genes = [g for g in genes if self.expr_data is not None and g in self.expr_data.index]

        return {
            "disease": disease_name,
            "category": category,
            "n_genes": len(genes),
            "n_genes_in_tcga": len(avail_genes),
            "genes": genes,
            "available_genes": avail_genes,
            "n_pathways": len(all_pathways),
            "pathways": sorted(all_pathways),
        }

    # ==========================================================
    # Layer 2: 通路影响力排名
    # ==========================================================

    def get_pathway_ranking(self, disease_name: str) -> pd.DataFrame:
        """
        Layer 2: 按通路影响力评分排序

        方法 (论文1: Weighted Gene Networks):
            PathwayScore(p, d) = |G_p ∩ G_d| × mean(CV(g))
            CV(g) = std(expr_g) / (mean(expr_g) + ε)
        """
        overview = self.get_disease_overview(disease_name)
        if "error" in overview:
            return pd.DataFrame()

        disease_genes = set(overview["available_genes"])
        if not disease_genes or self.expr_data is None:
            return pd.DataFrame()

        # 构建 基因 → 通路 映射
        rows = self.gene_disease[self.gene_disease['disease_name'] == disease_name]
        gene_to_pathways = {}
        for _, r in rows.iterrows():
            gene = str(r['gene_symbol']).split(',')[0].strip()
            if gene in disease_genes and pd.notna(r.get('gene_pathway')):
                pws = [p.strip() for p in str(r['gene_pathway']).split(';') if p.strip()]
                gene_to_pathways[gene] = pws

        # 反向: 通路 → 基因
        pathway_genes = {}
        for gene, pws in gene_to_pathways.items():
            for pw in pws:
                pathway_genes.setdefault(pw, []).append(gene)

        # 计算每个基因的 CV (表达变异系数)
        gene_cv = {}
        for gene in disease_genes:
            if gene in self.expr_data.index:
                expr = self.expr_data.loc[gene].values.astype(float)
                mean_val = np.mean(expr)
                std_val = np.std(expr)
                gene_cv[gene] = std_val / (abs(mean_val) + 1e-6)

        # 计算通路影响力
        results = []
        for pw, genes in pathway_genes.items():
            n_genes = len(genes)
            cvs = [gene_cv.get(g, 0) for g in genes]
            mean_cv = np.mean(cvs) if cvs else 0
            score = n_genes * mean_cv  # 论文1: 加权评分

            results.append({
                "通路": pw,
                "疾病基因数": n_genes,
                "平均变异系数(CV)": round(mean_cv, 4),
                "影响力评分": round(score, 4),
                "基因列表": ", ".join(sorted(genes)),
            })

        df = pd.DataFrame(results).sort_values("影响力评分", ascending=False)
        return df.reset_index(drop=True)

    # ==========================================================
    # Layer 3: 基因模块 + NAE 排序
    # ==========================================================

    def get_gene_module(self, disease_name: str, pathway_name: str) -> Tuple[pd.DataFrame, nx.Graph]:
        """
        Layer 3: 通路内基因的 NAE 排序 + 相关性网络

        方法 (论文2: Node Attribute Entropy):
            NAE(g) = -Σ p_k × log2(p_k)
            属性: 度中心性, 介数中心性, 表达方差, 通路数量
        """
        overview = self.get_disease_overview(disease_name)
        if "error" in overview:
            return pd.DataFrame(), nx.Graph()

        # 获取通路内的疾病基因
        rows = self.gene_disease[self.gene_disease['disease_name'] == disease_name]
        pathway_genes = []
        for _, r in rows.iterrows():
            gene = str(r['gene_symbol']).split(',')[0].strip()
            if pd.notna(r.get('gene_pathway')):
                pws = [p.strip() for p in str(r['gene_pathway']).split(';')]
                if pathway_name in pws and gene in overview["available_genes"]:
                    pathway_genes.append(gene)

        pathway_genes = list(set(pathway_genes))
        if len(pathway_genes) < 2:
            return pd.DataFrame(), nx.Graph()

        # 构建表达相关性网络
        expr_sub = self.expr_data.loc[[g for g in pathway_genes if g in self.expr_data.index]]
        corr_matrix = expr_sub.T.corr()
        G = nx.Graph()
        gene_names = list(expr_sub.index)
        for g in gene_names:
            G.add_node(g)
        for i, g1 in enumerate(gene_names):
            for j in range(i + 1, len(gene_names)):
                r = corr_matrix.iloc[i, j]
                if abs(r) > 0.3:
                    G.add_edge(g1, gene_names[j], weight=abs(r))

        # 计算网络拓扑指标
        degree_c = nx.degree_centrality(G) if G.number_of_nodes() > 0 else {}
        betweenness_c = nx.betweenness_centrality(G) if G.number_of_edges() > 0 else {}

        # 计算每个基因参与的通路数量
        gene_pathway_count = {}
        for _, r in rows.iterrows():
            gene = str(r['gene_symbol']).split(',')[0].strip()
            if pd.notna(r.get('gene_pathway')):
                gene_pathway_count[gene] = len(str(r['gene_pathway']).split(';'))

        # 计算 NAE
        results = []
        for gene in pathway_genes:
            if gene not in self.expr_data.index:
                continue

            expr = self.expr_data.loc[gene].values.astype(float)
            expr_var = float(np.var(expr))
            expr_mean = float(np.mean(expr))

            dc = degree_c.get(gene, 0)
            bc = betweenness_c.get(gene, 0)
            pw_count = gene_pathway_count.get(gene, 1)

            # NAE: 归一化各属性后计算信息熵
            attrs = np.array([dc, bc, expr_var / (expr_var + 1), pw_count / (pw_count + 10)])
            attrs = attrs / (attrs.sum() + 1e-10)  # 归一化为概率分布
            # 去掉 0 值避免 log(0)
            attrs = attrs[attrs > 0]
            nae = -np.sum(attrs * np.log2(attrs)) if len(attrs) > 0 else 0

            results.append({
                "基因": gene,
                "NAE评分": round(nae, 4),
                "度中心性": round(dc, 4),
                "介数中心性": round(bc, 4),
                "表达均值": round(expr_mean, 2),
                "表达方差": round(expr_var, 2),
                "参与通路数": pw_count,
            })

        df = pd.DataFrame(results).sort_values("NAE评分", ascending=False)
        return df.reset_index(drop=True), G

    # ==========================================================
    # Layer 4: 基因表达特征 + miRNA 调控
    # ==========================================================

    def get_gene_expression_profile(self, gene_name: str) -> dict:
        """
        Layer 4: 基因表达特征 + miRNA调控关系

        方法 (论文3: Colon Cancer miRNA-Gene):
            Pearson相关性筛选负调控 miRNA (r < -0.3, p < 0.05)
        """
        result = {"gene": gene_name}

        if self.expr_data is None or gene_name not in self.expr_data.index:
            result["error"] = f"{gene_name} 不在表达数据中"
            return result

        expr = self.expr_data.loc[gene_name].values.astype(float)
        result["expr_mean"] = round(float(np.mean(expr)), 2)
        result["expr_std"] = round(float(np.std(expr)), 2)
        result["expr_max"] = round(float(np.max(expr)), 2)
        result["expr_min"] = round(float(np.min(expr)), 2)

        # 表达排名
        all_means = self.expr_data.mean(axis=1).sort_values(ascending=False)
        result["expr_rank"] = list(all_means.index).index(gene_name) + 1 if gene_name in all_means.index else -1
        result["total_genes"] = len(all_means)

        # 按临床分组的表达
        if self.clinical_data is not None:
            common = list(set(self.expr_data.columns) & set(self.clinical_data.index))
            if common and 'Age_Group' in self.clinical_data.columns:
                group_expr = {}
                for sample in common:
                    group = self.clinical_data.loc[sample, 'Age_Group']
                    if pd.notna(group):
                        group_expr.setdefault(str(group), []).append(
                            float(self.expr_data.loc[gene_name, sample]))
                result["group_expression"] = group_expr

        # miRNA 调控 (论文3)
        if self.mirna_data is not None:
            common_samples = list(set(self.expr_data.columns) & set(self.mirna_data.columns))
            if len(common_samples) >= 10:
                gene_expr = self.expr_data.loc[gene_name, common_samples].values.astype(float)
                mirna_correlations = []
                for mirna in self.mirna_data.index[:200]:  # 限制前200个miRNA加速
                    mirna_expr = self.mirna_data.loc[mirna, common_samples].values.astype(float)
                    if np.std(mirna_expr) > 0 and np.std(gene_expr) > 0:
                        r, p = stats.pearsonr(gene_expr, mirna_expr)
                        if r < -0.3 and p < 0.05:
                            mirna_correlations.append({
                                "miRNA": mirna,
                                "相关系数": round(r, 4),
                                "p值": f"{p:.2e}",
                            })
                result["mirna_regulators"] = sorted(
                    mirna_correlations, key=lambda x: x["相关系数"])

        return result

    # ==========================================================
    # 可视化
    # ==========================================================

    def create_sankey_diagram(self, disease_name: str, top_pathways: int = 8,
                              top_genes: int = 3) -> go.Figure:
        """Sankey 图: 疾病 → 通路 → 基因"""
        pw_df = self.get_pathway_ranking(disease_name)
        if pw_df.empty:
            return go.Figure()

        pw_df = pw_df.head(top_pathways)
        labels = [disease_name]
        sources, targets, values, colors = [], [], [], []

        # 疾病 → 通路
        for i, row in pw_df.iterrows():
            pw_name = row["通路"]
            if len(pw_name) > 30:
                pw_name = pw_name[:28] + "..."
            labels.append(pw_name)
            sources.append(0)
            targets.append(len(labels) - 1)
            values.append(row["疾病基因数"])
            colors.append("rgba(102, 126, 234, 0.4)")

        # 通路 → 基因
        gene_set = set()
        for i, row in pw_df.iterrows():
            pw_idx = i + 1  # 通路在labels中的索引
            genes = row["基因列表"].split(", ")[:top_genes]
            for gene in genes:
                if gene not in gene_set:
                    labels.append(gene)
                    gene_set.add(gene)
                gene_idx = labels.index(gene)
                sources.append(pw_idx)
                targets.append(gene_idx)
                values.append(1)
                colors.append("rgba(245, 87, 108, 0.4)")

        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15, thickness=20,
                line=dict(color="black", width=0.5),
                label=labels,
                color=["#667eea"] + ["#f093fb"] * top_pathways + ["#43e97b"] * len(gene_set),
            ),
            link=dict(source=sources, target=targets, value=values, color=colors),
        )])
        fig.update_layout(
            title=dict(text=f"疾病溯源: {disease_name} → 通路 → 基因", x=0.5, font=dict(size=15)),
            height=500, font=dict(size=11))
        return fig

    def create_pathway_bar(self, pw_df: pd.DataFrame, top_n: int = 15) -> go.Figure:
        """通路影响力柱状图"""
        df = pw_df.head(top_n)
        if df.empty:
            return go.Figure()

        fig = go.Figure()
        fig.add_trace(go.Bar(
            y=df["通路"], x=df["影响力评分"],
            orientation='h',
            marker=dict(color=df["影响力评分"], colorscale="Viridis"),
            text=[f'{s:.2f}' for s in df["影响力评分"]],
            textposition='outside',
            hovertemplate='%{y}<br>评分: %{x:.3f}<br>基因数: %{customdata}<extra></extra>',
            customdata=df["疾病基因数"],
        ))
        fig.update_layout(
            title=dict(text="通路影响力排名 (CV加权)", x=0.5, font=dict(size=14)),
            xaxis_title="影响力评分 = 基因数 × 平均CV",
            yaxis=dict(autorange="reversed"),
            height=max(350, len(df) * 28 + 100),
            plot_bgcolor='#fafafa', margin=dict(l=200, r=60, t=50, b=40))
        return fig

    def create_gene_network_plot(self, G: nx.Graph, gene_df: pd.DataFrame) -> go.Figure:
        """基因模块网络图 + NAE 着色"""
        if G.number_of_nodes() == 0:
            return go.Figure()

        pos = nx.spring_layout(G, seed=42, k=1.0)
        nae_scores = dict(zip(gene_df["基因"], gene_df["NAE评分"])) if len(gene_df) > 0 else {}

        # 边
        ex, ey = [], []
        for e in G.edges():
            x0, y0 = pos[e[0]]; x1, y1 = pos[e[1]]
            ex.extend([x0, x1, None]); ey.extend([y0, y1, None])

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=ex, y=ey, mode='lines',
                                  line=dict(width=0.8, color='#ccc'),
                                  hoverinfo='none', showlegend=False))

        # 节点
        nodes = list(G.nodes())
        nx_arr = [pos[n][0] for n in nodes]
        ny_arr = [pos[n][1] for n in nodes]
        nae_vals = [nae_scores.get(n, 0) for n in nodes]
        sizes = [8 + 20 * v / (max(nae_vals) + 1e-6) for v in nae_vals]

        fig.add_trace(go.Scatter(
            x=nx_arr, y=ny_arr, mode='markers+text',
            marker=dict(size=sizes, color=nae_vals, colorscale='Viridis',
                        showscale=True, colorbar=dict(title="NAE"),
                        line=dict(width=0.5, color='white')),
            text=nodes, textposition='top center', textfont=dict(size=9),
            hovertemplate='%{text}<br>NAE: %{marker.color:.3f}<extra></extra>',
            showlegend=False,
        ))

        fig.update_layout(
            title=dict(text="基因相关性网络 (NAE评分着色)", x=0.5),
            xaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            yaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            plot_bgcolor='#fafafa', height=500)
        return fig

    def create_expression_boxplot(self, gene_name: str, profile: dict) -> go.Figure:
        """基因表达箱线图 (按临床分组)"""
        fig = go.Figure()
        group_expr = profile.get("group_expression", {})
        colors = {'young': '#43e97b', 'middle': '#667eea', 'old': '#f5576c'}

        for group, values in sorted(group_expr.items()):
            fig.add_trace(go.Box(
                y=values, name=group,
                marker_color=colors.get(group, '#888'),
                boxmean='sd',
            ))

        fig.update_layout(
            title=dict(text=f"{gene_name} 表达分布 (TCGA-COAD, 按年龄组)", x=0.5),
            yaxis_title="表达量", plot_bgcolor='#fafafa', height=400)
        return fig

    # ==========================================================
    # 溯源报告
    # ==========================================================

    def generate_tracing_report(self, disease: str, pathway: str, gene: str) -> str:
        """生成完整溯源摘要 HTML"""
        overview = self.get_disease_overview(disease)
        pw_df = self.get_pathway_ranking(disease)
        gene_df, G = self.get_gene_module(disease, pathway)
        profile = self.get_gene_expression_profile(gene)

        # 通路排名
        pw_rank = "N/A"
        if not pw_df.empty:
            pw_match = pw_df[pw_df["通路"] == pathway]
            pw_rank = f"#{pw_match.index[0] + 1}" if not pw_match.empty else "未排名"

        # 基因NAE排名
        gene_rank = "N/A"
        if not gene_df.empty:
            gene_match = gene_df[gene_df["基因"] == gene]
            gene_rank = f"#{gene_match.index[0] + 1}" if not gene_match.empty else "未排名"

        # miRNA
        mirna_list = profile.get("mirna_regulators", [])
        mirna_text = f"发现 {len(mirna_list)} 个负调控miRNA" if mirna_list else "未发现显著负调控miRNA"
        top_mirna = ", ".join([m["miRNA"] for m in mirna_list[:3]]) if mirna_list else "无"

        return f"""
        <div style="font-family:system-ui;padding:16px;max-width:700px;margin:auto;">
          <h3 style="text-align:center;">📋 溯源分析报告</h3>

          <div style="background:#f8f9fa;border-radius:10px;padding:16px;margin:8px 0;">
            <h4>🏥 疾病: {disease}</h4>
            <p>分类: {overview.get('category','未知')} | 关联基因: {overview.get('n_genes',0)} 个 |
               涉及通路: {overview.get('n_pathways',0)} 个</p>
          </div>

          <div style="text-align:center;font-size:1.5em;">⬇️</div>

          <div style="background:linear-gradient(135deg,#f093fb20,#f5576c20);border-radius:10px;
                      padding:16px;margin:8px 0;border-left:4px solid #f093fb;">
            <h4>🔬 通路: {pathway}</h4>
            <p>影响力排名: {pw_rank} |
               {'该通路是疾病的核心信号通路之一。' if pw_rank in ['#1','#2','#3'] else '该通路参与疾病的调控过程。'}</p>
          </div>

          <div style="text-align:center;font-size:1.5em;">⬇️</div>

          <div style="background:linear-gradient(135deg,#43e97b20,#38f9d720);border-radius:10px;
                      padding:16px;margin:8px 0;border-left:4px solid #43e97b;">
            <h4>🧬 基因: {gene}</h4>
            <p>NAE排名: {gene_rank} |
               表达均值: {profile.get('expr_mean','N/A')} |
               表达排名: #{profile.get('expr_rank','N/A')}/{profile.get('total_genes','N/A')}</p>
            <p>{mirna_text}。{'核心调控miRNA: ' + top_mirna if mirna_list else ''}</p>
          </div>

          <div style="background:#fff3cd;border-radius:10px;padding:12px;margin:12px 0;
                      border-left:4px solid #ffc107;">
            <b>📝 溯源结论:</b><br/>
            在 {disease} 中，{pathway} 通路{'是核心致病通路（影响力排名前3）' if pw_rank in ['#1','#2','#3'] else '参与疾病调控'}。
            基因 {gene} 在该通路中NAE评分{'排名靠前，是该通路的核心调控基因' if gene_rank in ['#1','#2','#3'] else '具有一定的调控作用'}，
            在TCGA-COAD数据中平均表达量为 {profile.get('expr_mean','N/A')}。
            {f'该基因受到 {len(mirna_list)} 个miRNA的负调控（{top_mirna}），提示存在表观遗传层面的调控机制。' if mirna_list else ''}
          </div>

          <div style="font-size:0.8em;color:#999;margin-top:12px;">
            方法参考: [1] Weighted Gene Networks (CV加权) [2] Node Attribute Entropy [3] miRNA-Gene Correlation
          </div>
        </div>
        """
