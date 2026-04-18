"""
疾病多尺度溯源分析模块 - Disease Multi-Scale Tracing

实现从疾病 → 信号通路 → 基因模块 → 基因表达特征的逐层下钻溯源。

方法:
    - CV加权通路影响力评分
    - NAE (节点属性熵) 基因核心度评分
    - miRNA-基因负相关调控关系映射
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

    # 支持的疾病
    PAPER_DISEASES = {
        "Colorectal cancer": {
            "method": "miRNA-基因整合网络分析",
        },
        "Non-small cell lung cancer": {
            "method": "多组学加权基因网络",
        },
        "Small cell lung cancer": {
            "method": "多组学加权基因网络",
        },
    }

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
        """获取支持的疾病列表"""
        if self.gene_disease is None:
            return []
        available = []
        for d in self.PAPER_DISEASES:
            if d in self.gene_disease['disease_name'].values:
                available.append(d)
        return available

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

        # 疾病信息
        disease_info = self.PAPER_DISEASES.get(disease_name, {})

        return {
            "disease": disease_name,
            "category": category,
            "n_genes": len(genes),
            "n_genes_in_tcga": len(avail_genes),
            "genes": genes,
            "available_genes": avail_genes,
            "n_pathways": len(all_pathways),
            "pathways": sorted(all_pathways),
            "method": disease_info.get("method", ""),
        }

    # ==========================================================
    # Layer 2: 通路影响力排名
    # ==========================================================

    def get_pathway_ranking(self, disease_name: str) -> pd.DataFrame:
        """
        Layer 2: 按通路影响力评分排序

        方法 (CV加权):
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
            score = n_genes * mean_cv  # 加权评分

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

        方法 (NAE):
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

        方法 (miRNA-Gene负相关):
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

        # miRNA 调控
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

    def create_pathway_universe(self, disease_name: str,
                                selected_pathway: str = None) -> go.Figure:
        """
        通路网络图 — 全部 ~350 条通路，按 KEGG 分类着色

        边 = 共享基因 >= 8，布局展开，疾病相关通路放大高亮。
        """
        np.random.seed(42)

        # 从 pathway_data 获取所有通路及其基因 + 分类
        pw_gene_sets = {}
        pw_classes = {}
        SUFFIX = ' - Homo sapiens (human)'
        if self.pathway_data is not None and 'Pathway_Name' in self.pathway_data.columns:
            for _, row in self.pathway_data.iterrows():
                pn = row['Pathway_Name']
                short = pn.replace(SUFFIX, '') if isinstance(pn, str) and pn.endswith(SUFFIX) else pn
                if pd.notna(row.get('Gene', '')):
                    genes = set(g.strip() for g in str(row['Gene']).replace(';', ',').split(',') if g.strip())
                    if genes:
                        pw_gene_sets[short] = genes
                cls = str(row.get('Class', '')) if pd.notna(row.get('Class', '')) else ''
                # 取大类（第一个分号前）
                major = cls.split(';')[0].strip() if cls else 'Other'
                pw_classes[short] = major

        if len(pw_gene_sets) < 5:
            return go.Figure()

        # 疾病相关通路
        pw_df = self.get_pathway_ranking(disease_name)
        disease_pw_set = set(pw_df["通路"].tolist()) if not pw_df.empty else set()

        # 构建图 — 共享基因 >= 8（更稀疏，结构更清晰）
        G = nx.Graph()
        for pn in pw_gene_sets:
            G.add_node(pn)

        pw_names = list(pw_gene_sets.keys())
        for i in range(len(pw_names)):
            for j in range(i + 1, len(pw_names)):
                shared = len(pw_gene_sets[pw_names[i]] & pw_gene_sets[pw_names[j]])
                if shared >= 50:
                    G.add_edge(pw_names[i], pw_names[j], weight=shared)

        # 布局 — k 大一些让节点更分散
        pos = nx.spring_layout(G, seed=42, k=3.5, iterations=80)

        # 如果选中通路，把它移到中心
        if selected_pathway and selected_pathway in pos:
            offset = np.array([0.0, 0.0]) - np.array(pos[selected_pathway])
            for n in pos:
                pos[n] = (pos[n][0] + offset[0], pos[n][1] + offset[1])

        fig = go.Figure()

        # ── 分类颜色 ──
        CLASS_COLORS = {
            'Metabolism':                      '#4ecdc4',
            'Environmental Information Processing': '#667eea',
            'Genetic Information Processing':   '#f093fb',
            'Cellular Processes':               '#ffa07a',
            'Organismal Systems':               '#43e97b',
            'Human Diseases':                   '#f5576c',
            'Drug Development':                 '#f7dc6f',
            'Other':                            '#888888',
        }

        # ── 边 ──
        sel_neighbors = set(G.neighbors(selected_pathway)) if selected_pathway and selected_pathway in G else set()
        highlight_ex, highlight_ey = [], []
        disease_ex, disease_ey = [], []
        bg_ex, bg_ey = [], []

        for u, v, d in G.edges(data=True):
            x0, y0 = pos[u]; x1, y1 = pos[v]
            involves_sel = selected_pathway and (u == selected_pathway or v == selected_pathway)
            both_disease = (u in disease_pw_set and v in disease_pw_set)

            if involves_sel:
                highlight_ex.extend([x0, x1, None])
                highlight_ey.extend([y0, y1, None])
            elif both_disease:
                disease_ex.extend([x0, x1, None])
                disease_ey.extend([y0, y1, None])
            else:
                bg_ex.extend([x0, x1, None])
                bg_ey.extend([y0, y1, None])

        if bg_ex:
            fig.add_trace(go.Scatter(
                x=bg_ex, y=bg_ey, mode='lines',
                line=dict(width=0.3, color='rgba(140,140,180,0.12)'),
                hoverinfo='skip', showlegend=False))
        if disease_ex:
            fig.add_trace(go.Scatter(
                x=disease_ex, y=disease_ey, mode='lines',
                line=dict(width=0.8, color='rgba(102,126,234,0.25)'),
                hoverinfo='skip', showlegend=False))
        if highlight_ex:
            fig.add_trace(go.Scatter(
                x=highlight_ex, y=highlight_ey, mode='lines',
                line=dict(width=2.5, color='rgba(255,215,0,0.5)'),
                hoverinfo='skip', showlegend=False))

        # ── 节点: 按分类分组绘制 ──
        # 先收集所有分类
        all_classes = set(pw_classes.values())
        pw_scores = dict(zip(pw_df["通路"], pw_df["影响力评分"])) if not pw_df.empty else {}
        max_score = max(pw_scores.values()) if pw_scores else 1

        # 按分类画节点（每类一个 trace → 自动图例）
        for cls in sorted(all_classes):
            cls_color = CLASS_COLORS.get(cls, '#888888')
            r, g, b = int(cls_color[1:3], 16), int(cls_color[3:5], 16), int(cls_color[5:7], 16)

            # 该分类下的背景通路
            bg_nodes = [n for n in G.nodes() if pw_classes.get(n, 'Other') == cls
                        and n not in disease_pw_set and n != selected_pathway]
            if bg_nodes:
                fig.add_trace(go.Scatter(
                    x=[pos[n][0] for n in bg_nodes],
                    y=[pos[n][1] for n in bg_nodes],
                    mode='markers',
                    marker=dict(size=4.5, color=f'rgba({r},{g},{b},0.30)',
                                line=dict(width=0.3, color=f'rgba({r},{g},{b},0.12)')),
                    hovertext=[n[:40] for n in bg_nodes], hoverinfo='text',
                    name=cls, showlegend=True,
                    legendgroup=cls,
                ))

            # 该分类下的疾病通路（更大更亮）
            dp_nodes = [n for n in G.nodes() if pw_classes.get(n, 'Other') == cls
                        and n in disease_pw_set and n != selected_pathway]
            if dp_nodes:
                dp_sizes = [6 + 12 * pw_scores.get(n, 0) / max_score for n in dp_nodes]
                fig.add_trace(go.Scatter(
                    x=[pos[n][0] for n in dp_nodes],
                    y=[pos[n][1] for n in dp_nodes],
                    mode='markers+text',
                    marker=dict(size=dp_sizes, color=f'rgba({r},{g},{b},0.7)',
                                line=dict(width=1, color=f'rgba(255,255,255,0.4)')),
                    text=[n[:18] if pw_scores.get(n, 0) / max_score > 0.3 else '' for n in dp_nodes],
                    textposition='top center',
                    textfont=dict(size=7, color=f'rgba({r},{g},{b},0.8)'),
                    hovertext=[f"<b>{n[:40]}</b><br>影响力: {pw_scores.get(n,0):.3f}" for n in dp_nodes],
                    hoverinfo='text',
                    showlegend=False, legendgroup=cls,
                ))

        # 选中通路（金色大发光）
        if selected_pathway and selected_pathway in pos:
            sx, sy = pos[selected_pathway]
            for sz, op in [(60, 0.04), (42, 0.10), (28, 0.22)]:
                fig.add_trace(go.Scatter(
                    x=[sx], y=[sy], mode='markers',
                    marker=dict(size=sz, color=f'rgba(255,215,0,{op})'),
                    hoverinfo='skip', showlegend=False))
            fig.add_trace(go.Scatter(
                x=[sx], y=[sy], mode='markers+text',
                marker=dict(size=20, color='#FFD700', line=dict(width=2, color='white')),
                text=[selected_pathway[:25]], textposition='bottom center',
                textfont=dict(size=11, color='#FFD700', family='Arial Black'),
                hovertext=[f"<b>★ {selected_pathway}</b><br>影响力: {pw_scores.get(selected_pathway,0):.3f}<br>"
                           f"邻居通路: {len(sel_neighbors)}"],
                hoverinfo='text', showlegend=False))

        n_disease = len(disease_pw_set & set(G.nodes()))
        fig.update_layout(
            title=dict(
                text=f"通路网络 — {G.number_of_nodes()} 条通路, {G.number_of_edges()} 条边 | "
                     f"{disease_name}: {n_disease} 条相关通路高亮",
                x=0.5, font=dict(size=13, color='#e0e0e0')),
            legend=dict(
                font=dict(size=9, color='#ccc'),
                bgcolor='rgba(15,16,41,0.7)',
                bordercolor='rgba(100,100,100,0.3)',
                borderwidth=1,
                orientation='v', x=1.01, y=1,
            ),
            xaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            yaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            plot_bgcolor='#0f1029', paper_bgcolor='#0f1029',
            font=dict(color='white'),
            height=700, margin=dict(l=15, r=120, t=50, b=15))
        return fig

    def create_gene_universe(self, disease_name: str,
                             selected_pathway: str = None,
                             selected_gene: str = None,
                             n_background: int = 400) -> go.Figure:
        """
        基因宇宙图 — 几百个基因构成的共表达网络

        取高方差基因 top N 构建共表达网络（|r|>0.6 有边）。
        选中通路的基因在网络中高亮，选中基因金色发光。
        """
        if self.expr_data is None:
            return go.Figure()
        np.random.seed(42)

        # 疾病直接关联基因
        pw_df = self.get_pathway_ranking(disease_name)
        disease_direct_genes = set()
        for _, row in pw_df.iterrows():
            disease_direct_genes.update(row["基因列表"].split(", "))

        # 疾病通路内的全部基因（从 pathway_data 获取，数量大得多）
        pathway_all_genes = set()
        disease_pw_names = set(pw_df["通路"].tolist())
        SUFFIX = ' - Homo sapiens (human)'
        if self.pathway_data is not None and 'Pathway_Name' in self.pathway_data.columns:
            for _, prow in self.pathway_data.iterrows():
                short = prow['Pathway_Name'].replace(SUFFIX, '') if isinstance(prow['Pathway_Name'], str) else ''
                if short in disease_pw_names:
                    if pd.notna(prow.get('Gene', '')):
                        for g in str(prow['Gene']).replace(';', ',').split(','):
                            g = g.strip()
                            if g and g in self.expr_data.index:
                                pathway_all_genes.add(g)

        # 选中通路的基因（子集）
        sel_pw_genes = set()
        if selected_pathway and not pw_df.empty:
            sel_rows = pw_df[pw_df["通路"] == selected_pathway]
            if not sel_rows.empty:
                sel_pw_genes = set(sel_rows.iloc[0]["基因列表"].split(", "))

        # 三层高亮:
        #   1. disease_direct_genes (绿色，最亮)
        #   2. pathway_all_genes - disease_direct_genes (青色，中亮)
        #   3. 其余 (灰色背景)
        # 确保通路基因能出现在网络中
        highlight_genes = disease_direct_genes | pathway_all_genes

        # 选取基因: 高亮基因(采样) + 高方差背景基因
        gene_var = self.expr_data.var(axis=1).sort_values(ascending=False)
        top_genes = list(gene_var.head(n_background).index)
        # 优先加入疾病直接基因
        for g in disease_direct_genes:
            if g in self.expr_data.index and g not in top_genes:
                top_genes.append(g)
        # 再加入通路基因（采样，不然太多）
        pw_only = list(pathway_all_genes - disease_direct_genes - set(top_genes))
        np.random.shuffle(pw_only)
        for g in pw_only[:80]:  # 最多额外加 80 个通路基因
            top_genes.append(g)
        top_genes = top_genes[:n_background + len(disease_direct_genes) + 80]

        # 计算相关性矩阵
        sub_expr = self.expr_data.loc[top_genes].values.astype(float)
        corr_matrix = np.corrcoef(sub_expr)

        # 构建网络 (|r| > 0.6)
        G = nx.Graph()
        for g in top_genes:
            G.add_node(g)

        gene_idx = {g: i for i, g in enumerate(top_genes)}
        for i in range(len(top_genes)):
            for j in range(i + 1, len(top_genes)):
                r = corr_matrix[i, j]
                if abs(r) > 0.5:
                    G.add_edge(top_genes[i], top_genes[j], weight=abs(r))

        # 删除无边的背景基因（保留通路基因即使无边）
        isolates = [n for n in nx.isolates(G) if n not in highlight_genes]
        keep_iso = isolates[:30]  # 保留一些孤立点做背景
        G.remove_nodes_from(isolates[80:])

        pos = nx.spring_layout(G, seed=42, k=1.2, iterations=50)

        # 把选中基因移到中心
        if selected_gene and selected_gene in pos:
            offset = np.array([0.0, 0.0]) - np.array(pos[selected_gene])
            for n in pos:
                pos[n] = (pos[n][0] + offset[0], pos[n][1] + offset[1])

        fig = go.Figure()

        # ── 边 ──
        sel_neighbors = set()
        if selected_gene and selected_gene in G:
            sel_neighbors = set(G.neighbors(selected_gene))

        highlight_ex, highlight_ey = [], []  # 选中基因的边
        pw_ex, pw_ey = [], []                # 高亮基因之间的边
        bg_ex, bg_ey = [], []                # 背景边

        for u, v, d in G.edges(data=True):
            x0, y0 = pos[u]; x1, y1 = pos[v]
            involves_sel = selected_gene and (u == selected_gene or v == selected_gene)
            both_highlight = (u in highlight_genes and v in highlight_genes)

            if involves_sel:
                highlight_ex.extend([x0, x1, None])
                highlight_ey.extend([y0, y1, None])
            elif both_highlight:
                pw_ex.extend([x0, x1, None])
                pw_ey.extend([y0, y1, None])
            else:
                bg_ex.extend([x0, x1, None])
                bg_ey.extend([y0, y1, None])

        if bg_ex:
            fig.add_trace(go.Scatter(
                x=bg_ex, y=bg_ey, mode='lines',
                line=dict(width=0.3, color='rgba(140,140,180,0.12)'),
                hoverinfo='skip', showlegend=False))
        if pw_ex:
            fig.add_trace(go.Scatter(
                x=pw_ex, y=pw_ey, mode='lines',
                line=dict(width=1.2, color='rgba(67,233,123,0.3)'),
                hoverinfo='skip', showlegend=False))
        if highlight_ex:
            fig.add_trace(go.Scatter(
                x=highlight_ex, y=highlight_ey, mode='lines',
                line=dict(width=2.5, color='rgba(255,215,0,0.5)'),
                hoverinfo='skip', showlegend=False))

        # ── 节点: 4类 ──
        bg_x, bg_y, bg_hover = [], [], []           # 完全无关
        pwg_x, pwg_y, pwg_labels = [], [], []       # 通路内基因（青色中亮）
        direct_x, direct_y, direct_labels = [], [], []  # 疾病直接基因（绿色亮）
        coex_x, coex_y, coex_labels = [], [], []    # 共表达邻居（橙色）

        for n in G.nodes():
            x, y = pos[n]
            if n == selected_gene:
                continue
            elif n in disease_direct_genes:
                direct_x.append(x); direct_y.append(y)
                direct_labels.append(n)
            elif n in pathway_all_genes:
                pwg_x.append(x); pwg_y.append(y)
                pwg_labels.append(n)
            elif n in sel_neighbors:
                coex_x.append(x); coex_y.append(y)
                coex_labels.append(n)
            else:
                bg_x.append(x); bg_y.append(y)
                bg_hover.append(n)

        # 背景基因（灰色）
        if bg_x:
            fig.add_trace(go.Scatter(
                x=bg_x, y=bg_y, mode='markers',
                marker=dict(size=3.5, color='rgba(150,150,190,0.30)',
                            line=dict(width=0.3, color='rgba(150,150,190,0.12)')),
                hovertext=bg_hover, hoverinfo='text', showlegend=False))

        # 共表达基因（橙色）
        if coex_x:
            fig.add_trace(go.Scatter(
                x=coex_x, y=coex_y, mode='markers',
                marker=dict(size=6, color='rgba(255,200,100,0.5)',
                            line=dict(width=0.5, color='rgba(255,200,100,0.25)')),
                hovertext=coex_labels, hoverinfo='text', showlegend=False))

        # 通路内基因（青色，中等亮度）
        if pwg_x:
            fig.add_trace(go.Scatter(
                x=pwg_x, y=pwg_y, mode='markers',
                marker=dict(size=5.5, color='rgba(78,205,196,0.50)',
                            line=dict(width=0.5, color='rgba(78,205,196,0.2)')),
                hovertext=pwg_labels, hoverinfo='text', showlegend=False))

        # 疾病直接基因（绿色亮 + 文字）
        if direct_x:
            fig.add_trace(go.Scatter(
                x=direct_x, y=direct_y, mode='markers',
                marker=dict(size=[14] * len(direct_x), color='rgba(67,233,123,0.08)'),
                hoverinfo='skip', showlegend=False))
            fig.add_trace(go.Scatter(
                x=direct_x, y=direct_y, mode='markers+text',
                marker=dict(size=9, color='#43e97b',
                            line=dict(width=1, color='rgba(255,255,255,0.4)')),
                text=direct_labels, textposition='top center',
                textfont=dict(size=8, color='rgba(67,233,123,0.8)'),
                hovertext=direct_labels, hoverinfo='text', showlegend=False))

        # 选中基因（金色 3 层光晕）
        if selected_gene and selected_gene in pos:
            sx, sy = pos[selected_gene]
            for sz, op in [(55, 0.04), (40, 0.10), (28, 0.22)]:
                fig.add_trace(go.Scatter(
                    x=[sx], y=[sy], mode='markers',
                    marker=dict(size=sz, color=f'rgba(255,215,0,{op})'),
                    hoverinfo='skip', showlegend=False))
            fig.add_trace(go.Scatter(
                x=[sx], y=[sy], mode='markers+text',
                marker=dict(size=18, color='#FFD700', symbol='star',
                            line=dict(width=2, color='white')),
                text=[f"★ {selected_gene}"], textposition='bottom center',
                textfont=dict(size=12, color='#FFD700', family='Arial Black'),
                hovertext=[f"<b>★ {selected_gene}</b><br>共表达邻居: {len(sel_neighbors)}<br>疾病基因: {'是' if selected_gene in highlight_genes else '否'}"],
                hoverinfo='text', showlegend=False))

        n_direct = len(disease_direct_genes & set(G.nodes()))
        n_pw_genes = len(pathway_all_genes & set(G.nodes()))
        fig.update_layout(
            title=dict(
                text=f"基因网络 — {G.number_of_nodes()} 个基因 | {disease_name}: {n_direct} 疾病基因(绿) + {n_pw_genes} 通路基因(青)",
                x=0.5, font=dict(size=13, color='#e0e0e0')),
            annotations=[dict(x=0.01, y=0.99, xref='paper', yref='paper', showarrow=False,
                              text="灰=背景  青=通路基因  绿=疾病基因  橙=共表达  金★=选中",
                              font=dict(size=9, color='rgba(200,200,200,0.5)'))],
            xaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            yaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            plot_bgcolor='#0d1117', paper_bgcolor='#0d1117',
            font=dict(color='white'),
            height=700, margin=dict(l=15, r=15, t=50, b=15))
        return fig

    def create_sankey_diagram(self, disease_name: str, selected_pathway: str = None,
                              selected_gene: str = None,
                              top_pathways: int = 20, top_genes: int = 6) -> go.Figure:
        """Sankey 溯源图: 疾病 → 通路 → 基因 — 选中路径高亮，其余灰暗"""
        pw_df = self.get_pathway_ranking(disease_name)
        if pw_df.empty:
            return go.Figure()

        pw_df = pw_df.head(top_pathways)
        labels = [f"🏥 {disease_name}"]
        sources, targets, values, link_colors = [], [], [], []

        # 颜色方案
        highlight_colors = ['#667eea', '#764ba2', '#f093fb', '#f5576c', '#43e97b',
                            '#ffa07a', '#4ecdc4', '#bb8fce', '#f7dc6f', '#45b7d1']
        BG_LINK = "rgba(160,160,160,0.06)"
        BG_NODE = "#555"

        # 找到选中通路的index
        selected_pw_idx = None
        selected_pw_genes = set()

        # 疾病 → 通路
        for i, row in pw_df.iterrows():
            pw_name = row["通路"]
            pw_short = pw_name[:35] + ("..." if len(pw_name) > 35 else "")
            labels.append(f"🔬 {pw_short}")
            sources.append(0)
            targets.append(len(labels) - 1)
            values.append(max(row["疾病基因数"], 1))

            is_selected = (selected_pathway and pw_name == selected_pathway)
            if is_selected:
                selected_pw_idx = i + 1
                selected_pw_genes = set(row["基因列表"].split(", "))
                c = highlight_colors[i % len(highlight_colors)]
                link_colors.append(f"rgba({int(c[1:3],16)},{int(c[3:5],16)},{int(c[5:7],16)},0.65)")
            else:
                link_colors.append(BG_LINK)

        # 通路 → 基因
        gene_set = set()
        for i, row in pw_df.iterrows():
            pw_idx = i + 1
            pw_name = row["通路"]
            is_selected_pw = (selected_pathway and pw_name == selected_pathway)
            genes = row["基因列表"].split(", ")[:top_genes]
            for gene in genes:
                if gene not in gene_set:
                    labels.append(f"🧬 {gene}")
                    gene_set.add(gene)
                gene_idx = labels.index(f"🧬 {gene}")
                sources.append(pw_idx)
                targets.append(gene_idx)
                values.append(1)
                if is_selected_pw:
                    link_colors.append("rgba(67, 233, 123, 0.55)")
                else:
                    link_colors.append(BG_LINK)

        # 节点颜色：选中通路彩色，选中基因金色，其余灰色
        n_pw = len(pw_df)
        node_colors = ["#667eea"]  # 疾病节点
        for i, row in pw_df.iterrows():
            is_sel = (selected_pathway and row["通路"] == selected_pathway)
            node_colors.append(highlight_colors[i % len(highlight_colors)] if is_sel else BG_NODE)

        for label in labels[n_pw + 1:]:
            gene = label.replace("🧬 ", "")
            if selected_gene and gene == selected_gene:
                node_colors.append("#FFD700")  # 金色
            elif gene in selected_pw_genes:
                node_colors.append("#43e97b")  # 绿色
            else:
                node_colors.append(BG_NODE)

        fig = go.Figure(data=[go.Sankey(
            arrangement="snap",
            node=dict(
                pad=15, thickness=25,
                line=dict(color="rgba(255,255,255,0.5)", width=1),
                label=labels,
                color=node_colors,
            ),
            link=dict(source=sources, target=targets, value=values, color=link_colors),
        )])

        subtitle = ""
        if selected_pathway:
            subtitle = f" | 聚焦: {selected_pathway[:30]}"
        fig.update_layout(
            title=dict(text=f"🔬 疾病多尺度溯源: {disease_name}{subtitle}", x=0.5,
                       font=dict(size=16, color='#e0e0e0')),
            height=600, font=dict(size=11, color='#ccc'),
            paper_bgcolor='#0f0f23', plot_bgcolor='#0f0f23')
        return fig

    def create_pathway_network(self, disease_name: str, selected_pathway: str = None,
                               top_n: int = 40) -> go.Figure:
        """通路星云图：全部通路铺开，选中通路高亮居中，关联通路中亮，无关灰暗"""
        pw_df = self.get_pathway_ranking(disease_name)
        if pw_df.empty or len(pw_df) < 2:
            return go.Figure()

        pw_df = pw_df.head(top_n)

        # 构建通路-通路网络（共享基因为边）
        pw_gene_sets = {}
        for _, row in pw_df.iterrows():
            pw_gene_sets[row["通路"]] = set(row["基因列表"].split(", "))

        G = nx.Graph()
        pw_names = list(pw_gene_sets.keys())
        scores = dict(zip(pw_df["通路"], pw_df["影响力评分"]))
        for pw in pw_names:
            G.add_node(pw, score=scores.get(pw, 0))

        for i in range(len(pw_names)):
            for j in range(i + 1, len(pw_names)):
                shared = pw_gene_sets[pw_names[i]] & pw_gene_sets[pw_names[j]]
                if shared:
                    G.add_edge(pw_names[i], pw_names[j], weight=len(shared),
                               shared_genes=", ".join(list(shared)[:10]))

        pos = nx.spring_layout(G, seed=42, k=2.5, iterations=80)

        # 如果有选中通路，把它移到中心
        if selected_pathway and selected_pathway in pos:
            sel_pos = pos[selected_pathway]
            offset = np.array([0.0, 0.0]) - np.array(sel_pos)
            for n in pos:
                pos[n] = (pos[n][0] + offset[0], pos[n][1] + offset[1])

        # 分类节点
        connected_to_selected = set()
        if selected_pathway and selected_pathway in G:
            connected_to_selected = set(G.neighbors(selected_pathway))

        fig = go.Figure()

        # === 边（批量绘制，减少trace数量）===
        sel_edge_x, sel_edge_y = [], []
        con_edge_x, con_edge_y = [], []
        bg_edge_x, bg_edge_y = [], []

        for u, v, d in G.edges(data=True):
            x0, y0 = pos[u]; x1, y1 = pos[v]
            involves_selected = selected_pathway and (u == selected_pathway or v == selected_pathway)
            involves_connected = (u in connected_to_selected or v in connected_to_selected)

            if involves_selected:
                sel_edge_x.extend([x0, x1, None])
                sel_edge_y.extend([y0, y1, None])
            elif involves_connected:
                con_edge_x.extend([x0, x1, None])
                con_edge_y.extend([y0, y1, None])
            else:
                if not selected_pathway:
                    bg_edge_x.extend([x0, x1, None])
                    bg_edge_y.extend([y0, y1, None])

        # 背景边（最淡）
        if bg_edge_x:
            fig.add_trace(go.Scatter(
                x=bg_edge_x, y=bg_edge_y, mode='lines',
                line=dict(width=0.5, color='rgba(80,80,80,0.08)'),
                hoverinfo='skip', showlegend=False,
            ))
        # 关联边
        if con_edge_x:
            fig.add_trace(go.Scatter(
                x=con_edge_x, y=con_edge_y, mode='lines',
                line=dict(width=1.2, color='rgba(102,126,234,0.25)'),
                hoverinfo='skip', showlegend=False,
            ))
        # 选中通路的边（最亮）
        if sel_edge_x:
            fig.add_trace(go.Scatter(
                x=sel_edge_x, y=sel_edge_y, mode='lines',
                line=dict(width=3, color='rgba(255,215,0,0.6)'),
                hoverinfo='skip', showlegend=False,
            ))

        # === 节点分三类渲染 ===
        max_s = max(scores.values()) if scores else 1

        # 1) 背景节点（灰色小点）
        bg_x, bg_y, bg_sizes, bg_hovers = [], [], [], []
        # 2) 关联节点（中等彩色）
        cn_x, cn_y, cn_sizes, cn_colors, cn_labels, cn_hovers = [], [], [], [], [], []
        # 3) 选中节点（金色大发光）
        sel_x, sel_y = [], []

        for n in G.nodes():
            x, y = pos[n]
            s = scores.get(n, 0)

            if selected_pathway and n == selected_pathway:
                sel_x.append(x); sel_y.append(y)
            elif n in connected_to_selected or not selected_pathway:
                cn_x.append(x); cn_y.append(y)
                cn_sizes.append(18 + 30 * s / max_s)
                cn_colors.append(s)
                lbl = n[:22] + "…" if len(n) > 22 else n
                cn_labels.append(lbl)
                cn_hovers.append(f"<b>{n}</b><br>影响力: {s:.3f}")
            else:
                bg_x.append(x); bg_y.append(y)
                bg_sizes.append(5 + 5 * s / max_s)
                bg_hovers.append(f"{n}<br>影响力: {s:.3f}")

        # 画背景节点
        if bg_x:
            fig.add_trace(go.Scatter(
                x=bg_x, y=bg_y, mode='markers',
                marker=dict(size=bg_sizes, color='rgba(100,100,100,0.25)',
                            line=dict(width=0)),
                hovertext=bg_hovers, hoverinfo='text', showlegend=False,
            ))

        # 画关联节点光晕
        if cn_x:
            fig.add_trace(go.Scatter(
                x=cn_x, y=cn_y, mode='markers',
                marker=dict(size=[s + 12 for s in cn_sizes],
                            color=cn_colors, colorscale='Viridis', opacity=0.12),
                hoverinfo='none', showlegend=False,
            ))
            # 关联节点实体
            fig.add_trace(go.Scatter(
                x=cn_x, y=cn_y, mode='markers+text',
                marker=dict(size=cn_sizes, color=cn_colors, colorscale='Viridis',
                            showscale=not selected_pathway,
                            colorbar=dict(title="影响力", tickfont=dict(color='white'),
                                          title_font=dict(color='white')) if not selected_pathway else None,
                            line=dict(width=1.5, color='rgba(255,255,255,0.6)')),
                text=cn_labels, textposition='bottom center',
                textfont=dict(size=9, color='rgba(200,200,200,0.8)'),
                hovertext=cn_hovers, hoverinfo='text', showlegend=False,
            ))

        # 画选中节点（3层发光）
        if sel_x and selected_pathway:
            sel_score = scores.get(selected_pathway, 0)
            # 第3层光晕（最大最淡）
            fig.add_trace(go.Scatter(
                x=sel_x, y=sel_y, mode='markers',
                marker=dict(size=[75], color='rgba(255,215,0,0.06)'),
                hoverinfo='none', showlegend=False,
            ))
            # 第2层光晕
            fig.add_trace(go.Scatter(
                x=sel_x, y=sel_y, mode='markers',
                marker=dict(size=[58], color='rgba(255,215,0,0.12)'),
                hoverinfo='none', showlegend=False,
            ))
            # 第1层光晕
            fig.add_trace(go.Scatter(
                x=sel_x, y=sel_y, mode='markers',
                marker=dict(size=[45], color='rgba(255,215,0,0.25)'),
                hoverinfo='none', showlegend=False,
            ))
            # 实体节点
            fig.add_trace(go.Scatter(
                x=sel_x, y=sel_y, mode='markers+text',
                marker=dict(size=[35], color='#FFD700',
                            line=dict(width=3, color='white')),
                text=[selected_pathway[:28]],
                textposition='bottom center',
                textfont=dict(size=13, color='#FFD700', family='Arial Black'),
                hovertext=[f"<b>★ {selected_pathway}</b><br>影响力: {sel_score:.3f}<br>关联通路: {len(connected_to_selected)}"],
                hoverinfo='text', showlegend=False,
            ))

        n_connected = len(connected_to_selected) if selected_pathway else 0
        title_text = f"通路星云图 — {disease_name}"
        if selected_pathway:
            title_text += f" | ★ {selected_pathway[:25]} ({n_connected} 关联通路)"

        fig.update_layout(
            title=dict(text=title_text, x=0.5, font=dict(size=14, color='#e0e0e0')),
            xaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            yaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            plot_bgcolor='#0a0a2e', paper_bgcolor='#0a0a2e',
            font=dict(color='white'),
            height=600, margin=dict(l=20, r=20, t=50, b=20))
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

    def create_gene_network_plot(self, G: nx.Graph, gene_df: pd.DataFrame,
                                selected_gene: str = None,
                                neighbor_genes: list = None) -> go.Figure:
        """基因宇宙图 — 选中基因金色发光，通路基因彩色，邻居基因灰色外围"""
        if G.number_of_nodes() == 0:
            return go.Figure()

        # 将邻居基因加入图中（灰色虚线连接）
        pathway_genes = set(G.nodes())
        if neighbor_genes:
            for ng in neighbor_genes:
                if ng not in G:
                    G.add_node(ng)
            # 邻居与通路基因之间添加虚连接（用于布局）
            for ng in neighbor_genes:
                if ng in G:
                    # 找一个通路内基因连上
                    for pg in pathway_genes:
                        if pg in G:
                            G.add_edge(ng, pg, weight=0.05, is_neighbor=True)
                            break

        pos = nx.spring_layout(G, seed=42, k=1.8, iterations=60)
        nae_scores = dict(zip(gene_df["基因"], gene_df["NAE评分"])) if len(gene_df) > 0 else {}
        neighbor_set = set(neighbor_genes) if neighbor_genes else set()

        fig = go.Figure()

        # === 边（批量绘制）===
        nb_edge_x, nb_edge_y = [], []
        sel_edge_x, sel_edge_y = [], []
        norm_edge_x, norm_edge_y = [], []

        for u, v, d in G.edges(data=True):
            x0, y0 = pos[u]; x1, y1 = pos[v]
            is_nb_edge = d.get('is_neighbor', False) or u in neighbor_set or v in neighbor_set

            if is_nb_edge:
                nb_edge_x.extend([x0, x1, None])
                nb_edge_y.extend([y0, y1, None])
            else:
                involves_sel = selected_gene and (u == selected_gene or v == selected_gene)
                if involves_sel:
                    sel_edge_x.extend([x0, x1, None])
                    sel_edge_y.extend([y0, y1, None])
                else:
                    norm_edge_x.extend([x0, x1, None])
                    norm_edge_y.extend([y0, y1, None])

        # 邻居边（灰色虚线）
        if nb_edge_x:
            fig.add_trace(go.Scatter(
                x=nb_edge_x, y=nb_edge_y, mode='lines',
                line=dict(width=0.5, color='rgba(120,120,120,0.15)', dash='dot'),
                hoverinfo='skip', showlegend=False,
            ))
        # 普通通路内边（白色半透明）
        if norm_edge_x:
            fig.add_trace(go.Scatter(
                x=norm_edge_x, y=norm_edge_y, mode='lines',
                line=dict(width=1.5, color='rgba(255,255,255,0.25)'),
                hoverinfo='skip', showlegend=False,
            ))
        # 选中基因的边（金色）
        if sel_edge_x:
            fig.add_trace(go.Scatter(
                x=sel_edge_x, y=sel_edge_y, mode='lines',
                line=dict(width=3, color='rgba(255,215,0,0.6)'),
                hoverinfo='skip', showlegend=False,
            ))

        # === 节点分三类 ===
        max_nae = max(nae_scores.values()) if nae_scores else 1

        # 1) 邻居基因（灰色小点外围）
        nb_x, nb_y, nb_hover = [], [], []
        # 2) 通路基因（Plasma彩色）
        pw_x, pw_y, pw_sizes, pw_nae, pw_labels, pw_degrees = [], [], [], [], [], []
        # 3) 选中基因
        sel_x, sel_y, sel_nae_val = [], [], 0

        for n in G.nodes():
            x, y = pos[n]
            if selected_gene and n == selected_gene:
                sel_x.append(x); sel_y.append(y)
                sel_nae_val = nae_scores.get(n, 0)
            elif n in neighbor_set:
                nb_x.append(x); nb_y.append(y)
                nb_hover.append(f"{n} (邻居通路基因)")
            else:
                pw_x.append(x); pw_y.append(y)
                nae = nae_scores.get(n, 0)
                pw_nae.append(nae)
                pw_sizes.append(12 + 25 * nae / max_nae)
                pw_labels.append(n)
                pw_degrees.append(G.degree(n))

        # 画邻居基因
        if nb_x:
            fig.add_trace(go.Scatter(
                x=nb_x, y=nb_y, mode='markers',
                marker=dict(size=6, color='rgba(100,100,100,0.3)',
                            symbol='diamond', line=dict(width=0)),
                hovertext=nb_hover, hoverinfo='text',
                showlegend=False, name='邻居基因',
            ))

        # 画通路基因光晕
        if pw_x:
            fig.add_trace(go.Scatter(
                x=pw_x, y=pw_y, mode='markers',
                marker=dict(size=[s + 15 for s in pw_sizes],
                            color=pw_nae, colorscale='Plasma', opacity=0.15),
                hoverinfo='none', showlegend=False,
            ))
            # 通路基因实体
            fig.add_trace(go.Scatter(
                x=pw_x, y=pw_y, mode='markers+text',
                marker=dict(
                    size=pw_sizes, color=pw_nae, colorscale='Plasma',
                    showscale=True,
                    colorbar=dict(title="NAE", tickfont=dict(color='white'),
                                  title_font=dict(color='white')),
                    line=dict(width=1.5, color='rgba(255,255,255,0.6)'),
                ),
                text=pw_labels, textposition='top center',
                textfont=dict(size=9, color='rgba(220,220,220,0.8)', family='monospace'),
                hovertemplate='<b>%{text}</b><br>NAE: %{marker.color:.3f}<br>度数: %{customdata}<extra></extra>',
                customdata=pw_degrees, showlegend=False,
            ))

        # 画选中基因（3层发光 + 金色）
        if sel_x and selected_gene:
            # 第3层（最大最淡）
            fig.add_trace(go.Scatter(
                x=sel_x, y=sel_y, mode='markers',
                marker=dict(size=[65], color='rgba(255,215,0,0.05)'),
                hoverinfo='none', showlegend=False,
            ))
            # 第2层
            fig.add_trace(go.Scatter(
                x=sel_x, y=sel_y, mode='markers',
                marker=dict(size=[50], color='rgba(255,215,0,0.12)'),
                hoverinfo='none', showlegend=False,
            ))
            # 第1层
            fig.add_trace(go.Scatter(
                x=sel_x, y=sel_y, mode='markers',
                marker=dict(size=[38], color='rgba(255,215,0,0.25)'),
                hoverinfo='none', showlegend=False,
            ))
            # 实体
            fig.add_trace(go.Scatter(
                x=sel_x, y=sel_y, mode='markers+text',
                marker=dict(size=[30], color='#FFD700',
                            line=dict(width=3, color='white')),
                text=[f"★ {selected_gene}"],
                textposition='bottom center',
                textfont=dict(size=13, color='#FFD700', family='Arial Black'),
                hovertext=[f"<b>★ {selected_gene}</b><br>NAE: {sel_nae_val:.3f}<br>度数: {G.degree(selected_gene)}"],
                hoverinfo='text', showlegend=False,
            ))

        title = "基因宇宙图 (NAE × Plasma)"
        if selected_gene:
            title += f" | ★ {selected_gene}"
        if neighbor_genes:
            title += f" | {len(neighbor_set)} 邻居基因"

        fig.update_layout(
            title=dict(text=title, x=0.5, font=dict(size=14, color='#e0e0e0')),
            xaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            yaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            plot_bgcolor='#0d1117', paper_bgcolor='#0d1117',
            font=dict(color='white'),
            height=600, margin=dict(l=20, r=20, t=50, b=20))
        return fig

    def create_expression_boxplot(self, gene_name: str, profile: dict,
                                  context_genes: list = None) -> go.Figure:
        """基因表达箱线图 — 选中基因彩色突出，同通路基因灰色背景"""
        fig = go.Figure()
        group_expr = profile.get("group_expression", {})
        colors = {'young': '#43e97b', 'middle': '#667eea', 'old': '#f5576c'}

        # 先画背景：同通路其他基因的表达分布（灰色小提琴）
        if context_genes and self.expr_data is not None:
            for cg in context_genes[:12]:  # 最多12个背景基因
                if cg == gene_name or cg not in self.expr_data.index:
                    continue
                vals = self.expr_data.loc[cg].dropna().values
                if len(vals) > 5:
                    fig.add_trace(go.Violin(
                        y=vals, name=cg,
                        line=dict(color='rgba(130,130,130,0.15)', width=0.5),
                        fillcolor='rgba(130,130,130,0.05)',
                        meanline=dict(visible=False),
                        showlegend=False, opacity=0.3,
                        scalemode='width', width=0.8,
                    ))

        # 再画选中基因（彩色box，叠在最上面）
        for group, values in sorted(group_expr.items()):
            fig.add_trace(go.Box(
                y=values, name=f"★ {gene_name} ({group})",
                marker_color=colors.get(group, '#888'),
                boxmean='sd', line=dict(width=2),
                fillcolor=colors.get(group, '#888'),
                opacity=0.9,
            ))

        bg_label = f" | 灰色背景: 同通路 {len(context_genes) if context_genes else 0} 个基因" if context_genes else ""
        fig.update_layout(
            title=dict(text=f"★ {gene_name} 表达分布 (TCGA-COAD){bg_label}", x=0.5,
                       font=dict(size=14, color='#e0e0e0')),
            yaxis_title="表达量",
            plot_bgcolor='#0d1117', paper_bgcolor='#0d1117',
            font=dict(color='#ccc'),
            height=450,
            xaxis=dict(tickfont=dict(color='#aaa')),
            yaxis=dict(tickfont=dict(color='#aaa'), gridcolor='rgba(100,100,100,0.2)'),
        )
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
