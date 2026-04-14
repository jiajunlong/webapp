"""
模型库模块 - Model Library Module
定义跨尺度仿真模型库中的所有模型目录和展示功能
"""

import pandas as pd


# ==================== 模型目录 ====================

MODEL_CATALOG = [
    {
        "id": "gene_interaction",
        "name": "基因互作网络模型",
        "name_en": "Gene Interaction Network Model",
        "scale": "分子尺度",
        "scale_en": "Molecular Scale",
        "scale_icon": "🧬",
        "scale_color": "#667eea",
        "description": "构建基因之间的无向相互作用网络，揭示基因在疾病中的协作关系。基于基因-疾病关联数据库，提取特定疾病相关基因及其互作关系，通过网络拓扑分析识别关键Hub基因。",
        "algorithm": "基于图论的网络构建 + 度中心性分析",
        "parameters": {
            "疾病名称": "目标疾病（如肺癌、乳腺癌等）",
            "网络类型": "互作网络（无向图）"
        },
        "input": "基因-疾病关联数据库（gene_disease.tsv）",
        "output": "基因互作网络图、Hub基因列表、网络拓扑指标（密度、聚类系数、平均度）",
        "metrics": ["网络密度", "聚类系数", "平均度", "度中心性"]
    },
    {
        "id": "gene_regulation",
        "name": "基因调控网络模型",
        "name_en": "Gene Regulatory Network Model",
        "scale": "分子尺度",
        "scale_en": "Molecular Scale",
        "scale_icon": "🧬",
        "scale_color": "#667eea",
        "description": "构建基因之间的有向调控网络，展示基因调控级联关系。调控基因（regulator）对目标基因（target）的调控方向通过有向边表示，用于发现关键调控因子。",
        "algorithm": "有向图构建 + 调控级联分析",
        "parameters": {
            "疾病名称": "目标疾病",
            "网络类型": "调控网络（有向图）"
        },
        "input": "基因-疾病关联数据库（含调控关系）",
        "output": "基因调控网络图（有向）、调控因子排名、调控级联路径",
        "metrics": ["入度", "出度", "调控级联深度"]
    },
    {
        "id": "is_coefficient",
        "name": "IS影响力系数模型",
        "name_en": "Influence Score (IS) Coefficient Model",
        "scale": "分子尺度",
        "scale_en": "Molecular Scale",
        "scale_icon": "🧬",
        "scale_color": "#667eea",
        "description": "量化特定生物学通路在疾病中的影响力。综合考虑疾病基因在通路中的覆盖率、通路基因在疾病中的占比、以及基因间的连接密度，计算通路对疾病的影响分数（IS系数）。",
        "algorithm": "IS = f(覆盖率, 占比, 连接密度)",
        "parameters": {
            "疾病名称": "目标疾病",
            "通路列表": "待评估的生物学通路集合"
        },
        "input": "基因-疾病网络 + 基因-通路映射",
        "output": "各通路IS系数排名、柱状图可视化、统计摘要",
        "metrics": ["IS系数", "覆盖率", "基因重叠度"]
    },
    {
        "id": "mrnetb",
        "name": "MRNetB网络推断模型",
        "name_en": "MRNetB Network Inference Model",
        "scale": "细胞/组织尺度",
        "scale_en": "Cellular/Tissue Scale",
        "scale_icon": "🔬",
        "scale_color": "#f093fb",
        "description": "基于互信息（Mutual Information）的基因表达网络反向工程推断算法。从TCGA-COAD真实临床基因表达数据出发，通过MRNetB算法推断基因之间的功能关联网络，支持按年龄、性别、疾病阶段分组分析。",
        "algorithm": "MRNetB: weight(i,j) = MI(i,j) - max(min(MI(i,:), MI(j,:)))",
        "parameters": {
            "表达数据": "基因/miRNA表达矩阵",
            "最大特征数": "参与计算的特征上限（默认200）",
            "特征选择方法": "按方差/均值/随机采样",
            "最小权重阈值": "过滤弱连接的阈值"
        },
        "input": "TCGA-COAD基因表达数据（filtered_hiseq_data.csv）/ miRNA表达数据",
        "output": "基因功能关联网络、边权重列表、网络拓扑统计",
        "metrics": ["互信息", "MRNetB权重", "网络密度", "聚类系数"]
    },
    {
        "id": "mfe",
        "name": "网络流熵(MFE)模型",
        "name_en": "Network Flow Entropy (MFE) Model",
        "scale": "细胞/组织尺度",
        "scale_en": "Cellular/Tissue Scale",
        "scale_icon": "🔬",
        "scale_color": "#f093fb",
        "description": "基于信息论的网络流熵度量模型。将基因表达值与网络拓扑结构结合，计算每个基因节点的信息流熵，用于衡量基因在网络中的信息传递能力和功能重要性。",
        "algorithm": "MFE(i) = -Σ(π_i × p_ij × log(π_i × p_ij))",
        "parameters": {
            "表达值": "基因表达量向量",
            "邻接矩阵": "网络结构矩阵"
        },
        "input": "基因表达网络 + 表达值",
        "output": "各节点MFE值、信息流分布",
        "metrics": ["网络流熵", "信息流分布"]
    },
    {
        "id": "sis_epidemic",
        "name": "SIS传播动力学模型",
        "name_en": "SIS Epidemic Dynamics Model",
        "scale": "群体尺度",
        "scale_en": "Population Scale",
        "scale_icon": "👥",
        "scale_color": "#f5576c",
        "description": "基于社区结构网络的SIS（易感-感染-易感）传播动力学模型。构建具有社区划分的社交网络，模拟疾病或信息在网络中的传播过程，分析社区结构对传播动态的影响。",
        "algorithm": "SIS: P(S→I) = 1-(1-β)^(感染邻居数), P(I→S) = γ",
        "parameters": {
            "节点总数(N)": "网络规模（50-200）",
            "社区数量(c)": "网络社区个数（2-10）",
            "感染率(β)": "每步感染概率（0.01-0.2）",
            "恢复率(γ)": "每步恢复概率（0.05-0.5）",
            "初始感染比例": "初始感染节点比例",
            "仿真步数": "模拟时间步长"
        },
        "input": "网络参数 + 传播参数",
        "output": "感染密度时间序列、社区网络结构图、感染状态快照",
        "metrics": ["感染密度", "传播速率", "稳态感染比例"]
    }
]


def get_model_catalog():
    """获取模型目录"""
    return MODEL_CATALOG


def create_model_cards_html(catalog=None):
    """
    生成模型库卡片HTML

    返回:
        HTML字符串，2列CSS Grid布局的模型卡片
    """
    if catalog is None:
        catalog = MODEL_CATALOG

    cards_html = ""
    for model in catalog:
        # 参数列表HTML
        params_html = ""
        for pname, pdesc in model["parameters"].items():
            params_html += f"<li><strong>{pname}</strong>：{pdesc}</li>"

        cards_html += f"""
        <div style="background: white; border-radius: 12px; padding: 20px;
                    box-shadow: 0 2px 8px rgba(0,0,0,0.1); border-left: 4px solid {model['scale_color']};">
            <div style="display: flex; align-items: center; margin-bottom: 10px;">
                <span style="background: {model['scale_color']}; color: white; padding: 3px 10px;
                            border-radius: 12px; font-size: 0.8em; margin-right: 8px;">
                    {model['scale_icon']} {model['scale']}
                </span>
            </div>
            <h3 style="margin: 5px 0; color: #333;">{model['name']}</h3>
            <p style="color: #888; font-size: 0.85em; margin: 2px 0 10px 0;">{model['name_en']}</p>
            <p style="color: #555; font-size: 0.9em; line-height: 1.5;">{model['description']}</p>
            <details style="margin-top: 10px;">
                <summary style="cursor: pointer; color: {model['scale_color']}; font-weight: bold;">
                    查看详细参数 ▼
                </summary>
                <div style="margin-top: 8px; background: #f9f9f9; padding: 12px; border-radius: 8px; font-size: 0.85em;">
                    <p><strong>算法：</strong><code>{model['algorithm']}</code></p>
                    <p><strong>参数：</strong></p>
                    <ul style="margin: 5px 0; padding-left: 20px;">{params_html}</ul>
                    <p><strong>输入：</strong>{model['input']}</p>
                    <p><strong>输出：</strong>{model['output']}</p>
                </div>
            </details>
        </div>
        """

    html = f"""
    <div style="display: grid; grid-template-columns: repeat(2, 1fr); gap: 16px; padding: 10px;">
        {cards_html}
    </div>
    """
    return html


def create_model_summary_table(catalog=None):
    """
    生成模型汇总表格DataFrame

    返回:
        pd.DataFrame 用于 gr.Dataframe 展示
    """
    if catalog is None:
        catalog = MODEL_CATALOG

    rows = []
    for model in catalog:
        rows.append({
            "模型名称": model["name"],
            "尺度": f"{model['scale_icon']} {model['scale']}",
            "算法": model["algorithm"],
            "输入": model["input"],
            "输出": model["output"]
        })

    return pd.DataFrame(rows)


def create_scale_distribution_data(catalog=None):
    """
    统计各尺度模型数量，用于饼图展示

    返回:
        dict: {尺度名称: 模型数量}
    """
    if catalog is None:
        catalog = MODEL_CATALOG

    distribution = {}
    for model in catalog:
        scale = f"{model['scale_icon']} {model['scale']}"
        distribution[scale] = distribution.get(scale, 0) + 1

    return distribution
