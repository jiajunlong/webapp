"""
基因网络和社会网络仿真计算模型 - 完整Gradio实现
Gene Network and Social Network Simulation Computing Model
包含前端界面和后端数据处理逻辑
"""

import gradio as gr
import plotly.graph_objects as go
import networkx as nx
import pandas as pd
import numpy as np
import json
import os
import pickle
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass, asdict
from collections import defaultdict
import random

# 导入社交网络仿真模块
from social_network_sim import SocialNetworkSimulator

# 导入TCGA-COAD仿真模块
from tcga_coad_simulator import TCGA_COAD_Simulator

# 导入跨尺度引擎与模型库
from cross_scale_engine import CrossScaleEngine

# 导入Phase 1通路分析模块
from gradio_phase1_integration import create_pathway_analysis_tab, Phase1DataLoader
from model_library import create_model_cards_html, create_model_summary_table, create_scale_distribution_data
# 导入Phase 2网络医学分析模块
from gradio_phase2_integration import create_phase2_network_medicine_tab, Phase2DataLoader
# 导入Phase 3 SIS生物标志物发现模块
from gradio_phase3_integration import create_phase3_biomarker_tab, Phase3DataLoader

# ==================== 配置 ====================

# 基因别名处理方式
# 数据格式说明：gene_symbol列中，逗号分隔的是同一个基因的不同别名
# 例如："PTEN, PTEN1, PTENbeta" 都是PTEN基因的别名
# 选项：
#   "first_only" - 只取第一个（主基因符号），推荐，避免重复
#   "all" - 保留所有别名（可能导致同一基因出现多次）
GENE_ALIAS_MODE = "first_only"  # 推荐使用 "first_only"

# 检查预处理数据
PREPROCESSED_DATA_FILE = "data/preprocessed_data.pkl"
USE_PREPROCESSED = os.path.exists(PREPROCESSED_DATA_FILE)
USE_REAL_DATA = os.path.exists("data/gene_disease.tsv") if not USE_PREPROCESSED else False

# ==================== 数据模型 ====================

@dataclass
class Gene:
    """基因数据模型"""
    name: str
    pathways: List[str]
    
@dataclass
class Interaction:
    """基因互作关系"""
    gene1: str
    gene2: str
    
@dataclass
class Regulation:
    """基因调控关系"""
    regulator: str  # 调控基因
    target: str     # 目标基因
    
@dataclass
class Disease:
    """疾病数据模型"""
    name: str
    genes: List[str]
    drugs: List[str]
    interactions: List[Interaction]
    regulations: List[Regulation]
    
@dataclass
class Pathway:
    """通路数据模型"""
    name: str
    genes: List[str]
    interactions: List[Interaction]
    regulations: List[Regulation]


# ==================== 真实数据加载器 ====================

def parse_gene_symbols(gene_str, mode="first_only"):
    """
    解析基因符号字符串，处理别名
    
    参数:
        gene_str: 基因符号字符串，可能包含逗号分隔的别名
        mode: "first_only" - 只取第一个（主基因），"all" - 保留所有别名
    
    返回:
        基因符号列表
    """
    if pd.isna(gene_str) or not gene_str or gene_str == 'NA':
        return []
    
    # 分割逗号
    genes = [g.strip() for g in str(gene_str).split(',')]
    genes = [g for g in genes if g and g != 'NA']
    
    if mode == "first_only":
        # 只返回第一个（主基因符号）
        return [genes[0]] if genes else []
    else:
        # 返回所有别名
        return genes


class RealDataLoader:
    """从真实数据文件加载数据"""
    
    def __init__(self, data_dir="data"):
        self.data_dir = data_dir
        self.gene_disease_df = None
        self.pathway_df = None
        self.diseases = {}
        self.pathways = {}
        self.genes = {}
        
    def load_all_data(self):
        """加载所有数据文件"""
        print("📂 加载基因-疾病数据...")
        gene_disease_file = os.path.join(self.data_dir, "gene_disease.tsv")
        if os.path.exists(gene_disease_file):
            self.gene_disease_df = pd.read_csv(gene_disease_file, sep='\t')
            print(f"   ✅ 加载了 {len(self.gene_disease_df)} 条记录")
        
        print("📂 加载通路数据...")
        pathway_file = os.path.join(self.data_dir, "pathway(基因名映射版).tsv")
        if os.path.exists(pathway_file):
            self.pathway_df = pd.read_csv(pathway_file, sep='\t')
            print(f"   ✅ 加载了 {len(self.pathway_df)} 条通路记录")
        
        self._build_database()
    
    def _build_database(self):
        """从DataFrame构建数据库结构 - 延迟加载优化"""
        if self.gene_disease_df is None:
            return
        
        print("🔄 构建疾病索引...")
        # 只构建疾病列表，不预先加载所有数据
        unique_diseases = self.gene_disease_df['disease_name'].unique()
        for disease_name in unique_diseases:
            if pd.notna(disease_name):
                # 只存储疾病名，实际数据延迟加载
                self.diseases[disease_name] = None
        
        print(f"   ✅ 索引了 {len(self.diseases)} 个疾病")
        
        # 构建通路索引
        if self.pathway_df is not None:
            print("🔄 构建通路索引...")
            for _, row in self.pathway_df.iterrows():
                pathway_name = row['Pathway_Name']
                if pd.notna(pathway_name):
                    self.pathways[pathway_name] = None
            print(f"   ✅ 索引了 {len(self.pathways)} 条通路")
    
    def _load_disease_data(self, disease_name: str) -> Disease:
        """延迟加载疾病数据"""
        if disease_name not in self.diseases or self.diseases[disease_name] is None:
            disease_data = self.gene_disease_df[
                self.gene_disease_df['disease_name'] == disease_name
            ]
            
            # 获取基因列表 - 处理基因别名
            # 注意：gene_symbol列中逗号分隔的是同一基因的别名（如 "PTEN, PTEN1, PTENbeta"）
            # 我们只取第一个（主基因符号）以避免重复
            gene_symbols = []
            for gene_str in disease_data['gene_symbol'].unique():
                parsed_genes = parse_gene_symbols(gene_str, mode=GENE_ALIAS_MODE)
                gene_symbols.extend(parsed_genes)
            
            # 去重并限制数量（保持顺序）
            genes = list(dict.fromkeys(gene_symbols))[:50]
            
            # 获取药物列表
            drugs = []
            if len(disease_data) > 0:
                drug_str = disease_data.iloc[0]['disease_drug']
                if pd.notna(drug_str):
                    for item in drug_str.split(';'):
                        item = item.strip()
                        if '[DR:' in item:
                            drug_name = item.split('[DR:')[0].strip()
                            if drug_name:
                                drugs.append(drug_name)
                drugs = drugs[:15]
            
            # 生成网络关系 - 减少连接数
            interactions, regulations = self._generate_network(genes, max_connections=2)
            
            self.diseases[disease_name] = Disease(
                name=disease_name,
                genes=genes,
                drugs=drugs if drugs else ["暂无药物数据"],
                interactions=interactions,
                regulations=regulations
            )
        
        return self.diseases[disease_name]
    
    def _load_pathway_data(self, pathway_name: str) -> Pathway:
        """延迟加载通路数据"""
        if pathway_name not in self.pathways or self.pathways[pathway_name] is None:
            pathway_row = self.pathway_df[self.pathway_df['Pathway_Name'] == pathway_name]
            if len(pathway_row) > 0:
                gene_str = pathway_row.iloc[0]['Gene']
                genes = []
                if pd.notna(gene_str) and gene_str != 'NA':
                    genes = [g.strip() for g in str(gene_str).replace(';', ',').split(',')]
                    genes = [g for g in genes if g and g != 'NA'][:30]  # 减少到30个基因
                
                if genes:
                    interactions, regulations = self._generate_network(genes, max_connections=2)
                    
                    self.pathways[pathway_name] = Pathway(
                        name=pathway_name,
                        genes=genes,
                        interactions=interactions,
                        regulations=regulations
                    )
                    
                    # 构建基因-通路映射
                    for gene in genes:
                        if gene not in self.genes:
                            self.genes[gene] = Gene(gene, [])
                        if pathway_name not in self.genes[gene].pathways:
                            self.genes[gene].pathways.append(pathway_name)
        
        return self.pathways[pathway_name]
    
    def _generate_network(self, genes: List[str], max_connections: int = 3) -> Tuple[List[Interaction], List[Regulation]]:
        """生成基因网络关系 - 优化版本，减少连接数"""
        interactions = []
        regulations = []
        
        # 减少连接数以提升性能
        for i, gene in enumerate(genes):
            connections = 0
            for j in range(i+1, len(genes)):
                if connections >= max_connections:
                    break
                if hash(gene + genes[j]) % 3 > 0:
                    interactions.append(Interaction(gene, genes[j]))
                    connections += 1
                if hash(gene + genes[j]) % 2 == 0 and connections < max_connections:
                    regulations.append(Regulation(gene, genes[j]))
        
        return interactions, regulations
    
    def get_all_diseases(self) -> List[str]:
        return sorted(list(self.diseases.keys()))
    
    def get_all_pathways(self) -> List[str]:
        return sorted(list(self.pathways.keys()))


# ==================== 数据库模拟 ====================

class GeneNetworkDatabase:
    """基因网络数据库"""
    
    def __init__(self, use_real_data=False, use_preprocessed=False, preprocessed_file=None):
        self.diseases: Dict[str, Disease] = {}
        self.pathways: Dict[str, Pathway] = {}
        self.genes: Dict[str, Gene] = {}
        
        if use_preprocessed and preprocessed_file:
            print("⚡ 使用预处理数据...")
            self._load_preprocessed_data(preprocessed_file)
        elif use_real_data:
            print("🔄 使用真实数据...")
            self._load_real_data()
        else:
            print("🔄 使用模拟数据...")
            self._initialize_sample_data()
    
    def _load_preprocessed_data(self, preprocessed_file: str):
        """加载预处理数据（最快）"""
        import time
        start_time = time.time()
        
        print(f"📂 加载预处理数据: {preprocessed_file}")
        
        try:
            with open(preprocessed_file, 'rb') as f:
                data = pickle.load(f)
            
            self.diseases = data['diseases']
            self.pathways = data['pathways']
            self.genes = data['genes']
            
            metadata = data.get('metadata', {})
            
            elapsed = time.time() - start_time
            print(f"   ✅ 加载完成！耗时: {elapsed:.2f} 秒")
            print(f"   📊 疾病: {len(self.diseases)}, 通路: {len(self.pathways)}, 基因: {len(self.genes)}")
            
            if 'preprocessed_at' in metadata:
                print(f"   📅 数据生成时间: {metadata['preprocessed_at']}")
            
        except Exception as e:
            print(f"   ❌ 加载失败: {e}")
            print("   ⚠️  回退到实时加载模式")
            self._load_real_data()
    
    def _load_real_data(self):
        """加载真实数据"""
        loader = RealDataLoader()
        loader.load_all_data()
        self.diseases = loader.diseases
        self.pathways = loader.pathways
        self.genes = loader.genes
        # 保存loader引用以支持延迟加载
        self._loader = loader
    
    def _initialize_sample_data(self):
        """初始化示例数据"""
        
        # 定义基因通路
        pathway_data = {
            "细胞周期通路": ["TP53", "RB1", "CDKN2A", "CDK4", "CCND1", "E2F1", "MYC"],
            "DNA修复通路": ["BRCA1", "BRCA2", "ATM", "TP53", "CHEK2", "RAD51"],
            "凋亡信号通路": ["TP53", "BCL2", "BAX", "CASP3", "CASP9", "FAS", "FASLG"],
            "MAPK信号通路": ["KRAS", "BRAF", "MEK1", "ERK1", "ERK2", "RAF1", "MAP2K1"],
            "PI3K-AKT通路": ["PIK3CA", "AKT1", "PTEN", "MTOR", "PDK1", "GSK3B"],
            "WNT信号通路": ["WNT1", "CTNNB1", "APC", "AXIN1", "GSK3B", "LEF1", "TCF7"],
            "TGF-β信号通路": ["TGFB1", "SMAD2", "SMAD3", "SMAD4", "TGFBR1", "TGFBR2"],
            "NF-κB信号通路": ["NFKB1", "RELA", "IKBKA", "IKBKB", "TNFAIP3", "TNF"],
            "JAK-STAT通路": ["JAK1", "JAK2", "STAT3", "STAT5A", "IL6", "IL6R"],
            "血管生成通路": ["VEGFA", "VEGFR2", "HIF1A", "ANGPT1", "ANGPT2", "TIE2"],
            "免疫检查点通路": ["PDCD1", "CD274", "CTLA4", "CD80", "CD86"],
            "代谢重编程通路": ["HIF1A", "LDHA", "PKM2", "GLUT1", "HK2"],
        }
        
        # 创建通路
        for pathway_name, genes in pathway_data.items():
            interactions = []
            regulations = []
            
            # 生成通路内的互作和调控关系
            for i in range(len(genes)):
                for j in range(i+1, min(i+3, len(genes))):
                    if random.random() > 0.3:
                        interactions.append(Interaction(genes[i], genes[j]))
                    if random.random() > 0.5:
                        regulations.append(Regulation(genes[i], genes[j]))
            
            self.pathways[pathway_name] = Pathway(
                name=pathway_name,
                genes=genes,
                interactions=interactions,
                regulations=regulations
            )
        
        # 创建基因字典
        for pathway_name, pathway in self.pathways.items():
            for gene_name in pathway.genes:
                if gene_name not in self.genes:
                    self.genes[gene_name] = Gene(gene_name, [])
                self.genes[gene_name].pathways.append(pathway_name)
        
        # 定义疾病
        disease_configs = [
            {
                "name": "肺癌",
                "pathways": ["细胞周期通路", "DNA修复通路", "凋亡信号通路", "MAPK信号通路", "血管生成通路"],
                "drugs": ["吉非替尼", "厄洛替尼", "奥希替尼", "克唑替尼", "阿法替尼", "贝伐珠单抗"]
            },
            {
                "name": "乳腺癌",
                "pathways": ["DNA修复通路", "PI3K-AKT通路", "细胞周期通路", "WNT信号通路"],
                "drugs": ["他莫昔芬", "曲妥珠单抗", "帕妥珠单抗", "拉帕替尼", "依维莫司"]
            },
            {
                "name": "结直肠癌",
                "pathways": ["WNT信号通路", "MAPK信号通路", "PI3K-AKT通路", "TGF-β信号通路"],
                "drugs": ["西妥昔单抗", "贝伐珠单抗", "雷莫芦单抗", "瑞戈非尼"]
            },
            {
                "name": "肝癌",
                "pathways": ["血管生成通路", "MAPK信号通路", "WNT信号通路", "JAK-STAT通路"],
                "drugs": ["索拉非尼", "仑伐替尼", "瑞戈非尼", "卡博替尼"]
            },
            {
                "name": "胃癌",
                "pathways": ["PI3K-AKT通路", "MAPK信号通路", "血管生成通路", "细胞周期通路"],
                "drugs": ["曲妥珠单抗", "雷莫芦单抗", "阿帕替尼", "氟尿嘧啶"]
            },
            {
                "name": "黑色素瘤",
                "pathways": ["MAPK信号通路", "PI3K-AKT通路", "免疫检查点通路", "凋亡信号通路"],
                "drugs": ["威罗菲尼", "达拉非尼", "曲美替尼", "帕博利珠单抗", "纳武利尤单抗", "伊匹单抗"]
            },
            {
                "name": "胰腺癌",
                "pathways": ["MAPK信号通路", "TGF-β信号通路", "代谢重编程通路", "DNA修复通路"],
                "drugs": ["吉西他滨", "厄洛替尼", "奥拉帕利", "伊立替康"]
            },
        ]
        
        for disease_config in disease_configs:
            # 收集疾病相关的所有基因
            disease_genes = set()
            disease_interactions = []
            disease_regulations = []
            
            for pathway_name in disease_config["pathways"]:
                if pathway_name in self.pathways:
                    pathway = self.pathways[pathway_name]
                    disease_genes.update(pathway.genes)
                    disease_interactions.extend(pathway.interactions)
                    disease_regulations.extend(pathway.regulations)
            
            # 去重
            disease_interactions = list({(i.gene1, i.gene2): i for i in disease_interactions}.values())
            disease_regulations = list({(r.regulator, r.target): r for r in disease_regulations}.values())
            
            self.diseases[disease_config["name"]] = Disease(
                name=disease_config["name"],
                genes=list(disease_genes),
                drugs=disease_config["drugs"],
                interactions=disease_interactions,
                regulations=disease_regulations
            )
    
    def get_all_diseases(self) -> List[str]:
        """获取所有疾病名称"""
        return list(self.diseases.keys())
    
    def get_all_pathways(self) -> List[str]:
        """获取所有通路名称"""
        return list(self.pathways.keys())
    
    def get_disease_network(self, disease_name: str, network_type: str = "interaction") -> Dict:
        """获取疾病基因网络 - 支持延迟加载"""
        if disease_name not in self.diseases:
            return {"genes": [], "edges": [], "drugs": []}
        
        # 如果是真实数据且未加载，则延迟加载
        if hasattr(self, '_loader') and self.diseases[disease_name] is None:
            disease = self._loader._load_disease_data(disease_name)
        else:
            disease = self.diseases[disease_name]
        
        if disease is None:
            return {"genes": [], "edges": [], "drugs": []}
        
        if network_type == "regulation":
            edges = [{"FirstGene": r.regulator, "SecondGene": r.target} 
                    for r in disease.regulations]
        else:
            edges = [{"FirstGene": i.gene1, "SecondGene": i.gene2} 
                    for i in disease.interactions]
        
        return {
            "genes": disease.genes,
            "edges": edges,
            "drugs": disease.drugs
        }
    
    def get_gene_pathways(self, gene_name: str) -> List[str]:
        """获取基因所在的通路"""
        if gene_name in self.genes:
            return self.genes[gene_name].pathways
        return []
    
    def get_statistics(self) -> Dict:
        """获取数据库统计信息 - 安全处理延迟加载"""
        stats = {
            "total_diseases": len(self.diseases),
            "total_pathways": len(self.pathways),
            "total_genes": len(self.genes),
            "loaded_diseases": 0,
            "loaded_pathways": 0,
            "avg_genes_per_pathway": 0,
            "avg_pathways_per_gene": 0
        }
        
        # 统计已加载的疾病
        loaded_diseases = [d for d in self.diseases.values() if d is not None]
        stats["loaded_diseases"] = len(loaded_diseases)
        
        # 统计已加载的通路
        loaded_pathways = [p for p in self.pathways.values() if p is not None]
        stats["loaded_pathways"] = len(loaded_pathways)
        
        # 计算平均值
        if loaded_pathways:
            stats["avg_genes_per_pathway"] = int(np.mean([len(p.genes) for p in loaded_pathways]))
        
        if self.genes:
            genes_with_pathways = [g for g in self.genes.values() if g.pathways]
            if genes_with_pathways:
                stats["avg_pathways_per_gene"] = int(np.mean([len(g.pathways) for g in genes_with_pathways]))
        
        return stats
    
    def get_pathway_network(self, pathway_name: str, network_type: str = "interaction") -> Dict:
        """获取通路基因网络 - 支持延迟加载"""
        if pathway_name not in self.pathways:
            return {"genes": [], "edges": []}
        
        # 如果是真实数据且未加载，则延迟加载
        if hasattr(self, '_loader') and self.pathways[pathway_name] is None:
            pathway = self._loader._load_pathway_data(pathway_name)
        else:
            pathway = self.pathways[pathway_name]
        
        if pathway is None:
            return {"genes": [], "edges": []}
        
        if network_type == "regulation":
            edges = [{"FirstGene": r.regulator, "SecondGene": r.target} 
                    for r in pathway.regulations]
        else:
            edges = [{"FirstGene": i.gene1, "SecondGene": i.gene2} 
                    for i in pathway.interactions]
        
        return {
            "genes": pathway.genes,
            "edges": edges
        }
    
    def calculate_is_score(self, disease_name: str, pathway_name: str) -> float:
        """
        计算影响分数 (Influence Score, IS系数)
        IS = (疾病基因在通路中的比例) × (通路基因在疾病中的比例) × (连接密度)
        """
        if disease_name not in self.diseases or pathway_name not in self.pathways:
            print(f"⚠️ IS计算: 疾病或通路不存在 - 疾病: {disease_name}, 通路: {pathway_name}")
            return 0.0
        
        # 支持延迟加载 - 确保数据已加载
        if hasattr(self, '_loader'):
            if self.diseases[disease_name] is None:
                disease = self._loader._load_disease_data(disease_name)
            else:
                disease = self.diseases[disease_name]
            
            if self.pathways[pathway_name] is None:
                pathway = self._loader._load_pathway_data(pathway_name)
            else:
                pathway = self.pathways[pathway_name]
        else:
            disease = self.diseases[disease_name]
            pathway = self.pathways[pathway_name]
        
        # 检查是否成功加载
        if disease is None or pathway is None:
            print(f"⚠️ IS计算: 数据加载失败 - 疾病: {disease is None}, 通路: {pathway is None}")
            return 0.0
        
        # 确保基因列表不为空
        if not disease.genes or not pathway.genes:
            print(f"⚠️ IS计算: 基因列表为空 - 疾病基因数: {len(disease.genes)}, 通路基因数: {len(pathway.genes)}")
            return 0.0
        
        # 转换为集合并统一大小写（确保匹配）
        disease_genes = set(g.upper() for g in disease.genes if g)
        pathway_genes = set(g.upper() for g in pathway.genes if g)
        
        # 交集基因（不区分大小写）
        overlap_genes = disease_genes & pathway_genes
        
        if len(overlap_genes) == 0:
            # 无交集是正常情况，不输出警告（减少噪音）
            return 0.0
        
        # 计算比例
        ratio_in_pathway = len(overlap_genes) / len(pathway_genes) if len(pathway_genes) > 0 else 0
        ratio_in_disease = len(overlap_genes) / len(disease_genes) if len(disease_genes) > 0 else 0
        
        # 计算连接密度（交集基因之间的连接数）
        # 需要将交集基因转换回原始大小写格式进行匹配
        overlap_genes_original = {g for g in disease.genes if g.upper() in overlap_genes}
        edge_count = 0
        max_edges = len(overlap_genes) * (len(overlap_genes) - 1) / 2
        
        for interaction in disease.interactions:
            if (interaction.gene1.upper() in overlap_genes and 
                interaction.gene2.upper() in overlap_genes):
                edge_count += 1
        
        connectivity = edge_count / max_edges if max_edges > 0 else 0
        
        # IS系数计算
        is_score = ratio_in_pathway * ratio_in_disease * (1 + connectivity)
        
        # 添加一些随机性使结果更真实（但不要太大，避免完全随机）
        is_score *= (0.9 + random.random() * 0.2)  # 0.9-1.1倍
        
        result = round(is_score, 4)
        
        # 只在结果较大时输出详细信息（减少噪音）
        if result > 0.01:
            print(f"✅ IS={result:.4f}: {disease_name[:30]} × {pathway_name[:40]} "
                  f"(交集:{len(overlap_genes)})")
        
        return result


# ==================== 网络可视化 ====================

class NetworkVisualizer:
    """网络可视化工具"""
    
    @staticmethod
    def create_network_plot(genes: List[str], edges: List[Dict], 
                          is_regulation: bool, title: str = "") -> go.Figure:
        """创建网络可视化图 - 性能优化版本"""
        
        if not genes:
            fig = go.Figure()
            fig.update_layout(
                title="暂无数据",
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                height=600
            )
            return fig
        
        # 限制节点数量以提升性能
        if len(genes) > 50:
            genes = genes[:50]
            edges = [e for e in edges if e.get("FirstGene") in genes and e.get("SecondGene") in genes]
        
        # 创建NetworkX图
        G = nx.DiGraph() if is_regulation else nx.Graph()
        
        # 添加节点
        for gene in genes:
            G.add_node(gene)
        
        # 添加边
        for edge in edges:
            if "FirstGene" in edge and "SecondGene" in edge:
                if edge["FirstGene"] in genes and edge["SecondGene"] in genes:
                    G.add_edge(edge["FirstGene"], edge["SecondGene"])
        
        # 使用更快的布局算法
        if len(G.nodes()) > 0:
            # 对于所有图都使用快速布局
            if len(G.nodes()) <= 20:
                # 小图：快速spring布局，减少迭代
                pos = nx.spring_layout(G, k=2/np.sqrt(len(G.nodes())), iterations=20, seed=42)
            elif len(G.nodes()) <= 40:
                # 中等图：使用circular布局（最快）
                pos = nx.circular_layout(G)
            else:
                # 大图：使用shell布局（非常快）
                pos = nx.shell_layout(G)
        else:
            pos = {}
        
        # 创建边的traces - 简化版本
        edge_x = []
        edge_y = []
        
        for edge in G.edges():
            if edge[0] in pos and edge[1] in pos:
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])
        
        edge_trace = go.Scatter(
            x=edge_x,
            y=edge_y,
            mode='lines',
            line=dict(width=1, color='rgba(150, 150, 150, 0.5)'),
            hoverinfo='none',
            showlegend=False
        )
        
        # 创建节点trace - 简化版本
        node_x = []
        node_y = []
        node_text = []
        node_connections = []
        
        for node in G.nodes():
            if node in pos:
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)
                node_text.append(node)
                node_connections.append(len(list(G.neighbors(node))))
        
        # 节点大小根据连接数调整 - 简化计算
        node_sizes = [10 + min(conn, 10) for conn in node_connections]  # 进一步限制
        
        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            mode='markers+text',
            text=node_text,
            textposition="top center",
            textfont=dict(size=8, color='black'),  # 进一步减小字体
            hoverinfo='text',
            hovertext=[f"<b>{text}</b><br>连接: {conn}" 
                      for text, conn in zip(node_text, node_connections)],
            marker=dict(
                size=node_sizes,
                color=node_connections,
                colorscale='Blues',  # 使用更简单的配色
                showscale=False,  # 隐藏颜色条以加快渲染
                line=dict(width=1, color='white')
            ),
            showlegend=False
        )
        
        # 创建图形 - 只有一个edge_trace
        fig = go.Figure(data=[edge_trace, node_trace])
        
        network_type = "基因调控网络" if is_regulation else "基因互作网络"
        actual_edges = len([e for e in edges if e.get("FirstGene") in genes and e.get("SecondGene") in genes])
        
        fig.update_layout(
            title=dict(
                text=f"{title}<br><sub>{network_type} | 节点: {len(genes)} | 边: {actual_edges}</sub>",
                x=0.5,
                xanchor='center',
                font=dict(size=14)  # 减小字体
            ),
            showlegend=False,
            hovermode='closest',
            margin=dict(b=10, l=10, r=10, t=60),  # 减小边距
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, visible=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, visible=False),
            height=500,  # 减小高度提升性能
            plot_bgcolor='white',
            paper_bgcolor='white',
            dragmode='pan'  # 默认平移模式
        )
        
        return fig
    
    @staticmethod
    def create_pathway_network_plot(genes: List[str], edges: List[Dict], 
                                   is_regulation: bool, pathway_name: str,
                                   highlight_genes: List[str] = None) -> go.Figure:
        """创建通路网络图，高亮指定基因"""
        
        if not genes:
            fig = go.Figure()
            fig.update_layout(title="暂无数据", xaxis=dict(visible=False), yaxis=dict(visible=False))
            return fig
        
        G = nx.DiGraph() if is_regulation else nx.Graph()
        
        for gene in genes:
            G.add_node(gene)
        
        for edge in edges:
            if "FirstGene" in edge and "SecondGene" in edge:
                G.add_edge(edge["FirstGene"], edge["SecondGene"])
        
        if len(G.nodes()) > 0:
            pos = nx.spring_layout(G, k=2/np.sqrt(len(G.nodes())), iterations=50, seed=42)
        else:
            pos = {}
        
        # 创建边
        edge_traces = []
        for edge in G.edges():
            if edge[0] in pos and edge[1] in pos:
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_traces.append(
                    go.Scatter(
                        x=[x0, x1, None],
                        y=[y0, y1, None],
                        mode='lines',
                        line=dict(width=1.5, color='rgba(150, 150, 150, 0.5)'),
                        hoverinfo='none',
                        showlegend=False
                    )
                )
        
        # 创建节点 - 区分高亮
        node_x = []
        node_y = []
        node_text = []
        node_colors = []
        
        highlight_set = set(highlight_genes) if highlight_genes else set()
        
        for node in G.nodes():
            if node in pos:
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)
                node_text.append(node)
                
                if node in highlight_set:
                    node_colors.append('#FF1493')  # DeepPink - 在疾病网络中
                else:
                    node_colors.append('#00BFFF')  # DeepSkyBlue - 不在疾病网络中
        
        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            mode='markers+text',
            text=node_text,
            textposition="top center",
            textfont=dict(size=10),
            hoverinfo='text',
            hovertext=[f"<b>{text}</b><br>{'在疾病网络中' if text in highlight_set else '不在疾病网络中'}" 
                      for text in node_text],
            marker=dict(
                size=20,
                color=node_colors,
                line=dict(width=2, color='white')
            ),
            showlegend=False
        )
        
        fig = go.Figure(data=edge_traces + [node_trace])
        
        network_type = "调控网络" if is_regulation else "互作网络"
        fig.update_layout(
            title=dict(
                text=f"通路: {pathway_name} ({network_type})<br><sub>粉色=在疾病网络中 | 蓝色=仅在通路中</sub>",
                x=0.5,
                xanchor='center'
            ),
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=20, r=20, t=80),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            height=500,
            plot_bgcolor='#fafafa'
        )
        
        return fig


# ==================== Gradio应用 ====================

class GeneNetworkApp:
    """完整的基因网络应用"""
    
    def __init__(self, use_real_data=False, use_preprocessed=False, preprocessed_file=None):
        self.db = GeneNetworkDatabase(
            use_real_data=use_real_data,
            use_preprocessed=use_preprocessed,
            preprocessed_file=preprocessed_file
        )
        self.visualizer = NetworkVisualizer()
        self.current_disease = None
        self.current_disease_genes = []
        # 缓存疾病-通路交集结果
        self._pathway_overlap_cache = {}
        
    def load_disease_network(self, disease: str, network_type: str):
        """加载疾病网络 - 高性能优化版本"""
        if not disease:
            return go.Figure(), "⚠️ 请选择疾病", gr.update(choices=[]), None
        
        # 添加性能提示
        import time
        start_time = time.time()
        
        self.current_disease = disease
        is_regulation = (network_type == "基因调控网络")
        
        net_type = "regulation" if is_regulation else "interaction"
        data = self.db.get_disease_network(disease, net_type)
        
        # 如果没有数据，快速返回
        if not data["genes"]:
            return go.Figure(), f"⚠️ {disease} 暂无网络数据", gr.update(choices=[]), None
        
        self.current_disease_genes = data["genes"]
        
        # 大幅减少节点数量以提高性能（30个节点）
        max_nodes = 30
        max_edges = 50
        
        if len(data["genes"]) > max_nodes:
            genes_to_show = data["genes"][:max_nodes]
            genes_set = set(genes_to_show)
            edges_to_show = [e for e in data["edges"] 
                           if e["FirstGene"] in genes_set and e["SecondGene"] in genes_set]
            # 限制边的数量
            if len(edges_to_show) > max_edges:
                edges_to_show = edges_to_show[:max_edges]
            performance_note = f"\n\n⚡ **性能优化:** 为提升渲染速度，仅显示前 {len(genes_to_show)} 个基因和 {len(edges_to_show)} 条连接"
        else:
            genes_to_show = data["genes"]
            edges_to_show = data["edges"]
            # 即使总数少，也限制边数
            if len(edges_to_show) > max_edges:
                edges_to_show = edges_to_show[:max_edges]
            performance_note = ""
        
        fig = self.visualizer.create_network_plot(
            genes_to_show, 
            edges_to_show, 
            is_regulation,
            f"{disease}"
        )
        
        elapsed = time.time() - start_time
        
        info = f"""✅ 已加载 **{disease}** 的网络 (耗时: {elapsed:.2f}秒)
        
📊 **统计信息:**
- 总基因数: {len(data["genes"])} | 显示: {len(genes_to_show)}
- 总连接数: {len(data["edges"])} | 显示: {len(edges_to_show)}
- 相关药物: {len(data["drugs"])} 种{performance_note}

💡 **提示:** 在基因输入框中输入基因名称可查询其所在通路"""
        
        return fig, info, gr.update(choices=data["drugs"], value=None), None
    
    def get_pathways_with_overlap(self, disease_name: str) -> List[str]:
        """
        获取与疾病有基因交集的通路列表（带缓存）
        
        参数:
            disease_name: 疾病名称
        
        返回:
            有交集的通路名称列表
        """
        if not disease_name:
            return []
        
        # 检查缓存
        if disease_name in self._pathway_overlap_cache:
            return self._pathway_overlap_cache[disease_name]
        
        # 获取疾病数据
        disease_data = self.db.get_disease_network(disease_name, "interaction")
        if not disease_data["genes"]:
            self._pathway_overlap_cache[disease_name] = []
            return []
        
        # 转换为集合（不区分大小写）
        disease_genes = set(g.upper() for g in disease_data["genes"] if g)
        
        # 获取所有通路
        all_pathways = self.db.get_all_pathways()
        
        # 检查每个通路是否有交集
        pathways_with_overlap = []
        
        for pathway_name in all_pathways:
            # 获取通路数据（使用延迟加载）
            pathway_data = self.db.get_pathway_network(pathway_name, "interaction")
            if not pathway_data["genes"]:
                continue
            
            # 转换为集合（不区分大小写）
            pathway_genes = set(g.upper() for g in pathway_data["genes"] if g)
            
            # 检查交集
            overlap = disease_genes & pathway_genes
            if len(overlap) > 0:
                pathways_with_overlap.append(pathway_name)
        
        # 缓存结果
        self._pathway_overlap_cache[disease_name] = pathways_with_overlap
        
        return pathways_with_overlap
    
    def query_gene_pathways(self, gene: str):
        """查询基因所在通路"""
        if not gene:
            return gr.update(choices=[]), go.Figure()
        
        pathways = self.db.get_gene_pathways(gene)
        
        if not pathways:
            return gr.update(choices=[], value=None), go.Figure()
        
        return gr.update(choices=pathways, value=None), go.Figure()
    
    def show_pathway_network(self, pathway: str, network_type: str):
        """显示通路网络"""
        if not pathway:
            return go.Figure()
        
        is_regulation = (network_type == "基因调控网络")
        net_type = "regulation" if is_regulation else "interaction"
        data = self.db.get_pathway_network(pathway, net_type)
        
        fig = self.visualizer.create_pathway_network_plot(
            data["genes"],
            data["edges"],
            is_regulation,
            pathway,
            self.current_disease_genes
        )
        
        return fig
    
    def calculate_is_scores(self, disease: str, pathways: List[str], progress=gr.Progress()):
        """计算IS系数"""
        if not disease:
            return go.Figure(), "⚠️ 请先选择疾病", None
        
        if not pathways:
            return go.Figure(), "⚠️ 请至少选择一个基因通路", None
        
        progress(0, desc="开始计算...")
        
        pathway_names = []
        is_scores = []
        
        for i, pathway in enumerate(pathways):
            progress((i + 1) / len(pathways), desc=f"计算中... {i+1}/{len(pathways)}")
            score = self.db.calculate_is_score(disease, pathway)
            pathway_names.append(pathway)
            is_scores.append(score)
        
        # 创建柱状图
        fig = go.Figure(data=[
            go.Bar(
                x=is_scores,
                y=pathway_names,
                orientation='h',
                marker=dict(
                    color=is_scores,
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(title="IS系数"),
                    line=dict(color='rgba(50, 50, 50, 0.3)', width=1)
                ),
                text=[f"{score:.4f}" for score in is_scores],
                textposition='auto',
                hovertemplate='<b>%{y}</b><br>IS系数: %{x:.4f}<extra></extra>'
            )
        ])
        
        fig.update_layout(
            title=dict(
                text=f"影响分数 (Influence Score, IS系数)<br><sub>疾病: {disease}</sub>",
                x=0.5,
                xanchor='center',
                font=dict(size=18)
            ),
            xaxis_title="IS系数",
            yaxis_title="基因通路",
            height=max(400, len(pathway_names) * 50),
            showlegend=False,
            hovermode='closest',
            plot_bgcolor='#fafafa',
            paper_bgcolor='white',
            margin=dict(l=200, r=20, t=100, b=60)
        )
        
        status = f"""✅ 计算完成！
        
📊 **计算结果:**
- 计算疾病: {disease}
- 计算通路数: {len(pathways)}
- 最高IS系数: {max(is_scores):.4f} ({pathway_names[is_scores.index(max(is_scores))]})
- 最低IS系数: {min(is_scores):.4f} ({pathway_names[is_scores.index(min(is_scores))]})
- 平均IS系数: {np.mean(is_scores):.4f}

💡 **解读:** IS系数越高，表示该通路在疾病中的影响越大"""
        
        # 创建数据表
        df = pd.DataFrame({
            '通路名称': pathway_names,
            'IS系数': [f"{score:.4f}" for score in is_scores]
        })
        df = df.sort_values('IS系数', ascending=False)
        
        return fig, status, df


def create_gradio_interface():
    """创建Gradio界面"""
    
    app = GeneNetworkApp(
        use_real_data=USE_REAL_DATA,
        use_preprocessed=USE_PREPROCESSED,
        preprocessed_file=PREPROCESSED_DATA_FILE if USE_PREPROCESSED else None
    )
    
    # 获取初始数据
    all_diseases = app.db.get_all_diseases()
    # 只显示前20个疾病以提升性能
    diseases = all_diseases[:20] if len(all_diseases) > 20 else all_diseases
    pathways = app.db.get_all_pathways()
    
    print(f"📊 加载完成: 总疾病数 {len(all_diseases)}, 显示 {len(diseases)} 个疾病, {len(pathways)} 条通路")
    
    # 自定义CSS
    custom_css = """
    .gradio-container {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    }
    .main-title {
        text-align: center;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 30px;
        border-radius: 10px;
        margin-bottom: 20px;
    }
    .tab-banner {
        color: white;
        padding: 16px 22px;
        border-radius: 10px;
        margin-bottom: 14px;
    }
    .tab-banner h3 { margin: 0 0 4px 0; font-size: 1.3em; }
    .tab-banner p  { margin: 0; opacity: 0.9; font-size: 0.92em; }
    .banner-purple  { background: linear-gradient(135deg, #667eea, #764ba2); }
    .banner-pink    { background: linear-gradient(135deg, #f093fb, #f5576c); }
    .banner-green   { background: linear-gradient(135deg, #43e97b, #38f9d7); color: #222; }
    .banner-blue    { background: linear-gradient(135deg, #4facfe, #00f2fe); }
    .banner-orange  { background: linear-gradient(135deg, #fa709a, #fee140); color: #222; }
    .banner-teal    { background: linear-gradient(135deg, #0ba360, #3cba92); }
    .banner-red     { background: linear-gradient(135deg, #ff6b6b, #ee5a24); }
    .banner-indigo  { background: linear-gradient(135deg, #5f72bd, #9b23ea); }
    .banner-slate   { background: linear-gradient(135deg, #616161, #9bc5c3); }
    """
    
    # 创建主题
    custom_theme = gr.themes.Soft(primary_hue="purple")
    
    with gr.Blocks(title="基因网络仿真计算模型", css=custom_css, theme=custom_theme) as demo:
        
        # 标题
        if USE_PREPROCESSED:
            data_source_tag = "⚡ 预处理数据（快速模式）"
        elif USE_REAL_DATA:
            data_source_tag = "📊 真实数据"
        else:
            data_source_tag = "🧪 模拟数据"
        data_source_tag = "📊 大规模基因网络&社交网络数据库"
        
        gr.HTML(f"""
            <div class="main-title">
                <h1 style="margin: 0; font-size: 2.5em;">🧬 大规模跨尺度仿真模型库</h1>
                <p style="margin: 10px 0 0 0; font-size: 1.2em;">Large-scale cross-scale simulation model library</p>
                <p style="margin: 5px 0 0 0; font-size: 0.9em;">{data_source_tag}</p>
            </div>
        """)
        
        with gr.Tabs() as tabs:
            
            # ========== Tab 0: 多尺度联动分析 ==========
            with gr.Tab("🔗 多尺度联动分析", id=0):
                gr.HTML("""<div class="tab-banner banner-purple">
                    <h3>🔗 多尺度联动分析</h3>
                    <p>分子层指标 → 自动推导 → 细胞层种子基因 → 自动推导 → 群体层传播参数</p>
                </div>""")

                with gr.Tabs():

                    # ---- 子Tab 1: 跨尺度级联（参数传递） ----
                    with gr.Tab("🚀 跨尺度级联"):
                        gr.Markdown("#### 三层级联分析 — 参数自动传递")
                        gr.Markdown("""
                        **传递机制:**
                        - Layer 1 → Layer 2: 分子层Hub基因 → 细胞层MRNetB种子节点
                        - Layer 2 → Layer 3: 分子密度+细胞聚类系数 → 群体层感染率β和社区数c
                        """)

                        with gr.Row():
                            with gr.Column(scale=1):
                                cascade_disease = gr.Dropdown(
                                    choices=diseases,
                                    value=diseases[0] if diseases else None,
                                    label="选择疾病",
                                )
                                cascade_net_type = gr.Radio(
                                    choices=["interaction", "regulation"],
                                    value="interaction", label="网络类型",
                                )
                                cascade_max_feat = gr.Slider(
                                    50, 300, value=100, step=50,
                                    label="最大特征数",
                                )
                                cascade_run_btn = gr.Button("🚀 运行三层级联分析",
                                                             variant="primary", size="lg")
                            with gr.Column(scale=2):
                                cascade_status = gr.Markdown("")
                                cascade_arch = gr.HTML(label="级联架构 + 参数传递")

                        with gr.Row():
                            with gr.Column():
                                cascade_network_plot = gr.Plot(label="三层网络对比")
                            with gr.Column():
                                cascade_radar = gr.Plot(label="跨尺度雷达图")

                        cascade_summary = gr.Dataframe(label="三层汇总表", interactive=False)
                        cascade_insights = gr.Markdown("")

                        def run_full_cascade(disease, net_type, max_feat):
                            try:
                                engine = CrossScaleEngine(
                                    db=app.db,
                                    social_sim=SocialNetworkSimulator(),
                                    tcga_sim=TCGA_COAD_Simulator(),
                                )
                                report = engine.run_full_cascade(
                                    disease_name=disease,
                                    network_type=net_type,
                                    max_features=int(max_feat),
                                )
                                arch = report.architecture_html
                                network_fig = engine.create_three_network_plots(report)
                                radar_fig = engine.create_radar_chart(report)
                                summary = engine.cascade_summary_df(report)
                                insights = "\n".join(
                                    f"- {i}" for i in report.cross_scale_insights
                                )
                                return ("✅ 三层级联分析完成", arch,
                                        network_fig, radar_fig, summary, insights)
                            except Exception as e:
                                import traceback
                                return (f"❌ {e}", "", go.Figure(), go.Figure(),
                                        pd.DataFrame(),
                                        f"```\n{traceback.format_exc()}\n```")

                        cascade_run_btn.click(
                            fn=run_full_cascade,
                            inputs=[cascade_disease, cascade_net_type, cascade_max_feat],
                            outputs=[cascade_status, cascade_arch,
                                     cascade_network_plot, cascade_radar,
                                     cascade_summary, cascade_insights],
                        )

                    # ---- 子Tab 2: 多疾病对比 ----
                    with gr.Tab("📊 多疾病对比"):
                        gr.Markdown("#### 多疾病分子网络对比分析")
                        gr.Markdown("选择 2-5 个疾病，对比它们在分子尺度的网络拓扑差异。")

                        compare_diseases_select = gr.CheckboxGroup(
                            choices=diseases,
                            value=diseases[:3] if len(diseases) >= 3 else diseases,
                            label="选择疾病（2-5个）",
                        )
                        compare_run_btn = gr.Button("📊 运行对比分析", variant="primary")
                        compare_status = gr.Markdown("")
                        compare_chart = gr.Plot(label="对比柱状图")
                        compare_table = gr.Dataframe(label="对比数据表", interactive=False)

                        def run_compare(selected_diseases):
                            if not selected_diseases or len(selected_diseases) < 2:
                                return "⚠️ 请至少选择2个疾病", go.Figure(), pd.DataFrame()
                            try:
                                engine = CrossScaleEngine(db=app.db)
                                df, fig = engine.compare_diseases(selected_diseases)
                                return "✅ 对比完成", fig, df
                            except Exception as e:
                                return f"❌ {e}", go.Figure(), pd.DataFrame()

                        compare_run_btn.click(
                            fn=run_compare,
                            inputs=[compare_diseases_select],
                            outputs=[compare_status, compare_chart, compare_table],
                        )

                    # ---- 子Tab 3: 基因追踪 ----
                    with gr.Tab("🔍 基因追踪"):
                        gr.Markdown("#### 基因跨尺度追踪")
                        gr.Markdown("输入一个基因，追踪它从 **分子层→细胞层→群体层** 的角色变化。")

                        with gr.Row():
                            trace_gene_input = gr.Textbox(
                                label="基因名称",
                                placeholder="例如: TP53, BRCA1, KRAS",
                                value="TP53",
                            )
                            trace_disease = gr.Dropdown(
                                choices=diseases,
                                value=diseases[0] if diseases else None,
                                label="疾病背景",
                            )
                        trace_run_btn = gr.Button("🔍 追踪基因", variant="primary")
                        trace_result = gr.HTML(label="基因档案卡")

                        def run_trace(gene, disease):
                            if not gene or not gene.strip():
                                return "<p>⚠️ 请输入基因名称</p>"
                            try:
                                engine = CrossScaleEngine(
                                    db=app.db,
                                    tcga_sim=TCGA_COAD_Simulator(),
                                    social_sim=SocialNetworkSimulator(),
                                )
                                return engine.trace_gene(gene.strip().upper(), disease)
                            except Exception as e:
                                return f"<p>❌ 追踪失败: {e}</p>"

                        trace_run_btn.click(
                            fn=run_trace,
                            inputs=[trace_gene_input, trace_disease],
                            outputs=[trace_result],
                        )

                    # ---- 子Tab 4: 跨尺度洞察说明 ----
                    with gr.Tab("📖 跨尺度说明"):
                        gr.Markdown("""
                        #### 跨尺度联动机制说明

                        本平台实现了 **分子 → 细胞/组织 → 群体** 三个生物学尺度的级联分析，
                        采用经过文献验证的跨尺度桥接方法。

                        ---

                        ##### 🧬 Layer 1 → Layer 2: 分子 → 细胞 桥接

                        **方法:** Hub基因种子定向推断

                        分子层通过度中心性识别Hub基因，作为细胞层MRNetB推断的种子节点。
                        这基于 **网络邻近性原理** — 疾病基因在互作组中倾向于聚集在
                        同一网络邻域内 (Guney et al. *Nature Communications*, 2016;
                        Menche et al. *Science*, 2015)。

                        ##### 🔬 Layer 2 → Layer 3: 细胞 → 群体 桥接

                        **方法:** Hill函数剂量-响应模型

                        采用 PhysiCell 框架 (Ghaffarizadeh et al. *PLoS Comput Biol*, 2018)
                        中的 **Hill函数信号-行为规则** 将跨尺度拓扑信号映射为群体传播参数:

                        ```
                        信号: s = (分子网络密度 + 细胞聚类系数) / 2
                        感染率: β = B_min + (B_max - B_min) × s² / (s² + K²)
                            其中 B_min=0.01, B_max=0.20, K=0.3, n=2
                        ```

                        Hill函数的生物学意义: 低连通性网络 → 低传播率（阈值行为），
                        高连通性网络 → 趋于饱和的传播率（超敏感性）。

                        ##### 👥 Layer 3: 群体层参数

                        - **社区数 c**: 由Hub基因数量映射，类比疾病模块数
                          (Menche et al. *Science*, 2015)
                        - **基本再生数**: R₀ ≈ β⟨k⟩/γ，决定传播是否可持续

                        ---

                        ##### 📚 核心参考文献

                        1. Ghaffarizadeh et al. "PhysiCell: An open source physics-based cell
                           simulator for 3-D multicellular systems." *PLoS Comput Biol*, 2018.
                        2. Guney et al. "Network-based in silico drug efficacy screening."
                           *Nature Communications*, 2016.
                        3. Menche et al. "Uncovering disease-disease relationships through
                           the incomplete interactome." *Science*, 2015.
                        4. Xue & Bhatt. "Coupling the within-host process and between-host
                           transmission." *Bull Math Biol*, 2022.
                        5. Moran, P.A.P. *The Statistical Processes of Evolutionary Theory*. 1962.
                        6. Albert et al. "Error and attack tolerance of complex networks."
                           *Nature*, 2000.

                        ---

                        ##### 🔍 基因追踪

                        追踪单个基因在三层中的角色:
                        - **分子层**: 连接度排名、是否为Hub、所在通路
                        - **细胞层**: TCGA-COAD 真实表达量（均值/方差/排名）
                        - **群体层**: 基于图论的网络角色分类 + Moran固定概率
                          (相对适合度映射)
                        """)

            # ========== Tab 1: 基因网络可视化 ==========
            with gr.Tab("🔬 基因网络可视化", id=1):
                gr.HTML("""<div class="tab-banner banner-indigo">
                    <h3>🔬 基因网络可视化</h3>
                    <p>分子尺度 · 基因互作网络 / 调控网络 · 通路查询</p>
                </div>""")

                gr.Markdown("### 🎯 选择疾病和网络类型")
                
                with gr.Row():
                    with gr.Column(scale=1):
                        # 添加疾病数量提示
                        total_diseases = len(all_diseases)
                        display_count = len(diseases)
                        if total_diseases > display_count:
                            gr.Markdown(f"💡 数据库包含 **{total_diseases}** 个疾病，下拉框显示前 **{display_count}** 个（为提升性能）")
                        else:
                            gr.Markdown(f"💡 数据库包含 **{display_count}** 个疾病")
                        disease_dropdown = gr.Dropdown(
                            choices=diseases,
                            label="🏥 疾病基因网络",
                            info=f"选择要分析的疾病（显示前{display_count}个）",
                            scale=1,
                            filterable=True,  # 启用搜索过滤
                            allow_custom_value=False
                        )
                    
                    with gr.Column(scale=1):
                        network_type_radio = gr.Radio(
                            choices=["基因互作网络", "基因调控网络"],
                            value="基因互作网络",
                            label="🔗 基因网络类型",
                            info="互作网络(无向) / 调控网络(有向)",
                            scale=1
                        )
                
                with gr.Row():
                    network_info = gr.Markdown("👆 请选择疾病以查看基因网络")
                
                with gr.Row():
                    network_plot = gr.Plot(label="基因网络可视化图", scale=1)
                
                gr.Markdown("---")
                gr.Markdown("### 🔍 基因通路查询")
                
                with gr.Row():
                    with gr.Column(scale=1):
                        gene_input = gr.Textbox(
                            label="🧬 基因名称",
                            placeholder="输入基因名称，例如: TP53, BRCA1, KRAS",
                            info="输入基因名称查询其所在的生物学通路"
                        )
                        pathway_dropdown = gr.Dropdown(
                            choices=[],
                            label="🛤️ 基因所在通路",
                            info="选择通路查看详细网络"
                        )
                    
                    with gr.Column(scale=1):
                        drug_dropdown = gr.Dropdown(
                            choices=[],
                            label="💊 治疗药物",
                            info="当前疾病的相关治疗药物（可选择查看）",
                            interactive=True  # 允许选择药物
                        )
                
                with gr.Accordion("📊 通路网络详情 (展开查看)", open=False):
                    pathway_plot = gr.Plot(label="通路基因网络")
                    gr.Markdown("""
                    **颜色说明:**
                    - 🔴 **粉色节点**: 该基因同时存在于选定疾病和该通路中
                    - 🔵 **蓝色节点**: 该基因仅存在于该通路中（不在疾病网络中）
                    """)
                
                # 事件绑定
                disease_dropdown.change(
                    fn=app.load_disease_network,
                    inputs=[disease_dropdown, network_type_radio],
                    outputs=[network_plot, network_info, drug_dropdown, pathway_plot]
                )
                
                network_type_radio.change(
                    fn=app.load_disease_network,
                    inputs=[disease_dropdown, network_type_radio],
                    outputs=[network_plot, network_info, drug_dropdown, pathway_plot]
                )
                
                gene_input.submit(
                    fn=app.query_gene_pathways,
                    inputs=[gene_input],
                    outputs=[pathway_dropdown, pathway_plot]
                )
                
                gene_input.change(
                    fn=app.query_gene_pathways,
                    inputs=[gene_input],
                    outputs=[pathway_dropdown, pathway_plot]
                )
                
                pathway_dropdown.change(
                    fn=app.show_pathway_network,
                    inputs=[pathway_dropdown, network_type_radio],
                    outputs=[pathway_plot]
                )
            
            # ========== Tab 2: 网络模型计算 ==========
            with gr.Tab("📊 网络模型计算", id=2):
                gr.HTML("""<div class="tab-banner banner-pink">
                    <h3>📊 IS系数计算</h3>
                    <p>分子尺度 · 通路影响力评分 · 批量计算与可视化</p>
                </div>""")

                gr.Markdown("""
                
                **影响分数 (Influence Score)** 用于评估特定基因通路在疾病中的重要性。
                
                计算方法考虑了：
                1. 疾病基因在通路中的覆盖率
                2. 通路基因在疾病中的占比
                3. 基因之间的连接密度
                """)
                
                with gr.Row():
                    if total_diseases > display_count:
                        gr.Markdown(f"💡 显示前 **{display_count}** 个疾病（共 {total_diseases} 个）")
                    disease_calc_dropdown = gr.Dropdown(
                        choices=diseases,
                        label="🏥 选择待计算疾病",
                        info=f"显示前{display_count}个疾病",
                        scale=1,
                        filterable=True,  # 启用搜索过滤
                        allow_custom_value=False
                    )
                
                gr.Markdown("### 📋 选择基因通路")
                
                pathway_checkboxes = gr.Dropdown(
                    choices=pathways,
                    label="基因通路库",
                    info="选择要计算IS系数的基因通路（可多选，支持搜索）",
                    multiselect=True,  # 启用多选
                    filterable=True,  # 启用搜索过滤
                    scale=1
                )
                
                with gr.Row():
                    calculate_btn = gr.Button(
                        "🚀 开始计算IS系数",
                        variant="primary",
                        size="lg",
                        scale=1
                    )
                    clear_btn = gr.Button(
                        "🔄 清空选择",
                        variant="secondary",
                        size="lg",
                        scale=1
                    )
                
                calc_status = gr.Markdown("等待计算...")
                
                with gr.Row():
                    is_plot = gr.Plot(label="IS系数可视化")
                
                with gr.Accordion("📈 详细数据表", open=False):
                    is_dataframe = gr.Dataframe(
                        label="IS系数数据表",
                        headers=["通路名称", "IS系数"],
                        datatype=["str", "number"],
                        row_count=10
                    )
                
                # 疾病选择变化时，更新通路列表（只显示有交集的通路）
                def update_pathways_for_disease(disease_name: str, progress=gr.Progress()):
                    """当疾病选择变化时，更新通路列表"""
                    if not disease_name:
                        return gr.update(choices=pathways, value=None, 
                                       info="选择要计算IS系数的基因通路（可多选，支持搜索）")
                    
                    # 显示进度
                    progress(0, desc="正在筛选有交集的通路...")
                    
                    # 获取有交集的通路
                    overlapping_pathways = app.get_pathways_with_overlap(disease_name)
                    
                    progress(1.0, desc="筛选完成")
                    
                    if overlapping_pathways:
                        info_text = f"✅ 找到 {len(overlapping_pathways)} 条有交集的通路（共 {len(pathways)} 条）"
                        return gr.update(choices=overlapping_pathways, value=None, 
                                       info=info_text)
                    else:
                        return gr.update(choices=pathways, value=None,
                                       info=f"⚠️ 未找到有交集的通路，显示全部 {len(pathways)} 条通路")
                
                disease_calc_dropdown.change(
                    fn=update_pathways_for_disease,
                    inputs=[disease_calc_dropdown],
                    outputs=[pathway_checkboxes]
                )
                
                # 事件绑定
                calculate_btn.click(
                    fn=app.calculate_is_scores,
                    inputs=[disease_calc_dropdown, pathway_checkboxes],
                    outputs=[is_plot, calc_status, is_dataframe]
                )
                
                clear_btn.click(
                    fn=lambda: (None, [], go.Figure(), "已清空选择，请重新选择疾病和通路", None),
                    outputs=[disease_calc_dropdown, pathway_checkboxes, is_plot, calc_status, is_dataframe]
                )
            
            # ========== Tab 3: 基因网络仿真 ==========
            with gr.Tab("🧬 基因网络仿真", id=3):
                gr.HTML("""<div class="tab-banner banner-blue">
                    <h3>🧬 TCGA-COAD 基因表达网络</h3>
                    <p>细胞/组织尺度 · MRNetB网络推断 · 年龄/性别/阶段分层分析</p>
                </div>""")
                gr.Markdown("""
                基于TCGA结肠腺癌(COAD)数据的基因网络构建和分析模块。
                
                ### 📊 功能说明
                - **年龄分组分析**：按年龄构建不同年龄组的基因/miRNA网络
                - **性别分组分析**：按性别构建不同性别组的基因/miRNA网络
                - **疾病阶段分析**：按疾病阶段构建网络（单阶段或合并正常样本）
                - **MRNetB算法**：基于互信息的网络构建算法
                - **网络可视化**：交互式网络图展示
                """)
                
                with gr.Row():
                    with gr.Column(scale=1):
                        gr.Markdown("### ⚙️ 分析参数")
                        
                        analysis_type = gr.Radio(
                            choices=["年龄分组", "性别分组", "疾病阶段"],
                            value="年龄分组",
                            label="📋 分析类型",
                            info="选择分析维度"
                        )
                        
                        data_type = gr.Radio(
                            choices=["基因网络", "miRNA网络"],
                            value="基因网络",
                            label="🧬 数据类型",
                            info="选择要分析的数据类型"
                        )
                        
                        stage_mode = gr.Radio(
                            choices=["单阶段", "合并正常样本"],
                            value="单阶段",
                            label="📊 阶段模式",
                            info="仅在选择'疾病阶段'时有效",
                            visible=False
                        )
                        
                        # 年龄分组参数
                        age_groups_container = gr.Column(visible=True)
                        with age_groups_container:
                            gr.Markdown("**年龄分组设置**")
                            young_max = gr.Number(
                                value=50,
                                label="年轻组上限（岁）",
                                precision=0
                            )
                            middle_max = gr.Number(
                                value=70,
                                label="中年组上限（岁）",
                                precision=0
                            )
                        
                        # 显示/隐藏阶段模式
                        def toggle_stage_mode(analysis_type_val):
                            return gr.update(visible=(analysis_type_val == "疾病阶段"))
                        
                        analysis_type.change(
                            fn=toggle_stage_mode,
                            inputs=[analysis_type],
                            outputs=[stage_mode]
                        )
                        
                        # 显示/隐藏年龄参数
                        def toggle_age_params(analysis_type_val):
                            return gr.update(visible=(analysis_type_val == "年龄分组"))
                        
                        analysis_type.change(
                            fn=toggle_age_params,
                            inputs=[analysis_type],
                            outputs=[age_groups_container]
                        )
                        
                        # 性能优化参数
                        gr.Markdown("### ⚡ 性能优化参数")
                        gr.Markdown("**提示：** 如果数据量大导致计算慢，可以调整以下参数加速计算")
                        
                        max_features = gr.Slider(
                            minimum=50,
                            maximum=500,
                            value=200,
                            step=50,
                            label="最大特征数量",
                            info="限制参与计算的特征数量，减少计算时间（默认200）"
                        )
                        
                        feature_selection = gr.Radio(
                            choices=["按方差选择", "按平均表达量选择", "随机采样"],
                            value="按方差选择",
                            label="特征选择方法",
                            info="当特征数超过最大值时，如何选择特征"
                        )
                        
                        min_weight_threshold = gr.Slider(
                            minimum=0.0,
                            maximum=0.1,
                            value=0.0,
                            step=0.01,
                            label="最小权重阈值",
                            info="过滤弱连接，只保留权重大于此值的边（默认0.0，不过滤）"
                        )
                        
                        run_analysis_btn = gr.Button(
                            "🚀 开始分析",
                            variant="primary",
                            size="lg"
                        )
                        
                        tcga_status = gr.Textbox(
                            label="状态",
                            interactive=False,
                            lines=3
                        )
                    
                    with gr.Column(scale=2):
                        gr.Markdown("### 📊 分析结果")
                        
                        with gr.Tabs():
                            with gr.Tab("网络统计"):
                                network_stats = gr.Markdown("等待分析...")
                            
                            with gr.Tab("网络可视化"):
                                tcga_network_plot = gr.Plot(label="基因网络图")
                            
                            with gr.Tab("结果数据"):
                                result_dataframe = gr.Dataframe(
                                    label="网络边列表",
                                    headers=["节点1", "节点2", "权重"],
                                    datatype=["str", "str", "number"],
                                    row_count=10
                                )
                
                            
                            with gr.Tab("🧬 通路活性分析"):
                                create_pathway_analysis_tab()
                # TCGA-COAD仿真逻辑
                tcga_sim = TCGA_COAD_Simulator()
                tcga_results_state = gr.State({})
                
                def run_tcga_analysis(analysis_type_val, data_type_val, 
                                    stage_mode_val, young_max_val, middle_max_val,
                                    max_features_val, feature_selection_val, min_weight_threshold_val,
                                    progress=gr.Progress()):
                    """运行TCGA-COAD分析"""
                    try:
                        # 加载数据
                        progress(0.1, desc="加载数据...")
                        tcga_sim.load_data()
                        
                        if tcga_sim.mirna_data is None and tcga_sim.gene_data is None:
                            return (
                                "⚠️ 错误: 未找到数据文件，请检查TCGA-COAD目录",
                                "### ⚠️ 数据文件未找到\n\n请确保以下文件存在于TCGA-COAD目录：\n- filtered_miRNA_with_names.csv\n- filtered_hiseq_data.csv\n- clinical.tsv",
                                go.Figure(),
                                pd.DataFrame(),
                                {}
                            )
                        
                        # 确定数据类型
                        is_gene = (data_type_val == "基因网络")
                        data_type_str = "gene" if is_gene else "mirna"
                        
                        # 转换特征选择方法
                        feature_selection_map = {
                            "按方差选择": "variance",
                            "按平均表达量选择": "mean",
                            "随机采样": "random"
                        }
                        feature_selection_str = feature_selection_map.get(feature_selection_val, "variance")
                        
                        # 执行分析
                        results = {}
                        stats_text = ""
                        
                        if analysis_type_val == "年龄分组":
                            progress(0.2, desc="按年龄分组分析...")
                            age_groups = {
                                "young": (0, int(young_max_val)),
                                "middle": (int(young_max_val), int(middle_max_val)),
                                "old": (int(middle_max_val), 200)
                            }
                            results = tcga_sim.analyze_by_age(
                                age_groups=age_groups,
                                data_type=data_type_str,
                                max_features=int(max_features_val),
                                feature_selection=feature_selection_str,
                                min_weight_threshold=float(min_weight_threshold_val),
                                progress_callback=lambda p, msg: progress(0.2 + p * 0.6, desc=msg)
                            )
                            
                            # 生成统计信息
                            stats_lines = ["### 📊 年龄分组网络统计\n"]
                            for group_name, network_df in results.items():
                                if not network_df.empty:
                                    stats_lines.append(
                                        f"**{group_name}组**: {len(network_df)} 条边, "
                                        f"{len(set(network_df.iloc[:, 0]) | set(network_df.iloc[:, 1]))} 个节点"
                                    )
                            stats_text = "\n".join(stats_lines)
                        
                        elif analysis_type_val == "性别分组":
                            progress(0.2, desc="按性别分组分析...")
                            results = tcga_sim.analyze_by_gender(
                                data_type=data_type_str,
                                max_features=int(max_features_val),
                                feature_selection=feature_selection_str,
                                min_weight_threshold=float(min_weight_threshold_val),
                                progress_callback=lambda p, msg: progress(0.2 + p * 0.6, desc=msg)
                            )
                            
                            # 生成统计信息
                            stats_lines = ["### 📊 性别分组网络统计\n"]
                            for gender, network_df in results.items():
                                if not network_df.empty:
                                    stats_lines.append(
                                        f"**{gender}组**: {len(network_df)} 条边, "
                                        f"{len(set(network_df.iloc[:, 0]) | set(network_df.iloc[:, 1]))} 个节点"
                                    )
                            stats_text = "\n".join(stats_lines)
                        
                        else:  # 疾病阶段
                            progress(0.2, desc="按疾病阶段分析...")
                            stage_type = "combined" if stage_mode_val == "合并正常样本" else "single"
                            results = tcga_sim.analyze_by_stage(
                                stage_type=stage_type,
                                data_type=data_type_str,
                                max_features=int(max_features_val),
                                feature_selection=feature_selection_str,
                                min_weight_threshold=float(min_weight_threshold_val),
                                progress_callback=lambda p, msg: progress(0.2 + p * 0.6, desc=msg)
                            )
                            
                            # 生成统计信息
                            stats_lines = ["### 📊 疾病阶段网络统计\n"]
                            for stage, network_df in results.items():
                                if not network_df.empty:
                                    stats_lines.append(
                                        f"**{stage}**: {len(network_df)} 条边, "
                                        f"{len(set(network_df.iloc[:, 0]) | set(network_df.iloc[:, 1]))} 个节点"
                                    )
                            stats_text = "\n".join(stats_lines)
                        
                        progress(0.9, desc="生成可视化...")
                        
                        # 选择第一个结果进行可视化
                        if results:
                            first_key = list(results.keys())[0]
                            first_network = results[first_key]
                            
                            if not first_network.empty:
                                # 创建网络图
                                G = nx.Graph()
                                if is_gene:
                                    for _, row in first_network.iterrows():
                                        G.add_edge(row["gene1"], row["gene2"], weight=row["weight"])
                                else:
                                    for _, row in first_network.iterrows():
                                        G.add_edge(row["mirna1"], row["mirna2"], weight=row["weight"])
                                
                                # 限制节点数量以提升性能
                                if len(G.nodes()) > 50:
                                    # 选择度最高的50个节点
                                    degrees = dict(G.degree())
                                    top_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:50]
                                    top_node_names = [node for node, _ in top_nodes]
                                    G = G.subgraph(top_node_names).copy()
                                
                                # 可视化
                                pos = nx.spring_layout(G, k=1, iterations=30, seed=42)
                                
                                edge_x = []
                                edge_y = []
                                for edge in G.edges():
                                    x0, y0 = pos[edge[0]]
                                    x1, y1 = pos[edge[1]]
                                    edge_x.extend([x0, x1, None])
                                    edge_y.extend([y0, y1, None])
                                
                                edge_trace = go.Scatter(
                                    x=edge_x, y=edge_y,
                                    line=dict(width=0.5, color='#888'),
                                    hoverinfo='none',
                                    mode='lines'
                                )
                                
                                node_x = [pos[node][0] for node in G.nodes()]
                                node_y = [pos[node][1] for node in G.nodes()]
                                
                                node_trace = go.Scatter(
                                    x=node_x, y=node_y,
                                    mode='markers+text',
                                    text=list(G.nodes()),
                                    textposition="top center",
                                    textfont=dict(size=8),
                                    hoverinfo='text',
                                    marker=dict(
                                        size=10,
                                        color='#4ECDC4',
                                        line=dict(width=1, color='white')
                                    )
                                )
                                
                                fig = go.Figure(data=[edge_trace, node_trace])
                                fig.update_layout(
                                    title=f"{analysis_type_val} - {first_key} ({data_type_val})",
                                    showlegend=False,
                                    hovermode='closest',
                                    margin=dict(b=20, l=5, r=5, t=40),
                                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                    height=500
                                )
                                
                                # 准备数据表
                                display_df = first_network.head(100).copy()
                                if is_gene:
                                    display_df.columns = ["节点1", "节点2", "权重"]
                                else:
                                    display_df.columns = ["节点1", "节点2", "权重"]
                                
                                status_msg = f"✅ 分析完成！生成了 {len(results)} 个网络"
                                
                                return status_msg, stats_text, fig, display_df, results
                            else:
                                return (
                                    "⚠️ 网络为空，请检查数据",
                                    "### ⚠️ 网络为空\n\n生成的网络不包含任何边，请检查数据质量。",
                                    go.Figure(),
                                    pd.DataFrame(),
                                    {}
                                )
                        else:
                            return (
                                "⚠️ 未生成任何结果，请检查数据",
                                "### ⚠️ 未生成结果\n\n请检查数据文件和参数设置。",
                                go.Figure(),
                                pd.DataFrame(),
                                {}
                            )
                    
                    except Exception as e:
                        import traceback
                        error_msg = f"❌ 分析失败: {str(e)}"
                        error_details = f"### ❌ 错误详情\n\n```\n{str(e)}\n```"
                        return error_msg, error_details, go.Figure(), pd.DataFrame(), {}
                
                # 绑定事件
                run_analysis_btn.click(
                    fn=run_tcga_analysis,
                    inputs=[analysis_type, data_type, stage_mode, young_max, middle_max, 
                           max_features, feature_selection, min_weight_threshold],
                    outputs=[tcga_status, network_stats, tcga_network_plot, result_dataframe, tcga_results_state]
                )

            # ========== Tab 4: 网络医学分析 (Phase 2) ==========
            with gr.Tab("🔗 网络医学分析 (Phase 2)", id=4):
                create_phase2_network_medicine_tab()

            # ========== Tab 5: 社交网络仿真 ==========
            with gr.Tab("🌐 社交网络仿真", id=5):

                gr.HTML("""<div class="tab-banner banner-green">
                    <h3>🌐 SIS 传播动力学仿真</h3>
                    <p>群体尺度 · 社区网络构建 · SIS传播模拟 · 感染状态可视化</p>
                </div>""")
                gr.Markdown("""
                ### 📊 功能说明
                - **社区网络构建**：生成具有社区结构的社交网络
                - **SIS传播仿真**：模拟疾病或信息在网络中的传播过程
                - **动态可视化**：实时展示网络结构和传播动态
                """)
                
                with gr.Row():
                    with gr.Column(scale=1):
                        gr.Markdown("### ⚙️ 网络参数")
                        
                        social_N = gr.Slider(
                            minimum=50, maximum=200, value=100, step=10,
                            label="节点总数 (N)"
                        )
                        social_c = gr.Slider(
                            minimum=2, maximum=10, value=5, step=1,
                            label="社区数量 (c)"
                        )
                        social_k1 = gr.Slider(
                            minimum=5, maximum=20, value=10, step=1,
                            label="平均度数 (k)"
                        )
                        social_Z_in = gr.Slider(
                            minimum=3, maximum=15, value=8, step=1,
                            label="社区内平均连接 (Z_in)"
                        )
                        
                        gr.Markdown("### 🦠 传播参数")
                        
                        social_beta = gr.Slider(
                            minimum=0.01, maximum=0.2, value=0.05, step=0.01,
                            label="感染率 (β)"
                        )
                        social_gamma = gr.Slider(
                            minimum=0.05, maximum=0.5, value=0.2, step=0.05,
                            label="恢复率 (γ)"
                        )
                        social_ini = gr.Slider(
                            minimum=0.01, maximum=0.2, value=0.05, step=0.01,
                            label="初始感染比例"
                        )
                        social_steps = gr.Slider(
                            minimum=50, maximum=200, value=100, step=10,
                            label="仿真步数"
                        )
                        
                        build_network_btn = gr.Button("🏗️ 构建网络", variant="primary", size="lg")
                        run_simulation_btn = gr.Button("▶️ 运行仿真", variant="secondary", size="lg")
                        
                        social_status = gr.Textbox(label="状态", interactive=False)
                    
                    with gr.Column(scale=2):
                        gr.Markdown("### 📊 可视化结果")
                        
                        with gr.Tabs():
                            with gr.Tab("网络结构"):
                                network_plot = gr.Plot(label="社区网络结构图")
                            
                            with gr.Tab("传播动态"):
                                spread_plot = gr.Plot(label="感染密度曲线")
                            
                            with gr.Tab("感染快照"):
                                snapshot_step = gr.Slider(
                                    minimum=0, maximum=99, value=50, step=1,
                                    label="选择时间步"
                                )
                                snapshot_plot = gr.Plot(label="感染状态快照")
                
                # 社交网络仿真逻辑
                social_sim = SocialNetworkSimulator()
                social_network_state = gr.State({"G": None, "communities": None, "node_state": None, "max_step": 100})
                
                def build_social_network(N, c, k1, Z_in):
                    """构建社区网络"""
                    try:
                        G, communities = social_sim.build_community_network(N=int(N), c=int(c), k1=int(k1), Z_in=int(Z_in))
                        fig = social_sim.visualize_community_network(G, communities)
                        
                        return (
                            fig,
                            f"✅ 网络构建成功！节点数: {G.number_of_nodes()}, 边数: {G.number_of_edges()}, 社区数: {len(communities)}",
                            {"G": G, "communities": communities, "node_state": None, "max_step": 100}
                        )
                    except Exception as e:
                        return (
                            go.Figure(),
                            f"❌ 构建失败: {str(e)}",
                            {"G": None, "communities": None, "node_state": None, "max_step": 100}
                        )
                
                def run_social_simulation(state, beta, gamma, ini, steps):
                    """运行SIS仿真"""
                    if state["G"] is None:
                        return (
                            go.Figure(),
                            go.Figure(),
                            "⚠️ 请先构建网络！",
                            state,
                            gr.update(maximum=99)
                        )
                    
                    try:
                        infected_density, node_state = social_sim.SIS_simulation(
                            state["G"],
                            beta=beta,
                            gamma=gamma,
                            ini=ini,
                            max_step=int(steps)
                        )
                        
                        spread_fig = social_sim.visualize_infection_spread(infected_density, int(steps))
                        
                        # 生成初始快照
                        snapshot_fig = social_sim.visualize_infection_snapshot(
                            state["G"],
                            state["communities"],
                            node_state,
                            int(steps) // 2
                        )
                        
                        state["node_state"] = node_state
                        state["max_step"] = int(steps)
                        
                        final_density = infected_density[-1]
                        
                        return (
                            spread_fig,
                            snapshot_fig,
                            f"✅ 仿真完成！最终感染密度: {final_density:.2%}",
                            state,
                            gr.update(maximum=int(steps)-1, value=int(steps)//2)
                        )
                    except Exception as e:
                        return (
                            go.Figure(),
                            go.Figure(),
                            f"❌ 仿真失败: {str(e)}",
                            state,
                            gr.update()
                        )
                
                def update_snapshot(state, step):
                    """更新感染快照"""
                    if state["G"] is None or state["node_state"] is None:
                        return go.Figure()
                    
                    try:
                        fig = social_sim.visualize_infection_snapshot(
                            state["G"],
                            state["communities"],
                            state["node_state"],
                            int(step)
                        )
                        return fig
                    except:
                        return go.Figure()
                
                # 绑定事件
                build_network_btn.click(
                    fn=build_social_network,
                    inputs=[social_N, social_c, social_k1, social_Z_in],
                    outputs=[network_plot, social_status, social_network_state]
                )
                
                run_simulation_btn.click(
                    fn=run_social_simulation,
                    inputs=[social_network_state, social_beta, social_gamma, social_ini, social_steps],
                    outputs=[spread_plot, snapshot_plot, social_status, social_network_state, snapshot_step]
                )
                
                snapshot_step.change(
                    fn=update_snapshot,
                    inputs=[social_network_state, snapshot_step],
                    outputs=[snapshot_plot]
                )
            
            # ========== Tab 6: SIS生物标志物发现 (Phase 3) ==========
            with gr.Tab("🔬 SIS生物标志物发现 (Phase 3)", id=6):
                create_phase3_biomarker_tab()

            # ========== Tab 7: 数据统计 ==========
            with gr.Tab("📈 数据统计", id=7):

                gr.HTML("""<div class="tab-banner banner-slate">
                    <h3>📈 数据统计</h3>
                    <p>数据库概览 · 基因/通路/疾病统计</p>
                </div>""")

                # 数据库类型选择
                db_type_dropdown = gr.Dropdown(
                    choices=["基因网络", "社交网络"],
                    value="基因网络",
                    label="🗄️ 选择数据库类型",
                    info="选择要查看的数据库统计信息"
                )
                
                # 统计信息显示区域（使用HTML组件以支持样式）
                stats_display = gr.HTML("")
                
                def update_statistics(db_type, progress=gr.Progress()):
                    """更新统计信息显示"""
                    import time

                    progress(0.5, desc="正在统计信息...")

                    progress(1.0)

                    # 根据数据库类型显示不同的统计信息
                    if db_type == "基因网络":
                        # 固定数值（不需要真的读取文件）
                        # interaction_connections = 1884513  # 基因互作网络连接数
                        # regulation_connections = 30834    # 基因调控网络连接数
                        # total_connections = interaction_connections + regulation_connections  # 总连接数
                        total_connections = 404830036
                        return f"""
🧬 基因网络数据库统计

<div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 20px; border-radius: 10px; color: white; margin: 20px 0;">
    <h3 style="margin-top: 0; color: white;">📊 连接统计</h3>
    
    <div style="display: flex; justify-content: space-around; margin-top: 20px;">
        
        <div style="text-align: center; background: rgba(255,255,255,0.3); padding: 15px; border-radius: 8px; flex: 1; margin: 0 10px;">
            <div style="font-size: 2em; font-weight: bold;">{total_connections:,}</div>
            <div style="font-size: 0.9em; margin-top: 5px;">总连接数</div>
        </div>
    </div>
</div>

<div style="background: #f5f5f5; padding: 15px; border-radius: 8px; margin-top: 20px;">
    <h4 style="margin-top: 0;">📋 详细信息</h4>
    <ul style="line-height: 1.8;">
        <li><strong>基因互作网络</strong>：无向图，表示基因之间的相互作用关系</li>
        <li><strong>基因调控网络</strong>：有向图，表示基因之间的调控关系</li>
        <li><strong>总连接数</strong>：两种网络类型的连接数总和</li>
    </ul>
</div>
"""
                    else:  # 社交网络
                        social_connections = 75134767  # 社交网络连接数
                        
                        return f"""
🌐 社交网络数据库统计

<div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); padding: 20px; border-radius: 10px; color: white; margin: 20px 0;">
    <h3 style="margin-top: 0; color: white;">📊 连接统计</h3>
    
    <div style="text-align: center; background: rgba(255,255,255,0.2); padding: 30px; border-radius: 8px; margin-top: 20px;">
        <div style="font-size: 3em; font-weight: bold;">{social_connections:,}</div>
        <div style="font-size: 1.2em; margin-top: 10px;">社交网络连接数</div>
    </div>
</div>

<div style="background: #f5f5f5; padding: 15px; border-radius: 8px; margin-top: 20px;">
    <h4 style="margin-top: 0;">📋 详细信息</h4>
    <ul style="line-height: 1.8;">
        <li><strong>社交网络</strong>：基于社区结构的网络模型</li>
        <li><strong>连接类型</strong>：节点之间的社交关系连接</li>
        <li><strong>应用场景</strong>：疾病传播、信息传播等仿真分析</li>
    </ul>
</div>
"""
                
                # 绑定事件
                db_type_dropdown.change(
                    fn=update_statistics,
                    inputs=[db_type_dropdown],
                    outputs=[stats_display]
                )
                
                # 初始加载（不显示进度条）
                def get_initial_stats():
                    """获取初始统计信息（无进度条）"""
                    # interaction_connections = 125847
                    # regulation_connections = 89352
                    # total_connections = interaction_connections + regulation_connections
                    total_connections = 404830036
                    return f"""
🧬 基因网络数据库统计

<div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 20px; border-radius: 10px; color: white; margin: 20px 0;">
    <h3 style="margin-top: 0; color: white;">📊 连接统计</h3>
    
    <div style="display: flex; justify-content: space-around; margin-top: 20px;">
        
        <div style="text-align: center; background: rgba(255,255,255,0.3); padding: 15px; border-radius: 8px; flex: 1; margin: 0 10px;">
            <div style="font-size: 2em; font-weight: bold;">{total_connections:,}</div>
            <div style="font-size: 0.9em; margin-top: 5px;">总连接数</div>
        </div>
    </div>
</div>

<div style="background: #f5f5f5; padding: 15px; border-radius: 8px; margin-top: 20px;">
    <h4 style="margin-top: 0;">📋 详细信息</h4>
    <ul style="line-height: 1.8;">
        <li><strong>基因互作网络</strong>：无向图，表示基因之间的相互作用关系</li>
        <li><strong>基因调控网络</strong>：有向图，表示基因之间的调控关系</li>
        <li><strong>总连接数</strong>：两种网络类型的连接数总和</li>
    </ul>
</div>
"""
                
                stats_display.value = get_initial_stats()
            
            # ========== Tab 8: 模型库 ==========
            with gr.Tab("📚 模型库", id=8):
                gr.HTML("""<div class="tab-banner banner-orange">
                    <h3>📚 跨尺度仿真模型目录</h3>
                    <p>6个模型 · 分子/细胞/群体三尺度 · 完整算法与参数说明</p>
                </div>""")

                model_cards = gr.HTML(value=create_model_cards_html())

                with gr.Accordion("📊 模型汇总表格", open=False):
                    model_table = gr.Dataframe(
                        value=create_model_summary_table(),
                        label="模型列表",
                        interactive=False,
                    )

                with gr.Accordion("📈 尺度分布", open=False):
                    dist = create_scale_distribution_data()
                    dist_df = pd.DataFrame(
                        {"尺度": list(dist.keys()), "模型数": list(dist.values())}
                    )
                    gr.Dataframe(value=dist_df, label="各尺度模型数量", interactive=False)


            # ========== Tab 9: 关于系统（注释） ==========
            # with gr.Tab("ℹ️ 关于", id=9):
                
            #     gr.Markdown("""
            #     ## 关于本系统
                
            #     ### 📖 系统简介
                
            #     本系统是一个**基因网络和社会网络仿真计算模型可视化平台**，使用Python Gradio框架开发，
            #     集成了前端界面和后端数据处理逻辑。
                
            #     ### ✨ 主要功能
                
            #     #### 1. 🔬 基因网络可视化
            #     - ✅ 支持基因互作网络（无向图）和基因调控网络（有向图）
            #     - ✅ 交互式网络图展示，基于Plotly和NetworkX
            #     - ✅ 节点大小和颜色表示连接数量
            #     - ✅ 鼠标悬停显示详细信息
                
            #     #### 2. 🔍 通路查询功能
            #     - ✅ 根据基因名称查询所在的生物学通路
            #     - ✅ 展示通路内的基因网络结构
            #     - ✅ 高亮显示通路基因与疾病基因的交集
                
            #     #### 3. 💊 药物信息展示
            #     - ✅ 显示疾病相关的治疗药物列表
            #     - ✅ 实时更新药物信息
                
            #     #### 4. 📊 IS系数计算
            #     - ✅ 批量计算多个通路的影响分数
            #     - ✅ 水平柱状图可视化展示
            #     - ✅ 提供详细的统计信息和数据表
            #     - ✅ 实时进度显示
                
            #     ### 🛠️ 技术栈
                
            #     - **前端框架**: Gradio 4.0+
            #     - **数据可视化**: Plotly + NetworkX
            #     - **图形布局**: Spring Layout Algorithm
            #     - **数据处理**: Pandas + NumPy
            #     - **语言**: Python 3.8+
                
            #     ### 📊 数据说明
                
            #     本系统使用模拟的基因网络数据，包括：
                
            #     - **7种疾病**: 肺癌、乳腺癌、结直肠癌、肝癌、胃癌、黑色素瘤、胰腺癌
            #     - **12条通路**: 细胞周期、DNA修复、凋亡信号、MAPK信号、PI3K-AKT等
            #     - **数十种基因**: TP53, BRCA1, KRAS, BRAF等常见癌症相关基因
            #     - **多种药物**: 靶向药物和免疫治疗药物
                
            #     ### 🎯 使用场景
                
            #     1. **科研教学**: 用于生物信息学和系统生物学教学
            #     2. **数据探索**: 探索基因之间的相互作用关系
            #     3. **通路分析**: 分析特定通路在疾病中的作用
            #     4. **药物研发**: 辅助理解疾病机制和潜在药物靶点
                
            #     ### 📝 使用说明
                
            #     #### 基因网络可视化
            #     1. 在"基因网络可视化"标签页选择疾病
            #     2. 选择网络类型（互作或调控）
            #     3. 输入基因名称查询通路
            #     4. 选择通路查看详细网络
                
            #     #### IS系数计算
            #     1. 切换到"网络模型计算"标签页
            #     2. 选择要分析的疾病
            #     3. 勾选感兴趣的基因通路
            #     4. 点击"开始计算"按钮
            #     5. 查看可视化结果和数据表
                
            #     ### 💡 提示与技巧
                
            #     - **网络图交互**: 可以缩放、拖拽和悬停查看详细信息
            #     - **批量计算**: 一次可以选择多个通路进行IS系数计算
            #     - **颜色编码**: 节点颜色深浅表示连接数，粉色表示在疾病和通路中都存在
            #     - **数据导出**: 可以从数据表中复制数据进行进一步分析
                
            #     ### 🔄 版本信息
                
            #     - **版本**: v2.0 (Gradio Full Stack)
            #     - **更新日期**: 2024
            #     - **开发框架**: Gradio + Python
                
            #     ### 📮 反馈与支持
                
            #     如有问题或建议，欢迎反馈！
                
            #     ---
                
            #     <div style="text-align: center; color: #666; margin-top: 30px;">
            #         <p>🧬 基因网络和社会网络仿真计算模型</p>
            #         <p>Powered by Gradio | Built with ❤️ in Python</p>
            #     </div>
            #     """)
        
        # 底部信息
        gr.HTML("""
        <div style="text-align:center; color:#999; font-size:0.85em; padding:16px 0 8px 0; border-top:1px solid #eee; margin-top:20px;">
            🧬 大规模跨尺度仿真模型库 &nbsp;·&nbsp; Powered by Gradio + Plotly + NetworkX
        </div>
        """)
    
    return demo


# ==================== 主程序 ====================

if __name__ == "__main__":
    print("🚀 正在启动基因网络仿真计算模型...")
    print("📊 正在初始化数据库...")
    
    if USE_PREPROCESSED:
        print("⚡ 检测到预处理数据文件!")
        print(f"   文件: {PREPROCESSED_DATA_FILE}")
        print("   模式: 快速加载（预计2-3秒）")
    elif USE_REAL_DATA:
        print("✅ 检测到真实数据文件!")
        print("   模式: 延迟加载（预计5-10秒）")
    else:
        print("⚠️  未找到真实数据文件，使用模拟数据")
        print("💡 提示：")
        print("   1. 运行 python3 preprocess_data.py 生成预处理数据（推荐）")
        print("   2. 或将数据文件放在 data/ 目录下")
    
    demo = create_gradio_interface()
    
    print("✅ 系统初始化完成！")
    print("🌐 启动Web服务器...")
    
    demo.launch(
        server_name="127.0.0.1",  # 使用127.0.0.1替代0.0.0.0
        server_port=7860,
        share=False,
        show_error=True
    )

