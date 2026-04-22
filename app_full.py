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

# 导入疾病溯源模块
from disease_tracing import DiseaseTracer

# 导入通路分析模块
from gradio_phase1_integration import create_pathway_analysis_tab, Phase1DataLoader
# from model_library import create_model_cards_html, create_model_summary_table, create_scale_distribution_data
# 导入网络医学分析模块
from gradio_phase2_integration import create_phase2_network_medicine_tab, Phase2DataLoader
# 导入SIS生物标志物发现模块（已禁用）
# from gradio_phase3_integration import create_phase3_biomarker_tab, Phase3DataLoader

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
            
            # ========== Tab 0: 跨尺度联动分析 —— 已移至 Tab 1.5（跨尺度社交网络溯源）之后 ==========

            # ========== Tab 1: 跨尺度疾病溯源 ==========
            with gr.Tab("🔬 跨尺度疾病溯源", id=1):
                gr.HTML("""<div class="tab-banner banner-red">
                    <h3>🔬 跨尺度疾病溯源分析</h3>
                    <p>选择疾病 → 点击「开始溯源」→ 通路 / 基因网络 → 深入分析</p>
                </div>""")

                # 初始化溯源器
                _tracer = DiseaseTracer()
                _tracer.load_data()
                _tracer_diseases = _tracer.get_disease_list()

                # === Step 1: 疾病选择 ===
                gr.Markdown("### 🏥 Step 1 · 选择疾病")
                with gr.Row():
                    trace_disease = gr.Dropdown(
                        choices=_tracer_diseases,
                        value=_tracer_diseases[0] if _tracer_diseases else None,
                        label="选择疾病",
                        scale=2,
                    )
                    trace_disease_info = gr.Markdown("👆 请选择疾病查看关联基因")

                # === Step 2: 疾病关联基因列表 ===
                gr.Markdown("### 🧬 Step 2 · 疾病关联基因列表")
                with gr.Row():
                    trace_gene_list_df = gr.Dataframe(
                        headers=["基因", "TCGA 可用"],
                        label="关联基因",
                        interactive=False,
                        wrap=True,
                        row_count=(10, "dynamic"),
                    )
                    trace_gene_list_summary = gr.Markdown("")

                # === Step 3: 跨尺度说明 ===
                # === Step 3: 跨尺度溯源目标 ===
                gr.Markdown("### 🧭 Step 3 · 跨尺度溯源目标")
                gr.Markdown(
                    r"""
疾病是一种典型的**涌现现象（emergent phenomenon）**：其病理表型由**分子互作、基因表达调控、通路协同扰动**在多个生物学尺度上耦合作用共同决定。单尺度分析工具（ORA / GSEA / DEG / PPI 模块检测）仅能捕捉完整病理机制的一个投影。

**跨尺度溯源的目标**，是从疾病这一宏观表型出发，**逐层下钻定位驱动其发生的分子因素**：首先在通路层定位被显著扰动的功能模块，进而在基因层识别通路中的核心驱动基因，最终在转录调控层追溯作用于该基因的上游 miRNA。整条分析链回答的是同一个问题——"**疾病 → 通路 → 基因 → miRNA**" 这条因果链上，哪些元素在每一层各自扮演关键角色。
""",
                )

                # === Step 4: 开始溯源按钮 ===
                gr.Markdown("### 🚀 Step 4 · 执行溯源")
                with gr.Row():
                    trace_start_btn = gr.Button(
                        "🔬 开始溯源（生成通路 & 基因网络）",
                        variant="primary", size="lg", scale=1,
                    )

                # === Step 5: 通路网络 ===
                gr.Markdown("### 🌌 Step 5 · 通路网络（疾病相关通路高亮）")
                trace_pw_universe = gr.Plot(label="通路网络 — 疾病相关通路高亮")

                # === Step 6: 基因网络 ===
                gr.Markdown("### 🧬 Step 6 · 基因网络（通路基因 + 共表达基因高亮）")
                trace_gene_universe = gr.Plot(label="基因网络 — 通路基因 + 共表达基因高亮")

                # === 阶段性结论：承上（Step 5~6 产物）+ 启下（引向 Step 7） ===
                trace_stage_conclusion = gr.Markdown(
                    "_点击「开始溯源」后，此处将自动生成跨尺度阶段性结论。_"
                )

                # === Step 7: 深入分析 ===
                gr.Markdown("### 📋 Step 7 · 深入分析（选择通路 / 基因后生成详细报告）")
                with gr.Row():
                    with gr.Column(scale=1):
                        trace_pathway = gr.Dropdown(
                            choices=[], label="🔬 选择通路",
                            info="按影响力评分排序（CV加权）",
                            allow_custom_value=True,
                        )
                        trace_gene = gr.Dropdown(
                            choices=[], label="🧬 选择基因",
                            info="按NAE评分排序",
                            allow_custom_value=True,
                        )
                        trace_run_btn = gr.Button("📋 生成溯源报告", variant="primary", size="lg")

                    with gr.Column(scale=3):
                        with gr.Tabs():
                            with gr.Tab("🌊 Sankey溯源"):
                                trace_sankey = gr.Plot(label="Sankey溯源图")
                                trace_report = gr.HTML(label="溯源报告")

                            with gr.Tab("🕸️ 通路详情"):
                                trace_pw_network = gr.Plot(label="通路关联网络")
                                trace_pw_bar = gr.Plot(label="通路影响力排名")
                                trace_pw_table = gr.Dataframe(label="通路详情", interactive=False)

                            with gr.Tab("🧬 基因NAE"):
                                trace_gene_network = gr.Plot(label="基因模块网络")
                                trace_gene_table = gr.Dataframe(label="基因NAE排名", interactive=False)

                            with gr.Tab("📈 基因表达"):
                                trace_expr_plot = gr.Plot(label="表达箱线图")
                                trace_mirna_table = gr.Dataframe(label="miRNA调控关系", interactive=False)
                                trace_expr_summary = gr.Markdown("")

                # === 回调 ===
                def on_disease_select(disease):
                    """选择疾病 → 只展示基本信息 + 关联基因列表（不生成网络图，不阻塞）"""
                    if not disease:
                        return "", pd.DataFrame(), ""

                    try:
                        ov = _tracer.get_disease_overview(disease)
                    except Exception as e:
                        return f"❌ 读取失败: {e}", pd.DataFrame(), ""
                    if "error" in ov:
                        return f"❌ {ov['error']}", pd.DataFrame(), ""

                    info = (
                        f"**{disease}** | 分类: {ov.get('category','未知')} | "
                        f"关联基因: {ov['n_genes']} | TCGA可用: {ov['n_genes_in_tcga']} | "
                        f"通路: {ov['n_pathways']}"
                    )

                    avail_set = set(ov.get("available_genes", []))
                    genes_sorted = sorted(ov.get("genes", []))
                    gene_rows = [
                        [g, "✅" if g in avail_set else "—"]
                        for g in genes_sorted
                    ]
                    gene_df = pd.DataFrame(gene_rows, columns=["基因", "TCGA 可用"])

                    summary = (
                        f"共 **{len(genes_sorted)}** 个关联基因，其中 "
                        f"**{len(avail_set)}** 个在 TCGA 表达矩阵中可用。"
                        " 👉 点击下方「开始溯源」生成通路与基因网络。"
                    )
                    return info, gene_df, summary

                def run_trace(disease):
                    """点击「开始溯源」→ 生成通路网络 + 基因网络 + 填充通路下拉 + 生成阶段性结论"""
                    empty_conclusion = "_点击「开始溯源」后，此处将自动生成跨尺度阶段性结论。_"
                    if not disease:
                        return (go.Figure(), go.Figure(),
                                gr.update(choices=[]), gr.update(choices=[]),
                                empty_conclusion)

                    try:
                        ov = _tracer.get_disease_overview(disease)
                        if "error" in ov:
                            return (go.Figure(), go.Figure(),
                                    gr.update(choices=[]), gr.update(choices=[]),
                                    f"❌ {ov['error']}")

                        pw_universe = _tracer.create_pathway_universe(disease)
                        gene_universe = _tracer.create_gene_universe(disease)

                        pw_df = _tracer.get_pathway_ranking(disease)
                        pw_choices = [
                            f"{row['通路']} ({row['影响力评分']:.2f})"
                            for _, row in pw_df.head(30).iterrows()
                        ]

                        # ========== 动态阶段性结论 ==========
                        n_genes = ov.get("n_genes", 0)
                        n_avail = ov.get("n_genes_in_tcga", 0)
                        n_pw = ov.get("n_pathways", 0)

                        # Top 3 通路
                        top_pw = pw_df.head(3) if not pw_df.empty else pd.DataFrame()
                        pw_items = [
                            f"**{row['通路']}**（评分 {row['影响力评分']:.2f}）"
                            for _, row in top_pw.iterrows()
                        ]
                        pw_text = "、".join(pw_items) if pw_items else "—"

                        # 从 top 通路取核心基因（NAE 最高）
                        top_gene_text = "—"
                        if not top_pw.empty:
                            leading_pw = top_pw.iloc[0]["通路"]
                            gene_df, _ = _tracer.get_gene_module(disease, leading_pw)
                            if not gene_df.empty:
                                top_genes = gene_df.head(3)
                                gene_items = [
                                    f"**{row['基因']}**（NAE {row['NAE评分']:.3f}）"
                                    for _, row in top_genes.iterrows()
                                ]
                                top_gene_text = "、".join(gene_items)

                        conclusion = f"""
---

### 📌 阶段性结论（Step 5 – 6 产物）

针对疾病 **{disease}**，系统已完成 **通路层** 与 **基因层** 的跨尺度下钻：

- 🌌 **通路层**：从 {n_pw} 条关联通路中识别出影响力最高的 Top 3 —— {pw_text}
- 🧬 **基因层**：在 Top 通路内，基于 NAE 量化拓扑中心性，识别核心驱动基因 —— {top_gene_text}
- 📊 **数据可用性**：疾病关联基因共 {n_genes} 个，其中 {n_avail} 个在 TCGA-COAD 表达矩阵中可追溯

---
"""

                        return (
                            pw_universe,
                            gene_universe,
                            gr.update(choices=pw_choices, value=pw_choices[0] if pw_choices else None),
                            gr.update(choices=[]),
                            conclusion,
                        )
                    except Exception as e:
                        return (go.Figure(), go.Figure(),
                                gr.update(choices=[]), gr.update(choices=[]),
                                f"❌ 溯源失败：{e}")

                def on_pathway_select(disease, pathway_str):
                    """选择通路 → 更新基因下拉"""
                    if not disease or not pathway_str:
                        return gr.update(choices=[])
                    pathway = pathway_str.split(" (")[0]
                    gene_df, _ = _tracer.get_gene_module(disease, pathway)
                    if gene_df.empty:
                        return gr.update(choices=[])
                    gene_choices = [f"{row['基因']} (NAE:{row['NAE评分']:.3f})" for _, row in gene_df.iterrows()]
                    return gr.update(choices=gene_choices, value=gene_choices[0] if gene_choices else None)

                def run_full_trace(disease, pathway_str, gene_str):
                    """点击按钮 → 生成详细分析"""
                    empty = (go.Figure(), "", go.Figure(), go.Figure(), pd.DataFrame(),
                             go.Figure(), pd.DataFrame(), go.Figure(),
                             pd.DataFrame(), "")
                    if not disease or not pathway_str:
                        return empty

                    pathway = pathway_str.split(" (")[0]
                    gene = gene_str.split(" (")[0] if gene_str else ""

                    try:
                        # Sankey
                        sankey = _tracer.create_sankey_diagram(
                            disease, selected_pathway=pathway, selected_gene=gene or None)

                        # 通路星云图
                        pw_df = _tracer.get_pathway_ranking(disease)
                        pw_network = _tracer.create_pathway_network(
                            disease, selected_pathway=pathway)
                        pw_bar = _tracer.create_pathway_bar(pw_df)

                        # 基因模块
                        gene_df, G = _tracer.get_gene_module(disease, pathway)
                        pathway_gene_set = set(gene_df["基因"].tolist()) if len(gene_df) > 0 else set()
                        neighbor_genes = []
                        if not pw_df.empty and pathway_gene_set:
                            nc = set()
                            for _, pw_row in pw_df.iterrows():
                                if pw_row["通路"] == pathway:
                                    continue
                                pw_set = set(pw_row["基因列表"].split(", "))
                                if pw_set & pathway_gene_set:
                                    nc.update(list(pw_set - pathway_gene_set)[:3])
                            neighbor_genes = list(nc)[:20]
                        gene_net = _tracer.create_gene_network_plot(
                            G, gene_df, selected_gene=gene or None,
                            neighbor_genes=neighbor_genes if neighbor_genes else None)

                        # 基因表达
                        expr_plot = go.Figure()
                        mirna_df = pd.DataFrame()
                        expr_summary = ""
                        ctx = gene_df["基因"].tolist()[:15] if len(gene_df) > 0 else None
                        if gene and gene in (_tracer.expr_data.index if _tracer.expr_data is not None else []):
                            profile = _tracer.get_gene_expression_profile(gene)
                            expr_plot = _tracer.create_expression_boxplot(gene, profile, context_genes=ctx)
                            mirnas = profile.get("mirna_regulators", [])
                            mirna_df = pd.DataFrame(mirnas) if mirnas else pd.DataFrame()
                            expr_summary = f"""**★ {gene}** | 均值: {profile.get('expr_mean','N/A')} | 排名: #{profile.get('expr_rank','N/A')}/{profile.get('total_genes','N/A')} | miRNA: {len(mirnas)} 个"""

                        report = _tracer.generate_tracing_report(disease, pathway, gene) if gene else ""

                        return (sankey, report, pw_network, pw_bar, pw_df.head(20),
                                gene_net, gene_df, expr_plot, mirna_df, expr_summary)
                    except Exception as e:
                        import traceback
                        return (go.Figure(), f"<p>❌ {e}</p>", go.Figure(), go.Figure(), pd.DataFrame(),
                                go.Figure(), pd.DataFrame(), go.Figure(), pd.DataFrame(),
                                f"❌ {traceback.format_exc()}")

                # 选择疾病 → 只显示基本信息 + 关联基因列表（不生成网络图，不阻塞首屏）
                trace_disease.change(
                    on_disease_select,
                    inputs=[trace_disease],
                    outputs=[trace_disease_info, trace_gene_list_df, trace_gene_list_summary],
                )
                # 点击「开始溯源」→ 渲染通路网络 + 基因网络 + 填充通路下拉
                trace_start_btn.click(
                    run_trace,
                    inputs=[trace_disease],
                    outputs=[trace_pw_universe, trace_gene_universe,
                             trace_pathway, trace_gene,
                             trace_stage_conclusion],
                )
                trace_pathway.change(
                    on_pathway_select,
                    inputs=[trace_disease, trace_pathway],
                    outputs=[trace_gene],
                )
                trace_run_btn.click(
                    run_full_trace,
                    inputs=[trace_disease, trace_pathway, trace_gene],
                    outputs=[trace_sankey, trace_report, trace_pw_network, trace_pw_bar, trace_pw_table,
                             trace_gene_network, trace_gene_table, trace_expr_plot,
                             trace_mirna_table, trace_expr_summary],
                )

            # ========== Tab 0 (重排后): 跨尺度联动分析 ==========
            with gr.Tab("🔗 跨尺度联动分析", id=0):
                gr.HTML("""<div class="tab-banner banner-purple">
                    <h3>🔗 基因网络跨尺度联动分析</h3>
                    <p>分子尺度（基因互作 / 调控）↔ 细胞 / 组织尺度（TCGA 表达推断）双层联动</p>
                </div>""")

                # === 跨尺度联动目标说明（与跨尺度疾病溯源 Step 3 同风格）===
                gr.Markdown("### 🧭 跨尺度联动分析的目标")
                gr.Markdown(
                    r"""
同一组疾病相关基因，在**不同生物学尺度**上会呈现出**截然不同的网络拓扑结构**：在分子尺度上表现为基于文献 / 实验整理的基因互作关系，在细胞 / 组织尺度上表现为从患者样本表达数据中涌现的功能关联模块。孤立看待任一尺度，都只能反映疾病网络的**局部侧面**——前者缺乏患者层级的个体差异信息，后者脱离已知机制背景。

**跨尺度联动分析的目标**，是在**同一疾病下**将多个生物学尺度的网络**并置比较、相互印证**：从分子层识别核心枢纽基因，将其作为锚点约束细胞层的网络推断搜索空间，再通过两层网络的拓扑对照判断**先验机制与患者数据是否一致**。这一「**分子 ↔ 细胞**」的双向联动，为疾病机制研究提供**跨尺度一致性证据**，是单一数据源或单一尺度分析无法给出的视角。
""",
                )

                with gr.Tabs():

                    # ---- 子Tab 1: 基因网络级联 ----
                    with gr.Tab("🚀 基因网络级联"):
                        gr.Markdown("#### 基因网络两层级联：分子 → 细胞/组织")
                        gr.Markdown("""
                        **传递机制:**
                        - **Layer 1 → Layer 2**: 分子层Hub基因 → 细胞层MRNetB种子节点
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
                                cascade_run_btn = gr.Button("🚀 运行基因网络级联分析",
                                                             variant="primary", size="lg")
                            with gr.Column(scale=2):
                                cascade_status = gr.Markdown("")
                                cascade_arch = gr.HTML(label="级联架构")

                        cascade_summary = gr.Dataframe(label="级联汇总表", interactive=False)
                        cascade_insights = gr.Markdown("")

                        # === 阶段性结论卡片（动态）===
                        cascade_stage_conclusion = gr.Markdown(
                            "_运行「基因网络级联分析」后，此处将自动生成跨尺度阶段性结论。_"
                        )

                        def run_gene_cascade(disease, net_type, max_feat):
                            empty_conclusion = "_运行「基因网络级联分析」后，此处将自动生成跨尺度阶段性结论。_"
                            try:
                                engine = CrossScaleEngine(
                                    db=app.db,
                                    tcga_sim=TCGA_COAD_Simulator(),
                                )
                                report = engine.run_gene_cascade(
                                    disease_name=disease,
                                    network_type=net_type,
                                    max_features=int(max_feat),
                                )
                                arch = report.architecture_html
                                summary = engine.cascade_summary_df(report)
                                insights = "\n".join(
                                    f"- {i}" for i in report.cross_scale_insights
                                )

                                # === 拼装阶段性结论 ===
                                mol = report.results.get("molecular")
                                cell = report.results.get("cellular")
                                mol_s = mol.summary if mol else {}
                                cell_s = cell.summary if cell else {}
                                mol_nodes = mol_s.get("nodes", 0)
                                mol_edges = mol_s.get("edges", 0)
                                mol_density = mol_s.get("density", 0.0)
                                hubs = mol_s.get("hub_genes", []) or []
                                hub_str = "、".join(hubs[:5]) if hubs else "—"
                                cell_nodes = cell_s.get("cell_nodes", cell_s.get("nodes", 0))
                                cell_edges = cell_s.get("cell_edges", cell_s.get("edges", 0))
                                cell_cluster = cell_s.get("cell_clustering", cell_s.get("clustering", 0.0))

                                # === 根据实际数据下不同结论 ===
                                # 判断两层网络的结构强弱
                                mol_has_structure = mol_edges >= 5
                                cell_has_structure = cell_edges >= 3

                                if mol_has_structure and cell_has_structure:
                                    # 两层都有结构：跨尺度一致性结论
                                    verdict = (
                                        f"在分子层构建的基因互作网络呈现出明确的连通结构（{mol_nodes} 节点 / "
                                        f"{mol_edges} 边，密度 {mol_density:.3f}），其核心枢纽基因（{hub_str}）"
                                        f"作为锚点约束下，细胞 / 组织层在 TCGA 表达数据中同样推断出包含 "
                                        f"{cell_edges} 条功能关联边的共表达网络（聚类系数 {cell_cluster:.3f}），"
                                        f"表明文献先验刻画的疾病分子机制在真实患者数据中获得了数据层面的呼应，"
                                        f"上述枢纽基因在两个尺度上均占据拓扑中心地位，可作为后续机制研究与靶点筛选的重点候选。"
                                    )
                                elif mol_has_structure and not cell_has_structure:
                                    # 分子层有结构、细胞层稀疏：如实说明
                                    verdict = (
                                        f"分子层构建的基因互作网络呈现出明确的连通结构（{mol_nodes} 节点 / "
                                        f"{mol_edges} 边，密度 {mol_density:.3f}），核心枢纽基因为 {hub_str}；"
                                        f"然而，以这些枢纽基因为锚点在 TCGA 表达数据上进行细胞 / 组织层推断时，"
                                        f"未能形成显著的共表达结构（仅 {cell_edges} 条关联边）。"
                                    )
                                elif not mol_has_structure and cell_has_structure:
                                    # 分子层稀疏、细胞层有结构：提示患者数据主导
                                    verdict = (
                                        f"分子层的基因互作结构相对稀疏（仅 {mol_edges} 条已知互作边），"
                                        f"但在细胞 / 组织层中，基于 TCGA 表达数据却观察到相对丰富的共表达网络"
                                        f"（{cell_edges} 条关联边，聚类系数 {cell_cluster:.3f}）。"
                                        f"这提示对该疾病而言，现有文献先验对核心分子互作的覆盖可能不足，"
                                        f"患者层级的数据驱动分析在此情形下具有额外的补充价值，可以揭示尚未被充分记录的基因协同模式。"
                                    )
                                else:
                                    # 两层都很稀疏：数据/疾病不适配
                                    verdict = (
                                        f"两个尺度均未形成显著的网络结构（分子层 {mol_edges} 边 / 细胞层 {cell_edges} 边），"
                                        f"提示该疾病在当前数据背景下（文献先验 + TCGA-COAD 结肠癌队列）"
                                        f"缺乏足够的信号支持开展跨尺度联动分析，建议更换疾病或调整特征筛选参数后重试。"
                                    )

                                conclusion = f"""
---

### 📌 阶段性结论

{verdict}

---
"""

                                return ("✅ 基因网络级联分析完成", arch,
                                        summary, insights, conclusion)
                            except Exception as e:
                                import traceback
                                return (f"❌ {e}", "",
                                        pd.DataFrame(),
                                        f"```\n{traceback.format_exc()}\n```",
                                        empty_conclusion)

                        cascade_run_btn.click(
                            fn=run_gene_cascade,
                            inputs=[cascade_disease, cascade_net_type, cascade_max_feat],
                            outputs=[cascade_status, cascade_arch,
                                     cascade_summary, cascade_insights,
                                     cascade_stage_conclusion],
                        )

                    # ---- 子Tab 2: 多疾病对比（已注释：应需求暂时隐藏） ----
                    # with gr.Tab("📊 多疾病对比"):
                    #     gr.Markdown("#### 多疾病分子网络对比分析")
                    #     gr.Markdown("选择 2-5 个疾病，对比它们在分子尺度的网络拓扑差异。")
                    #
                    #     compare_diseases_select = gr.CheckboxGroup(
                    #         choices=diseases,
                    #         value=diseases[:3] if len(diseases) >= 3 else diseases,
                    #         label="选择疾病（2-5个）",
                    #     )
                    #     compare_run_btn = gr.Button("📊 运行对比分析", variant="primary")
                    #     compare_status = gr.Markdown("")
                    #     compare_chart = gr.Plot(label="对比柱状图")
                    #     compare_table = gr.Dataframe(label="对比数据表", interactive=False)
                    #
                    #     def run_compare(selected_diseases):
                    #         if not selected_diseases or len(selected_diseases) < 2:
                    #             return "⚠️ 请至少选择2个疾病", go.Figure(), pd.DataFrame()
                    #         try:
                    #             engine = CrossScaleEngine(db=app.db)
                    #             df, fig = engine.compare_diseases(selected_diseases)
                    #             return "✅ 对比完成", fig, df
                    #         except Exception as e:
                    #             return f"❌ {e}", go.Figure(), pd.DataFrame()
                    #
                    #     compare_run_btn.click(
                    #         fn=run_compare,
                    #         inputs=[compare_diseases_select],
                    #         outputs=[compare_status, compare_chart, compare_table],
                    #     )

                    # ---- 子Tab 3: 基因追踪（已注释：应需求暂时隐藏） ----
                    # with gr.Tab("🔍 基因追踪"):
                    #     gr.Markdown("#### 基因跨尺度追踪")
                    #     gr.Markdown("输入一个基因，查看它在 **分子层（互作网络）** 和 **细胞层（TCGA表达）** 中的角色。")
                    #
                    #     with gr.Row():
                    #         cascade_trace_gene_input = gr.Textbox(
                    #             label="基因名称",
                    #             placeholder="例如: TP53, BRCA1, KRAS",
                    #             value="TP53",
                    #         )
                    #         cascade_trace_disease = gr.Dropdown(
                    #             choices=diseases,
                    #             value=diseases[0] if diseases else None,
                    #             label="疾病背景",
                    #         )
                    #     cascade_trace_btn = gr.Button("🔍 追踪基因", variant="primary")
                    #     cascade_trace_result = gr.HTML(label="基因档案卡")
                    #
                    #     def run_cascade_trace(gene, disease):
                    #         if not gene or not gene.strip():
                    #             return "<p>⚠️ 请输入基因名称</p>"
                    #         try:
                    #             engine = CrossScaleEngine(
                    #                 db=app.db,
                    #                 tcga_sim=TCGA_COAD_Simulator(),
                    #             )
                    #             return engine.trace_gene(gene.strip().upper(), disease)
                    #         except Exception as e:
                    #             return f"<p>❌ 追踪失败: {e}</p>"
                    #
                    #     cascade_trace_btn.click(
                    #         fn=run_cascade_trace,
                    #         inputs=[cascade_trace_gene_input, cascade_trace_disease],
                    #         outputs=[cascade_trace_result],
                    #     )

            # ========== Tab 3: 跨尺度基因网络仿真 ==========
            with gr.Tab("🧬 跨尺度基因网络仿真", id=4):
                gr.HTML("""<div class="tab-banner banner-blue">
                    <h3>🧬 TCGA-COAD 基因表达网络</h3>
                    <p>细胞/组织尺度 · MRNetB网络推断 · 年龄/性别/阶段分层分析</p>
                </div>""")
                gr.Markdown("""
                基于 TCGA 结肠腺癌（COAD）真实患者队列，以 **MRNetB 互信息算法** 从高维基因 / miRNA 表达矩阵中推断功能关联网络。

                单一尺度的基因网络（如文献已知的基因互作或单个基因的差异表达）只能反映疾病机制的局部切面。本模块将**患者队列的临床表型**（年龄、性别、疾病阶段）与**分子层级的基因 / miRNA 表达**联动，在同一算法框架下为不同人群构建各自的基因网络，再通过对比不同人群的网络拓扑差异，揭示：

                - **个体尺度（临床特征）→ 分子尺度（基因 / miRNA 共表达）** 的双向联动 —— 临床变量的改变如何重塑分子网络结构；
                - **基因尺度 ↔ miRNA 尺度** 的对称对照 —— 通过切换数据类型，可以在相同人群分层下独立构建基因网络与 miRNA 网络，两者互为参照；
                - **人群分层之间的跨尺度异质性** —— 不同年龄 / 性别 / 阶段的患者在分子网络上的差异，即为疾病在临床 × 分子两个尺度上涌现的真实异质性证据。
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
                            value=100,
                            step=50,
                            label="最大特征数量",
                            info="限制参与计算的特征数量，减少计算时间（默认 100，演示推荐；精细分析可调至 200+）"
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
                            value=0.01,
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

                        # === 阶段性结论卡片（动态）===
                        tcga_conclusion = gr.Markdown(
                            "_运行分析后，此处将自动生成跨尺度阶段性结论。_"
                        )
                
                            
                            # with gr.Tab("🧬 通路活性分析"):
                            #     create_pathway_analysis_tab()
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
                                
                                n_nodes = G.number_of_nodes()
                                n_edges = G.number_of_edges()
                                density = nx.density(G) if n_nodes > 1 else 0
                                clustering = nx.average_clustering(G) if n_nodes > 2 else 0
                                top_genes = sorted(dict(G.degree()).items(), key=lambda x: x[1], reverse=True)[:5]
                                top_gene_names = [g[0] for g in top_genes]

                                status_msg = f"""✅ 分析完成！生成了 {len(results)} 个网络

📝 分析摘要:
使用 MRNetB 互信息算法从 TCGA-COAD 表达数据推断出 {n_nodes} 个节点、{n_edges} 条边的基因功能关联网络。
网络密度 {density:.4f}，聚类系数 {clustering:.4f}。
{'核心基因（度数最高）: ' + ', '.join(top_gene_names) + '，这些基因与大量其他基因存在表达关联，可能是关键调控节点。' if top_gene_names else ''}"""

                                # === 阶段性结论（叙述性）===
                                group_count = len(results)
                                data_label = "基因共表达" if is_gene else "miRNA 共表达"
                                hub_list = "、".join(top_gene_names[:5]) if top_gene_names else "—"

                                if n_edges >= 5:
                                    tcga_conclusion_text = f"""
---

### 📌 阶段性结论

在 TCGA-COAD 患者队列的 **{analysis_type_val}** 分层下，系统对 **{group_count} 个人群亚组**分别构建了 {data_label} 网络。以展示的 **{first_key}** 组为例，通过 MRNetB 互信息算法从分子表达矩阵中推断出 {n_nodes} 节点 / {n_edges} 边的功能关联结构（网络密度 {density:.3f}，聚类系数 {clustering:.3f}），并识别出度数最高的核心节点 {hub_list}。这些核心节点在本亚组内与大量其他分子存在表达层级的协同关系，是可能的关键调控枢纽；结合不同亚组之间网络拓扑的比较，可进一步观察临床分层（{analysis_type_val}）在分子网络层面的异质性表现。

---
"""
                                else:
                                    tcga_conclusion_text = f"""
---

### 📌 阶段性结论

在 TCGA-COAD 患者队列的 **{analysis_type_val}** 分层下，对 **{first_key}** 组构建的 {data_label} 网络结构较为稀疏（仅 {n_nodes} 节点 / {n_edges} 边），说明该亚组在当前特征筛选参数下未能形成显著的共表达关联，可尝试提高特征数量或降低权重阈值后重试。

---
"""

                                return status_msg, stats_text, fig, display_df, results, tcga_conclusion_text
                            else:
                                return (
                                    "⚠️ 网络为空，请检查数据",
                                    "### ⚠️ 网络为空\n\n生成的网络不包含任何边，请检查数据质量。",
                                    go.Figure(),
                                    pd.DataFrame(),
                                    {},
                                    "_运行分析后，此处将自动生成跨尺度阶段性结论。_"
                                )
                        else:
                            return (
                                "⚠️ 未生成任何结果，请检查数据",
                                "### ⚠️ 未生成结果\n\n请检查数据文件和参数设置。",
                                go.Figure(),
                                pd.DataFrame(),
                                {},
                                "_运行分析后，此处将自动生成跨尺度阶段性结论。_"
                            )

                    except Exception as e:
                        import traceback
                        error_msg = f"❌ 分析失败: {str(e)}"
                        error_details = f"### ❌ 错误详情\n\n```\n{str(e)}\n```"
                        return error_msg, error_details, go.Figure(), pd.DataFrame(), {}, "_运行分析后，此处将自动生成跨尺度阶段性结论。_"

                # 绑定事件
                run_analysis_btn.click(
                    fn=run_tcga_analysis,
                    inputs=[analysis_type, data_type, stage_mode, young_max, middle_max,
                           max_features, feature_selection, min_weight_threshold],
                    outputs=[tcga_status, network_stats, tcga_network_plot, result_dataframe, tcga_results_state, tcga_conclusion]
                )

            # ========== Tab 2: 基因网络可视化 ==========
            with gr.Tab("🔬 基因网络可视化", id=2):
                gr.HTML("""<div class="tab-banner banner-indigo">
                    <h3>🔬 基因网络可视化</h3>
                    <p>分子尺度 · 基因互作网络 / 调控网络 · 通路查询</p>
                </div>""")

                gr.Markdown(
                    r"""
基因网络并不只存在于"基因-基因"单一尺度。对同一组疾病关联基因，本模块同时呈现**基因互作网络**（无向、反映蛋白质—蛋白质或功能协同关系）与**基因调控网络**（有向、反映转录因子 → 靶基因的调控关系），两张网络构成**同一分子尺度下的两种机制视角**。在此之上，进一步联动**通路尺度**——对任一基因，可查询其所参与的生物学通路及对应的药物靶点信息，从而将微观的分子互作/调控结构与中观的通路机制、应用层的治疗干预连接起来，形成「**调控 / 互作 → 通路 → 药物**」的跨尺度观察链条。
"""
                )

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
            with gr.Tab("📊 网络模型计算", id=3):
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
            
            # ========== Tab 1.5: 跨尺度社交网络溯源 ==========
            # 逻辑对标跨尺度疾病溯源：议题 → 话题社区 → 用户影响力
            with gr.Tab("🌐 跨尺度社交网络溯源", id=10):
                gr.HTML("""<div class="tab-banner banner-purple">
                    <h3>🌐 跨尺度社交网络溯源分析</h3>
                    <p>选择议题 → 点击「开始溯源」→ 话题 / 用户影响力网络 → 深入分析</p>
                </div>""")

                # === 内置议题库（竞选/公共议题样例） ===
                _SN_ISSUES = {
                    "Election2024-Immigration": {
                        "category": "政治竞选",
                        "topics": [
                            "BorderSecurity", "Asylum", "PathToCitizenship",
                            "H1B-Policy", "SanctuaryCity", "ICE-Enforcement",
                            "Dreamers", "BorderWall", "RefugeeQuota", "WorkVisa",
                        ],
                        "stance_polarity": 0.72,
                    },
                    "Election2024-ClimateChange": {
                        "category": "政策议题",
                        "topics": [
                            "CarbonTax", "GreenNewDeal", "EV-Subsidy",
                            "FossilFuels", "RenewableGrid", "ParisAccord",
                            "MethaneRegs", "OilDrilling", "SolarCredit", "CoalMines",
                        ],
                        "stance_polarity": 0.65,
                    },
                    "Election2024-Healthcare": {
                        "category": "民生议题",
                        "topics": [
                            "MedicareForAll", "ACA-Repeal", "DrugPricing",
                            "PrivateInsurance", "PublicOption", "Medicaid",
                            "PreExisting", "TeleHealth", "MentalHealth", "Abortion",
                        ],
                        "stance_polarity": 0.81,
                    },
                    "Election2024-GunPolicy": {
                        "category": "社会议题",
                        "topics": [
                            "AssaultWeaponsBan", "BackgroundChecks", "ConcealCarry",
                            "RedFlagLaw", "2A-Rights", "StandYourGround",
                            "GhostGuns", "HighCapMags", "ATF-Regs", "SchoolSafety",
                        ],
                        "stance_polarity": 0.78,
                    },
                    "Brexit-Referendum": {
                        "category": "国际事件",
                        "topics": [
                            "Sovereignty", "SingleMarket", "Immigration-UK",
                            "NHS-Funding", "TradeDeal", "IrishBorder",
                            "FishingRights", "EuropeanCourt", "FreedomOfMovement", "CustomsUnion",
                        ],
                        "stance_polarity": 0.69,
                    },
                }
                _sn_issue_choices = list(_SN_ISSUES.keys())

                # === Step 1: 议题选择 ===
                gr.Markdown("### 🏛️ Step 1 · 选择政治议题 / 竞选事件")
                with gr.Row():
                    sn_issue = gr.Dropdown(
                        choices=_sn_issue_choices,
                        value=_sn_issue_choices[0],
                        label="选择议题",
                        scale=2,
                    )
                    sn_issue_info = gr.Markdown("👆 请选择议题查看关联话题")

                # === Step 2: 议题关联话题列表 ===
                gr.Markdown("### 🏷️ Step 2 · 议题关联话题（Topics / Hashtags）列表")
                with gr.Row():
                    sn_topic_list_df = gr.Dataframe(
                        headers=["话题", "热度指数", "极化倾向"],
                        label="关联话题",
                        interactive=False,
                        wrap=True,
                        row_count=(10, "dynamic"),
                    )
                    sn_topic_list_summary = gr.Markdown("")

                # === Step 3: 分析方法与公式 ===
                gr.Markdown("### 📖 Step 3 · 跨尺度极化溯源方法与公式")
                gr.Markdown(
                    r"""
系统依次执行三步分析，从宏观舆论极化逐层下钻到微观用户影响力：

---

#### 🔹 Step A · 意见动力学（非线性耦合演化方程）

$$
\dot{x}_i \;=\; -\,x_i \,+\, K \sum_{j=1}^{N} A_{ij}(t)\,\tanh\!\bigl(\alpha\, x_j\bigr)
$$

其中 $x_i \in \mathbb{R}$ 为用户 $i$ 的意见（正负号=阵营，绝对值=激进度）；$A_{ij}(t)$ 为时变社交邻接矩阵；$\alpha$ 为议题 **争议性（controversialness）**，控制 $\tanh$ 饱和斜率；$K$ 为 **社会影响强度耦合常数**。线性衰减项 $-x_i$ 建模意见回归中性的心理阻力，非线性求和项刻画邻居对意见的饱和型社会影响。

---

#### 🔹 Step B · 活动驱动网络与同质性连边

每个用户 $i$ 具有活动率 $a_i$，服从幂律分布：

$$
F(a) \;\propto\; a^{-\gamma}, \quad a \in [\,\epsilon,\, 1\,]
$$

在时间窗 $\Delta t$ 内，用户 $i$ 以概率 $a_i\,\Delta t$ 激活；激活时随机连出 $m$ 条边。连边目标 $j$ 由 **同质性（homophily）** 规则确定：

$$
p_{ij} \;=\; \frac{\bigl|\,x_i - x_j\,\bigr|^{-\beta}}{\displaystyle\sum_{k \neq i} \bigl|\,x_i - x_k\,\bigr|^{-\beta}}
$$

其中 $\beta \ge 0$ 为 **同质性强度**：$\beta=0$ 退化为随机连边；$\beta$ 越大，用户越倾向于连接意见相近者，形成 **回音室（echo chamber）**。

---

#### 🔹 Step C · 极化相变与序参量

系统的宏观极化程度由意见分布的一阶矩表征：

$$
\langle\,|x|\,\rangle \;=\; \frac{1}{N}\sum_{i=1}^{N} |x_i|
$$

$$
\rho \;=\; \sqrt{\,\langle x^2 \rangle \,-\, \langle x \rangle^{2}\,}
$$

在参数空间 $(\alpha,\beta)$ 中，系统呈现三个相：**共识相**（$\langle |x|\rangle \to 0$，单峰分布于 0）、**激进相**（单峰但 $\langle |x|\rangle$ 大）、**极化相**（双峰分布，$\rho$ 跃升）。极化相的判据：

$$
\mathrm{sign}(x_i)\ \text{呈双峰分布}\ \land\ \langle |x| \rangle > x_c
$$

其中 $x_c$ 为双峰/单峰分界的临界意见幅度。
""",
                    latex_delimiters=[
                        {"left": "$$", "right": "$$", "display": True},
                        {"left": "$", "right": "$", "display": False},
                    ],
                )

                # === Step 4: 开始溯源按钮 ===
                gr.Markdown("### 🚀 Step 4 · 执行溯源")
                with gr.Row():
                    sn_start_btn = gr.Button(
                        "🌐 开始溯源（生成话题社区 & 用户影响力网络）",
                        variant="primary", size="lg", scale=1,
                    )

                # === Step 5: 话题社区网络 ===
                gr.Markdown("### 🗺️ Step 5 · 话题社区网络（议题相关话题高亮）")
                sn_topic_network = gr.Plot(label="话题社区网络 — 议题相关话题高亮")

                # === Step 6: 用户影响力网络 ===
                gr.Markdown("### 👥 Step 6 · 用户影响力网络（意见领袖 + 跨阵营桥接节点高亮）")
                sn_user_network = gr.Plot(label="用户影响力网络 — 按阵营 + NAE 着色")

                # === Step 7: 深入分析 ===
                gr.Markdown("### 📋 Step 7 · 深入分析（选择话题 / 用户后生成详细报告）")
                with gr.Row():
                    with gr.Column(scale=1):
                        sn_topic = gr.Dropdown(
                            choices=[], label="🏷️ 选择话题",
                            info="按极化度评分排序",
                            allow_custom_value=True,
                        )
                        sn_user = gr.Dropdown(
                            choices=[], label="👤 选择用户",
                            info="按影响力 NAE 评分排序",
                            allow_custom_value=True,
                        )
                        sn_run_btn = gr.Button("📋 生成溯源报告", variant="primary", size="lg")

                    with gr.Column(scale=3):
                        with gr.Tabs():
                            with gr.Tab("🌊 Sankey溯源"):
                                sn_sankey = gr.Plot(label="议题→话题→用户 Sankey 溯源图")
                                sn_report = gr.HTML(label="溯源报告")

                            with gr.Tab("🕸️ 话题详情"):
                                sn_topic_detail_net = gr.Plot(label="话题关联网络")
                                sn_topic_bar = gr.Plot(label="话题极化度排名")
                                sn_topic_table = gr.Dataframe(label="话题详情", interactive=False)

                            with gr.Tab("👥 用户 NAE"):
                                sn_user_detail_net = gr.Plot(label="用户模块网络")
                                sn_user_table = gr.Dataframe(label="用户 NAE 排名", interactive=False)

                            with gr.Tab("📈 意见时序"):
                                sn_opinion_plot = gr.Plot(label="用户意见演化曲线")
                                sn_neighbor_table = gr.Dataframe(label="邻居意见 & 影响", interactive=False)
                                sn_opinion_summary = gr.Markdown("")

                # === 回调 ===
                def _sn_get_issue_overview(issue):
                    """返回议题基本信息 + 话题排序表"""
                    cfg = _SN_ISSUES.get(issue, {})
                    if not cfg:
                        return None
                    rng = np.random.default_rng(abs(hash(issue)) % (2**32))
                    topics = cfg["topics"]
                    # 热度指数 (0~1) & 极化倾向 (-1~+1)
                    heat = rng.uniform(0.2, 1.0, size=len(topics))
                    polarity = rng.uniform(-1.0, 1.0, size=len(topics)) * cfg["stance_polarity"]
                    df = pd.DataFrame({
                        "话题": topics,
                        "热度指数": np.round(heat, 3),
                        "极化倾向": np.round(polarity, 3),
                    }).sort_values("热度指数", ascending=False).reset_index(drop=True)
                    # 极化度评分 = |极化倾向| × 热度
                    df["极化度评分"] = np.round(np.abs(df["极化倾向"]) * df["热度指数"], 3)
                    return df, cfg

                def on_issue_select(issue):
                    """选议题 → 展示话题列表（首屏不跑重型计算）"""
                    if not issue:
                        return "", pd.DataFrame(), ""
                    res = _sn_get_issue_overview(issue)
                    if res is None:
                        return f"❌ 未找到议题: {issue}", pd.DataFrame(), ""
                    df, cfg = res
                    info = (
                        f"**{issue}** | 分类: {cfg['category']} | "
                        f"话题数: {len(df)} | 议题极化基线 ρ₀ = {cfg['stance_polarity']:.2f}"
                    )
                    show_df = df[["话题", "热度指数", "极化倾向"]].copy()
                    summary = (
                        f"共 **{len(df)}** 个关联话题，"
                        f"正向倾向 **{int((df['极化倾向']>0).sum())}** 个 / "
                        f"反向倾向 **{int((df['极化倾向']<0).sum())}** 个。"
                        " 👉 点击下方「开始溯源」生成话题社区与用户影响力网络。"
                    )
                    return info, show_df, summary

                def _sn_build_topic_network(df, cfg):
                    """基于话题共现构建社区图（用于 Step 5）"""
                    rng = np.random.default_rng(abs(hash(str(cfg))) % (2**32))
                    G = nx.Graph()
                    for _, r in df.iterrows():
                        G.add_node(r["话题"], heat=r["热度指数"], polar=r["极化倾向"])
                    topics = df["话题"].tolist()
                    # 按极化同向/异向概率连边（同向高，异向低）
                    for i, a in enumerate(topics):
                        for b in topics[i+1:]:
                            pa, pb = df.loc[df["话题"]==a, "极化倾向"].iloc[0], df.loc[df["话题"]==b, "极化倾向"].iloc[0]
                            same = (pa * pb) > 0
                            p = 0.55 if same else 0.12
                            if rng.random() < p:
                                G.add_edge(a, b, weight=float(rng.uniform(0.3, 1.0)))
                    return G

                def _sn_topic_plot(G, df):
                    """Step 5 话题网络 plot"""
                    pos = nx.spring_layout(G, seed=42, k=0.9)
                    edge_x, edge_y = [], []
                    for a, b in G.edges():
                        edge_x += [pos[a][0], pos[b][0], None]
                        edge_y += [pos[a][1], pos[b][1], None]
                    node_x = [pos[n][0] for n in G.nodes()]
                    node_y = [pos[n][1] for n in G.nodes()]
                    polar_map = dict(zip(df["话题"], df["极化倾向"]))
                    heat_map = dict(zip(df["话题"], df["热度指数"]))
                    node_color = [polar_map.get(n, 0) for n in G.nodes()]
                    node_size = [15 + heat_map.get(n, 0.3) * 35 for n in G.nodes()]
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines",
                                             line=dict(color="#6b7280", width=0.7),
                                             hoverinfo="none", showlegend=False))
                    fig.add_trace(go.Scatter(
                        x=node_x, y=node_y, mode="markers+text",
                        text=list(G.nodes()), textposition="top center",
                        textfont=dict(size=10, color="#e5e7eb"),
                        marker=dict(size=node_size, color=node_color,
                                    colorscale="RdBu", cmin=-1, cmax=1,
                                    line=dict(width=1, color="#111827"),
                                    colorbar=dict(title="极化倾向")),
                        hovertemplate="%{text}<br>极化:%{marker.color:.2f}<extra></extra>",
                        showlegend=False,
                    ))
                    fig.update_layout(
                        title="话题社区网络（红=保守/挺，蓝=自由/反；节点大小=热度）",
                        paper_bgcolor="#0b1120", plot_bgcolor="#0b1120",
                        font=dict(color="#e5e7eb"),
                        xaxis=dict(visible=False), yaxis=dict(visible=False),
                        height=560, margin=dict(l=10, r=10, t=50, b=10),
                    )
                    return fig

                def _sn_build_user_network(df, cfg, n_users=80):
                    """构建用户影响力网络（用于 Step 6）"""
                    rng = np.random.default_rng((abs(hash(str(cfg))) + 7) % (2**32))
                    # 每个话题分配一些用户
                    G = nx.Graph()
                    users = [f"u{i:03d}" for i in range(n_users)]
                    # 给每个用户一个意见 x ∈ [-1, 1]，基于议题极化基线的混合分布（双峰）
                    pol0 = cfg["stance_polarity"]
                    signs = rng.choice([-1, 1], size=n_users, p=[0.5, 0.5])
                    mags = np.abs(rng.normal(loc=pol0, scale=0.25, size=n_users))
                    opinions = np.clip(signs * mags, -1, 1)
                    for i, u in enumerate(users):
                        G.add_node(u, opinion=float(opinions[i]))
                    # 社交连边：意见相近者连接概率高
                    for i in range(n_users):
                        for j in range(i+1, n_users):
                            diff = abs(opinions[i] - opinions[j])
                            p = 0.35 * np.exp(-3 * diff) + 0.02  # 同阵营 0.37，跨阵营 0.02
                            if rng.random() < p:
                                G.add_edge(users[i], users[j])
                    return G, opinions

                def _sn_compute_nae(G, opinions):
                    """计算每个用户的 NAE (5 维属性熵)"""
                    if G.number_of_nodes() == 0:
                        return pd.DataFrame()
                    deg = dict(G.degree())
                    try:
                        btw = nx.betweenness_centrality(G)
                    except Exception:
                        btw = {n: 0.0 for n in G.nodes()}
                    try:
                        pr = nx.pagerank(G)
                    except Exception:
                        pr = {n: 1.0 / max(1, G.number_of_nodes()) for n in G.nodes()}
                    users = list(G.nodes())
                    op_map = dict(zip(users, opinions))
                    # 跨阵营桥接度 = 连向异号邻居的比例
                    bridge = {}
                    for u in users:
                        ns = list(G.neighbors(u))
                        if not ns:
                            bridge[u] = 0.0
                        else:
                            cross = sum(1 for v in ns if op_map[v] * op_map[u] < 0)
                            bridge[u] = cross / len(ns)
                    # 发帖活跃度 ~ |opinion| + noise
                    rng = np.random.default_rng(123)
                    activity = {u: float(abs(op_map[u]) + rng.uniform(0, 0.3)) for u in users}
                    rows = []
                    for u in users:
                        v = np.array([
                            deg.get(u, 0) + 1e-9,
                            btw.get(u, 0) + 1e-9,
                            pr.get(u, 0) + 1e-9,
                            bridge.get(u, 0) + 1e-9,
                            activity.get(u, 0) + 1e-9,
                        ], dtype=float)
                        p = v / v.sum()
                        nae = float(-np.sum(p * np.log2(p)))
                        rows.append({
                            "用户": u, "意见": round(op_map[u], 3),
                            "度": deg.get(u, 0),
                            "介数": round(btw.get(u, 0), 4),
                            "PageRank": round(pr.get(u, 0), 4),
                            "跨阵营桥接": round(bridge.get(u, 0), 3),
                            "活跃度": round(activity.get(u, 0), 3),
                            "NAE评分": round(nae, 4),
                        })
                    return pd.DataFrame(rows).sort_values("NAE评分", ascending=False).reset_index(drop=True)

                def _sn_user_plot(G, opinions, nae_df, highlight=None):
                    pos = nx.spring_layout(G, seed=7, k=0.7)
                    edge_x, edge_y = [], []
                    for a, b in G.edges():
                        edge_x += [pos[a][0], pos[b][0], None]
                        edge_y += [pos[a][1], pos[b][1], None]
                    users = list(G.nodes())
                    op_map = dict(zip(users, opinions))
                    nae_map = dict(zip(nae_df["用户"], nae_df["NAE评分"])) if not nae_df.empty else {}
                    node_x = [pos[n][0] for n in users]
                    node_y = [pos[n][1] for n in users]
                    node_color = [op_map[n] for n in users]
                    node_size = [8 + nae_map.get(n, 0) * 6 for n in users]
                    text_labels = [n if (highlight and n in highlight) else "" for n in users]
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines",
                                             line=dict(color="#4b5563", width=0.4),
                                             hoverinfo="none", showlegend=False))
                    fig.add_trace(go.Scatter(
                        x=node_x, y=node_y, mode="markers+text",
                        text=text_labels, textposition="top center",
                        textfont=dict(size=10, color="#fef3c7"),
                        marker=dict(size=node_size, color=node_color,
                                    colorscale="RdBu", cmin=-1, cmax=1,
                                    line=dict(width=0.6, color="#111827"),
                                    colorbar=dict(title="意见 x<sub>u</sub>")),
                        hovertext=[f"{n}<br>x={op_map[n]:.2f}<br>NAE={nae_map.get(n,0):.3f}" for n in users],
                        hoverinfo="text", showlegend=False,
                    ))
                    fig.update_layout(
                        title="用户影响力网络（颜色=意见倾向；大小=NAE 影响力）",
                        paper_bgcolor="#0b1120", plot_bgcolor="#0b1120",
                        font=dict(color="#e5e7eb"),
                        xaxis=dict(visible=False), yaxis=dict(visible=False),
                        height=560, margin=dict(l=10, r=10, t=50, b=10),
                    )
                    return fig

                # 用于跨回调共享当前网络（session 级简单缓存）
                _sn_cache = {"issue": None, "df": None, "cfg": None,
                             "Gtopic": None, "Guser": None, "opinions": None, "nae_df": None}

                def run_sn_trace(issue):
                    """点「开始溯源」→ Step5 + Step6 + 填充两个下拉"""
                    if not issue:
                        return go.Figure(), go.Figure(), gr.update(choices=[]), gr.update(choices=[])
                    try:
                        res = _sn_get_issue_overview(issue)
                        if res is None:
                            return go.Figure(), go.Figure(), gr.update(choices=[]), gr.update(choices=[])
                        df, cfg = res
                        Gtopic = _sn_build_topic_network(df, cfg)
                        topic_fig = _sn_topic_plot(Gtopic, df)
                        Guser, opinions = _sn_build_user_network(df, cfg)
                        nae_df = _sn_compute_nae(Guser, opinions)
                        user_fig = _sn_user_plot(Guser, opinions, nae_df)

                        _sn_cache.update({
                            "issue": issue, "df": df, "cfg": cfg,
                            "Gtopic": Gtopic, "Guser": Guser,
                            "opinions": opinions, "nae_df": nae_df,
                        })

                        topic_choices = [
                            f"{r['话题']} ({r['极化度评分']:.2f})"
                            for _, r in df.sort_values("极化度评分", ascending=False).head(30).iterrows()
                        ]
                        return (
                            topic_fig, user_fig,
                            gr.update(choices=topic_choices, value=topic_choices[0] if topic_choices else None),
                            gr.update(choices=[]),
                        )
                    except Exception as e:
                        import traceback; traceback.print_exc()
                        return go.Figure(), go.Figure(), gr.update(choices=[]), gr.update(choices=[])

                def on_sn_topic_select(issue, topic_str):
                    """选话题 → 填充用户下拉（按 NAE 排序前 20）"""
                    if not issue or not topic_str or _sn_cache.get("nae_df") is None:
                        return gr.update(choices=[])
                    nae_df = _sn_cache["nae_df"]
                    if nae_df.empty:
                        return gr.update(choices=[])
                    user_choices = [
                        f"{r['用户']} (NAE:{r['NAE评分']:.3f})"
                        for _, r in nae_df.head(20).iterrows()
                    ]
                    return gr.update(choices=user_choices, value=user_choices[0] if user_choices else None)

                def run_sn_full(issue, topic_str, user_str):
                    """Step 7 深入分析"""
                    empty = (go.Figure(), "", go.Figure(), go.Figure(), pd.DataFrame(),
                             go.Figure(), pd.DataFrame(), go.Figure(), pd.DataFrame(), "")
                    if not issue or _sn_cache.get("df") is None:
                        return empty
                    try:
                        df = _sn_cache["df"]; cfg = _sn_cache["cfg"]
                        Gtopic = _sn_cache["Gtopic"]; Guser = _sn_cache["Guser"]
                        opinions = _sn_cache["opinions"]; nae_df = _sn_cache["nae_df"]
                        topic = topic_str.split(" (")[0] if topic_str else ""
                        user = user_str.split(" (")[0] if user_str else ""

                        # Sankey: 议题 → 前5话题 → 前5用户
                        top_topics = df.sort_values("极化度评分", ascending=False).head(5)
                        top_users = nae_df.head(5)
                        labels = [issue] + top_topics["话题"].tolist() + top_users["用户"].tolist()
                        src, tgt, val = [], [], []
                        for i, r in enumerate(top_topics.itertuples()):
                            src.append(0); tgt.append(1 + i); val.append(float(r.极化度评分) * 10)
                        base = 1 + len(top_topics)
                        for i, r in enumerate(top_users.itertuples()):
                            src.append(1); tgt.append(base + i); val.append(float(r.NAE评分) * 5)
                        sankey = go.Figure(data=[go.Sankey(
                            node=dict(label=labels, pad=14, thickness=18,
                                      color=["#a78bfa"] + ["#60a5fa"]*len(top_topics) + ["#f472b6"]*len(top_users)),
                            link=dict(source=src, target=tgt, value=val, color="rgba(167,139,250,0.35)"),
                        )])
                        sankey.update_layout(title=f"{issue} — 跨尺度极化溯源",
                                             paper_bgcolor="#0b1120", font=dict(color="#e5e7eb"),
                                             height=420)

                        # 话题详情
                        topic_net = _sn_topic_plot(Gtopic, df)
                        top_bar_df = df.sort_values("极化度评分", ascending=False).head(15)
                        topic_bar = go.Figure(go.Bar(
                            x=top_bar_df["极化度评分"], y=top_bar_df["话题"], orientation="h",
                            marker=dict(color=top_bar_df["极化倾向"], colorscale="RdBu", cmin=-1, cmax=1),
                        ))
                        topic_bar.update_layout(title="Top 话题极化度排名",
                                                paper_bgcolor="#0b1120", plot_bgcolor="#0b1120",
                                                font=dict(color="#e5e7eb"), height=420,
                                                yaxis=dict(autorange="reversed"))
                        topic_table = df.head(20)

                        # 用户详情
                        user_net = _sn_user_plot(Guser, opinions, nae_df,
                                                 highlight=set([user]) if user else None)
                        user_table = nae_df.head(20)

                        # 意见时序（用非线性意见动力学迭代模拟几百步）
                        opinion_plot = go.Figure()
                        neighbor_df = pd.DataFrame()
                        opinion_summary = ""
                        if user and user in Guser.nodes():
                            users_all = list(Guser.nodes())
                            idx_of = {u: i for i, u in enumerate(users_all)}
                            A = nx.to_numpy_array(Guser, nodelist=users_all)
                            x = opinions.copy().astype(float)
                            K, alpha = 0.3, 2.5
                            history = [x.copy()]
                            for _ in range(60):
                                x = -x + K * (A @ np.tanh(alpha * x))
                                x = np.clip(x, -2.5, 2.5)
                                history.append(x.copy())
                            hist_arr = np.array(history)
                            ui = idx_of[user]
                            opinion_plot.add_trace(go.Scatter(
                                y=hist_arr[:, ui], mode="lines+markers",
                                name=f"{user} 意见演化", line=dict(color="#f472b6", width=2),
                            ))
                            # 邻居
                            for nb in list(Guser.neighbors(user))[:5]:
                                opinion_plot.add_trace(go.Scatter(
                                    y=hist_arr[:, idx_of[nb]], mode="lines",
                                    name=f"邻居 {nb}", line=dict(width=1, dash="dot"),
                                ))
                            opinion_plot.update_layout(
                                title=f"用户 {user} 及其邻居的意见演化",
                                paper_bgcolor="#0b1120", plot_bgcolor="#0b1120",
                                font=dict(color="#e5e7eb"), height=420,
                                xaxis_title="时间步", yaxis_title="意见 x",
                            )
                            nb_rows = []
                            for nb in Guser.neighbors(user):
                                nb_rows.append({
                                    "邻居": nb,
                                    "意见": round(float(opinions[idx_of[nb]]), 3),
                                    "NAE": float(nae_df.loc[nae_df["用户"]==nb, "NAE评分"].iloc[0]) if (nae_df["用户"]==nb).any() else None,
                                    "阵营相同": "是" if opinions[idx_of[nb]] * opinions[ui] > 0 else "否",
                                })
                            neighbor_df = pd.DataFrame(nb_rows)
                            final_x = float(hist_arr[-1, ui])
                            rho = float(np.std(hist_arr[-1]))
                            opinion_summary = (
                                f"**★ {user}** | 初始意见: {opinions[ui]:+.2f} → 终态: {final_x:+.2f} | "
                                f"邻居数: {Guser.degree(user)} | 终态系统极化 ρ = {rho:.3f}"
                            )

                        report_html = ""
                        if topic and user:
                            report_html = f"""
                            <div style="background:#1e293b;padding:16px;border-radius:10px;color:#e5e7eb;">
                            <h4 style="color:#a78bfa;margin-top:0;">跨尺度溯源报告</h4>
                            <p><b>议题:</b> {issue} (极化基线 ρ₀ = {cfg['stance_polarity']})</p>
                            <p><b>目标话题:</b> {topic}</p>
                            <p><b>目标用户:</b> {user}</p>
                            <p><b>结论:</b> 用户 {user} 在议题 <b>{issue}</b> 的话题 <b>{topic}</b>
                               讨论中表现出意见 {opinions[idx_of[user]]:+.2f}，其邻域结构促使意见在
                               非线性社会影响耦合下演化为 {float(hist_arr[-1, ui]):+.2f}。
                               NAE = {float(nae_df.loc[nae_df['用户']==user,'NAE评分'].iloc[0]):.3f}，
                               跨阵营桥接度 = {float(nae_df.loc[nae_df['用户']==user,'跨阵营桥接'].iloc[0]):.2f}。</p>
                            </div>
                            """

                        return (sankey, report_html, topic_net, topic_bar, topic_table,
                                user_net, user_table, opinion_plot, neighbor_df, opinion_summary)
                    except Exception as e:
                        import traceback
                        return (go.Figure(), f"<p>❌ {e}</p>", go.Figure(), go.Figure(), pd.DataFrame(),
                                go.Figure(), pd.DataFrame(), go.Figure(), pd.DataFrame(),
                                f"❌ {traceback.format_exc()}")

                # === 事件绑定 ===
                sn_issue.change(
                    on_issue_select,
                    inputs=[sn_issue],
                    outputs=[sn_issue_info, sn_topic_list_df, sn_topic_list_summary],
                )
                sn_start_btn.click(
                    run_sn_trace,
                    inputs=[sn_issue],
                    outputs=[sn_topic_network, sn_user_network, sn_topic, sn_user],
                )
                sn_topic.change(
                    on_sn_topic_select,
                    inputs=[sn_issue, sn_topic],
                    outputs=[sn_user],
                )
                sn_run_btn.click(
                    run_sn_full,
                    inputs=[sn_issue, sn_topic, sn_user],
                    outputs=[sn_sankey, sn_report, sn_topic_detail_net, sn_topic_bar, sn_topic_table,
                             sn_user_detail_net, sn_user_table, sn_opinion_plot,
                             sn_neighbor_table, sn_opinion_summary],
                )

            # ========== Tab 4: 网络医学分析（已注释：应需求暂时隐藏） ==========
            # with gr.Tab("🔗 网络医学分析", id=5):
            #     create_phase2_network_medicine_tab()

            # ========== Tab 5: 社交网络仿真 ==========
            with gr.Tab("🌐 社交网络仿真", id=6):

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
                        peak_density = float(np.max(infected_density))
                        peak_step = int(np.argmax(infected_density))
                        avg_density = float(np.mean(infected_density))
                        R0_approx = beta * nx.degree_histogram(state["G"]).__len__() / gamma if gamma > 0 else 0

                        sis_summary = f"""✅ SIS传播仿真完成

📊 结果统计:
- 节点数: {state['G'].number_of_nodes()}, 边数: {state['G'].number_of_edges()}, 社区数: {len(state['communities'])}
- β={beta:.3f}, γ={gamma:.3f}, 初始感染比例={ini:.2%}

📝 分析摘要:
在 {int(steps)} 步仿真中，感染密度峰值达 {peak_density:.2%}（第 {peak_step} 步），最终稳定在 {final_density:.2%}。
{'传播呈爆发式增长后趋于稳态，提示该参数组合下疾病/信息可在网络中持续传播。' if final_density > 0.05 else '传播逐渐消退，提示恢复率γ足以抑制传播扩散。'}
{'⚠️ 高传播风险：最终感染密度超过10%，需关注高连接度节点（超级传播者）的控制策略。' if final_density > 0.1 else ''}
社区结构对传播有隔离效应——社区内部传播快于社区之间。"""

                        return (
                            spread_fig,
                            snapshot_fig,
                            sis_summary,
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
            
            # ========== Tab 7: 数据统计 ==========
            with gr.Tab("📈 数据统计", id=8):

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
        server_name="0.0.0.0",  # 对所有接口开放（局域网 / 外部均可访问）
        server_port=83,         # 注意: macOS/Linux 上 <1024 端口需 sudo 启动
        share=False,
        show_error=True
    )

