"""
真实数据加载器
从TSV文件加载基因、疾病和通路数据
"""

import pandas as pd
from typing import Dict, List, Set
from collections import defaultdict
import os


class RealDataLoader:
    """从真实数据文件加载数据"""
    
    def __init__(self, data_dir="data"):
        self.data_dir = data_dir
        self.gene_disease_df = None
        self.pathway_df = None
        self.related_pathway_df = None
        
    def load_all_data(self):
        """加载所有数据文件"""
        print("📂 加载基因-疾病数据...")
        self.load_gene_disease_data()
        
        print("📂 加载通路数据...")
        self.load_pathway_data()
        
        print("📂 加载通路关系数据...")
        self.load_related_pathway_data()
        
        print("✅ 所有数据加载完成！")
    
    def load_gene_disease_data(self):
        """加载gene_disease.tsv"""
        file_path = os.path.join(self.data_dir, "gene_disease.tsv")
        if not os.path.exists(file_path):
            print(f"⚠️  文件不存在: {file_path}")
            return
        
        self.gene_disease_df = pd.read_csv(file_path, sep='\t')
        print(f"   - 加载了 {len(self.gene_disease_df)} 条记录")
    
    def load_pathway_data(self):
        """加载pathway数据"""
        file_path = os.path.join(self.data_dir, "pathway(基因名映射版).tsv")
        if not os.path.exists(file_path):
            print(f"⚠️  文件不存在: {file_path}")
            return
        
        self.pathway_df = pd.read_csv(file_path, sep='\t')
        print(f"   - 加载了 {len(self.pathway_df)} 条通路记录")
    
    def load_related_pathway_data(self):
        """加载Related Pathway.txt"""
        file_path = os.path.join(self.data_dir, "Related Pathway.txt")
        if not os.path.exists(file_path):
            print(f"⚠️  文件不存在: {file_path}")
            return
        
        self.related_pathway_df = pd.read_csv(file_path, sep='\t')
        print(f"   - 加载了 {len(self.related_pathway_df)} 条通路关系")
    
    def get_all_diseases(self) -> List[str]:
        """获取所有疾病名称"""
        if self.gene_disease_df is None:
            return []
        
        diseases = self.gene_disease_df['disease_name'].unique().tolist()
        return sorted(diseases)
    
    def get_all_pathways(self) -> List[str]:
        """获取所有通路名称"""
        if self.pathway_df is None:
            return []
        
        pathways = self.pathway_df['Pathway_Name'].unique().tolist()
        return sorted([p for p in pathways if pd.notna(p)])
    
    def get_disease_genes(self, disease_name: str) -> List[str]:
        """获取疾病相关的基因"""
        if self.gene_disease_df is None:
            return []
        
        disease_data = self.gene_disease_df[
            self.gene_disease_df['disease_name'] == disease_name
        ]
        
        genes = disease_data['gene_symbol'].unique().tolist()
        return [g for g in genes if pd.notna(g)]
    
    def get_disease_drugs(self, disease_name: str) -> List[str]:
        """获取疾病相关的药物"""
        if self.gene_disease_df is None:
            return []
        
        disease_data = self.gene_disease_df[
            self.gene_disease_df['disease_name'] == disease_name
        ]
        
        # 获取第一行的药物数据（因为同一疾病药物信息相同）
        if len(disease_data) > 0:
            drug_str = disease_data.iloc[0]['disease_drug']
            if pd.notna(drug_str):
                # 解析药物字符串
                drugs = []
                for item in drug_str.split(';'):
                    item = item.strip()
                    if '[DR:' in item:
                        # 提取药物名称（去除DR编号）
                        drug_name = item.split('[DR:')[0].strip()
                        if drug_name:
                            drugs.append(drug_name)
                return drugs[:20]  # 限制返回前20个药物
        
        return []
    
    def get_pathway_genes(self, pathway_name: str) -> List[str]:
        """获取通路包含的基因"""
        if self.pathway_df is None:
            return []
        
        pathway_data = self.pathway_df[
            self.pathway_df['Pathway_Name'] == pathway_name
        ]
        
        if len(pathway_data) > 0:
            gene_str = pathway_data.iloc[0]['Gene']
            if pd.notna(gene_str):
                # 基因用逗号或分号分隔
                genes = [g.strip() for g in gene_str.replace(';', ',').split(',')]
                return [g for g in genes if g and g != 'NA']
        
        return []
    
    def get_gene_pathways(self, gene_name: str) -> List[str]:
        """获取基因所在的通路"""
        if self.gene_disease_df is None:
            return []
        
        gene_data = self.gene_disease_df[
            self.gene_disease_df['gene_symbol'] == gene_name
        ]
        
        if len(gene_data) > 0:
            pathway_str = gene_data.iloc[0]['gene_pathway']
            if pd.notna(pathway_str):
                # 解析通路字符串
                pathways = [p.strip() for p in pathway_str.split(';')]
                # 去重并过滤空值
                pathways = list(set([p for p in pathways if p and p != 'NA']))
                return sorted(pathways)[:15]  # 限制返回前15个
        
        return []
    
    def get_related_pathways(self, pathway_id: str) -> List[str]:
        """获取相关通路"""
        if self.related_pathway_df is None:
            return []
        
        related = self.related_pathway_df[
            self.related_pathway_df['Pathway'] == pathway_id
        ]['Related Pathway'].tolist()
        
        return related
    
    def build_gene_network_from_pathway(self, genes: List[str]) -> Dict:
        """从通路基因构建网络（生成互作关系）"""
        # 简化版：创建部分基因之间的连接
        # 实际应用中应该从PPI数据库获取真实的互作关系
        
        interactions = []
        regulations = []
        
        # 为每个基因创建2-3个连接
        for i, gene in enumerate(genes):
            # 与后续2-3个基因创建连接
            for j in range(i+1, min(i+4, len(genes))):
                if hash(gene + genes[j]) % 3 > 0:  # 随机性
                    interactions.append({
                        "FirstGene": gene,
                        "SecondGene": genes[j]
                    })
                if hash(gene + genes[j]) % 2 == 0:  # 随机性
                    regulations.append({
                        "FirstGene": gene,
                        "SecondGene": genes[j]
                    })
        
        return {
            "interactions": interactions,
            "regulations": regulations
        }


# 测试代码
if __name__ == "__main__":
    loader = RealDataLoader()
    loader.load_all_data()
    
    print("\n📊 数据统计:")
    print(f"疾病数量: {len(loader.get_all_diseases())}")
    print(f"通路数量: {len(loader.get_all_pathways())}")
    
    # 测试获取一个疾病的数据
    diseases = loader.get_all_diseases()
    if diseases:
        test_disease = diseases[0]
        print(f"\n测试疾病: {test_disease}")
        genes = loader.get_disease_genes(test_disease)
        drugs = loader.get_disease_drugs(test_disease)
        print(f"  - 基因数: {len(genes)}")
        print(f"  - 药物数: {len(drugs)}")
        if genes:
            print(f"  - 示例基因: {genes[:5]}")
        if drugs:
            print(f"  - 示例药物: {drugs[:3]}")

