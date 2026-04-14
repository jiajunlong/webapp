"""
TCGA-COAD 基因网络仿真模块
基于年龄、性别、疾病阶段的基因网络构建和分析
"""

import pandas as pd
import numpy as np
from sklearn.feature_selection import mutual_info_regression
from itertools import combinations
from joblib import Parallel, delayed
from typing import Dict, List, Tuple, Optional
import os
import time


class TCGA_COAD_Simulator:
    """TCGA-COAD 基因网络仿真器"""
    
    def __init__(self, data_dir: str = "TCGA-COAD"):
        """
        初始化仿真器
        
        参数:
            data_dir: 数据文件目录
        """
        self.data_dir = data_dir
        self.mirna_data = None
        self.gene_data = None
        self.clinical_data = None
        
    def load_data(self, 
                  mirna_file: str = "filtered_miRNA_with_names.csv",
                  gene_file: str = "filtered_hiseq_data.csv",
                  clinical_file: str = "clinical.tsv"):
        """
        加载数据文件
        
        参数:
            mirna_file: miRNA表达数据文件
            gene_file: 基因表达数据文件
            clinical_file: 临床信息文件
        """
        print("📂 正在加载数据文件...")
        
        # 加载miRNA数据
        mirna_path = os.path.join(self.data_dir, mirna_file)
        if os.path.exists(mirna_path):
            self.mirna_data = pd.read_csv(mirna_path, index_col=0)
            print(f"   ✅ miRNA数据: {self.mirna_data.shape[0]} 个miRNA, {self.mirna_data.shape[1]} 个样本")
        else:
            print(f"   ⚠️  未找到文件: {mirna_path}")
        
        # 加载基因数据
        gene_path = os.path.join(self.data_dir, gene_file)
        if os.path.exists(gene_path):
            self.gene_data = pd.read_csv(gene_path, index_col=0)
            print(f"   ✅ 基因数据: {self.gene_data.shape[0]} 个基因, {self.gene_data.shape[1]} 个样本")
        else:
            print(f"   ⚠️  未找到文件: {gene_path}")
        
        # 加载临床数据
        clinical_path = os.path.join(self.data_dir, clinical_file)
        if os.path.exists(clinical_path):
            self.clinical_data = pd.read_csv(clinical_path, sep="\t")
            print(f"   ✅ 临床数据: {len(self.clinical_data)} 条记录")
        else:
            # 尝试其他格式
            clinical_path_txt = os.path.join(self.data_dir, "clinical.txt")
            if os.path.exists(clinical_path_txt):
                self.clinical_data = pd.read_csv(clinical_path_txt, sep="\t")
                print(f"   ✅ 临床数据: {len(self.clinical_data)} 条记录")
            else:
                print(f"   ⚠️  未找到临床数据文件")
        
        print("✅ 数据加载完成")
    
    def _calculate_mi(self, i: int, j: int, data: pd.DataFrame) -> float:
        """计算互信息"""
        return mutual_info_regression(data.iloc[:, [i]], data.iloc[:, j])[0]
    
    def _calculate_mrnetb(self, i: int, j: int, mi_matrix: np.ndarray) -> float:
        """计算MRNetB网络权重"""
        if i != j:
            min_mi = np.minimum(mi_matrix[i, :], mi_matrix[j, :])
            return mi_matrix[i, j] - np.max(min_mi)
        return 0
    
    def build_network_mrnetb(self, 
                            expression_data: pd.DataFrame,
                            data_type: str = "gene",
                            max_features: int = 200,
                            feature_selection: str = "variance",
                            min_weight_threshold: float = 0.0,
                            progress_callback=None) -> pd.DataFrame:
        """
        使用MRNetB算法构建网络（优化版）
        
        参数:
            expression_data: 表达数据（行为特征，列为样本）
            data_type: 数据类型 ("gene" 或 "mirna")
            max_features: 最大特征数量（用于加速计算，默认200）
            feature_selection: 特征选择方法 ("variance", "mean", "random")
            min_weight_threshold: 最小权重阈值（过滤弱连接，默认0.0）
            progress_callback: 进度回调函数
        
        返回:
            边列表DataFrame (gene1/gene2, weight)
        """
        if expression_data is None or expression_data.empty:
            return pd.DataFrame()
        
        # 处理缺失值
        expression_data = expression_data.fillna(0)
        
        # 删除全为0的特征
        zero_features = expression_data.sum(axis=1) == 0
        expression_data_filtered = expression_data[~zero_features]
        
        if len(expression_data_filtered) == 0:
            return pd.DataFrame()
        
        # 特征筛选：如果特征太多，进行筛选
        n_total = len(expression_data_filtered)
        if n_total > max_features:
            if progress_callback:
                progress_callback(0.05, f"特征筛选: {n_total} -> {max_features}...")
            
            if feature_selection == "variance":
                # 按方差选择（选择变化最大的特征）
                variances = expression_data_filtered.var(axis=1)
                top_features = variances.nlargest(max_features).index
                expression_data_filtered = expression_data_filtered.loc[top_features]
            elif feature_selection == "mean":
                # 按平均表达量选择
                means = expression_data_filtered.mean(axis=1)
                top_features = means.nlargest(max_features).index
                expression_data_filtered = expression_data_filtered.loc[top_features]
            else:  # random
                # 随机采样
                expression_data_filtered = expression_data_filtered.sample(n=max_features, random_state=42)
        
        # 转置数据（样本为行，特征为列）
        data_t = expression_data_filtered.T
        
        n_features = data_t.shape[1]
        
        if progress_callback:
            progress_callback(0.1, "计算互信息矩阵...")
        
        # 计算互信息矩阵
        mi_matrix = np.zeros((n_features, n_features))
        results = Parallel(n_jobs=-1)(
            delayed(self._calculate_mi)(i, j, data_t)
            for i, j in combinations(range(n_features), 2)
        )
        
        for idx, (i, j) in enumerate(combinations(range(n_features), 2)):
            mi_matrix[i, j] = results[idx]
            mi_matrix[j, i] = results[idx]
        
        if progress_callback:
            progress_callback(0.6, "计算MRNetB网络...")
        
        # 优化：直接计算MRNetB网络（避免重复计算）
        mrnetb_network = np.zeros((n_features, n_features))
        
        # 只计算上三角矩阵（对称矩阵）
        for i in range(n_features):
            for j in range(i + 1, n_features):
                if i != j:
                    min_mi = np.minimum(mi_matrix[i, :], mi_matrix[j, :])
                    weight = mi_matrix[i, j] - np.max(min_mi)
                    # 只保留大于阈值的权重
                    if weight > min_weight_threshold:
                        mrnetb_network[i, j] = weight
                        mrnetb_network[j, i] = weight  # 对称矩阵
        
        if progress_callback:
            progress_callback(0.9, "生成边列表...")
        
        # 转换为边列表（只保留非零边）
        edges = np.where(mrnetb_network > min_weight_threshold)
        
        if data_type == "gene":
            edge_list = pd.DataFrame({
                "gene1": expression_data_filtered.index[edges[0]],
                "gene2": expression_data_filtered.index[edges[1]],
                "weight": mrnetb_network[edges]
            })
        else:  # mirna
            edge_list = pd.DataFrame({
                "mirna1": expression_data_filtered.index[edges[0]],
                "mirna2": expression_data_filtered.index[edges[1]],
                "weight": mrnetb_network[edges]
            })
        
        if progress_callback:
            progress_callback(1.0, "完成")
        
        return edge_list
    
    def analyze_by_age(self, 
                      age_groups: Dict[str, Tuple[float, float]] = None,
                      data_type: str = "gene",
                      max_features: int = 200,
                      feature_selection: str = "variance",
                      min_weight_threshold: float = 0.0,
                      progress_callback=None) -> Dict[str, pd.DataFrame]:
        """
        按年龄分组分析
        
        参数:
            age_groups: 年龄分组字典，例如 {"young": (0, 50), "middle": (50, 70), "old": (70, 200)}
            data_type: 数据类型 ("gene" 或 "mirna")
            progress_callback: 进度回调函数
        
        返回:
            每个年龄组的网络边列表字典
        """
        if age_groups is None:
            age_groups = {
                "young": (0, 50),
                "middle": (50, 70),
                "old": (70, 200)
            }
        
        if self.clinical_data is None or self.gene_data is None:
            return {}
        
        # 选择数据
        data = self.gene_data if data_type == "gene" else self.mirna_data
        if data is None:
            return {}
        
        # 提取肿瘤样本
        tumor_samples = [col for col in data.columns if col.endswith("01")]
        data_tumor = data[tumor_samples]
        
        # 处理年龄信息
        if "age_at_index" in self.clinical_data.columns:
            self.clinical_data["age_at_index"] = pd.to_numeric(
                self.clinical_data["age_at_index"], errors="coerce"
            )
            age_info = self.clinical_data.set_index("case_submitter_id")["age_at_index"].to_dict()
        else:
            return {}
        
        results = {}
        total_groups = len(age_groups)
        
        for idx, (group_name, (min_age, max_age)) in enumerate(age_groups.items()):
            if progress_callback:
                progress_callback(
                    idx / total_groups,
                    f"处理 {group_name} 组 ({min_age}-{max_age}岁)..."
                )
            
            # 筛选样本
            if max_age == 200:  # old组
                group_samples = [
                    col for col in data_tumor.columns
                    if age_info.get(col[:12], np.nan) > min_age
                ]
            else:
                group_samples = [
                    col for col in data_tumor.columns
                    if min_age <= age_info.get(col[:12], np.nan) < max_age
                ]
            
            if len(group_samples) == 0:
                continue
            
            group_data = data_tumor[group_samples]
            
            # 构建网络
            network = self.build_network_mrnetb(
                group_data,
                data_type=data_type,
                max_features=max_features,
                feature_selection=feature_selection,
                min_weight_threshold=min_weight_threshold,
                progress_callback=None
            )
            
            results[group_name] = network
        
        return results
    
    def analyze_by_gender(self,
                         data_type: str = "gene",
                         max_features: int = 200,
                         feature_selection: str = "variance",
                         min_weight_threshold: float = 0.0,
                         progress_callback=None) -> Dict[str, pd.DataFrame]:
        """
        按性别分组分析
        
        参数:
            data_type: 数据类型 ("gene" 或 "mirna")
            progress_callback: 进度回调函数
        
        返回:
            每个性别组的网络边列表字典
        """
        if self.clinical_data is None:
            return {}
        
        # 选择数据
        data = self.gene_data if data_type == "gene" else self.mirna_data
        if data is None:
            return {}
        
        # 提取肿瘤样本
        tumor_samples = [col for col in data.columns if col.endswith("01")]
        data_tumor = data[tumor_samples]
        
        # 获取性别信息
        if "gender" not in self.clinical_data.columns:
            return {}
        
        gender_info = self.clinical_data.set_index("case_submitter_id")["gender"].to_dict()
        
        results = {}
        genders = ["male", "female"]
        
        for idx, gender in enumerate(genders):
            if progress_callback:
                progress_callback(
                    idx / len(genders),
                    f"处理 {gender} 组..."
                )
            
            # 筛选样本
            gender_samples = [
                col for col in data_tumor.columns
                if gender_info.get(col[:12]) == gender
            ]
            
            if len(gender_samples) == 0:
                continue
            
            gender_data = data_tumor[gender_samples]
            
            # 构建网络
            network = self.build_network_mrnetb(
                gender_data,
                data_type=data_type,
                max_features=max_features,
                feature_selection=feature_selection,
                min_weight_threshold=min_weight_threshold,
                progress_callback=None
            )
            
            results[gender] = network
        
        return results
    
    def analyze_by_stage(self,
                        stage_type: str = "single",
                        data_type: str = "gene",
                        max_features: int = 200,
                        feature_selection: str = "variance",
                        min_weight_threshold: float = 0.0,
                        progress_callback=None) -> Dict[str, pd.DataFrame]:
        """
        按疾病阶段分组分析
        
        参数:
            stage_type: 阶段类型 ("single" 或 "combined")
            data_type: 数据类型 ("gene" 或 "mirna")
            progress_callback: 进度回调函数
        
        返回:
            每个阶段的网络边列表字典
        """
        if self.clinical_data is None:
            return {}
        
        # 选择数据
        data = self.gene_data if data_type == "gene" else self.mirna_data
        if data is None:
            return {}
        
        # 阶段映射
        stage_mapping = {
            "Stage I": ["Stage I", "Stage IA", "Stage IB"],
            "Stage II": ["Stage II", "Stage IIA", "Stage IIB", "Stage IIC"],
            "Stage III": ["Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"],
            "Stage IV": ["Stage IV", "Stage IVA", "Stage IVB"]
        }
        
        if "ajcc_pathologic_stage" not in self.clinical_data.columns:
            return {}
        
        results = {}
        stages = list(stage_mapping.keys())
        
        # 提取正常样本（如果combined模式）
        normal_samples = []
        if stage_type == "combined":
            normal_samples = [col for col in data.columns if col.endswith("11")]
        
        # 提取肿瘤样本
        tumor_samples = [col for col in data.columns if col.endswith("01")]
        data_tumor = data[tumor_samples]
        
        for idx, stage in enumerate(stages):
            if progress_callback:
                progress_callback(
                    idx / len(stages),
                    f"处理 {stage}..."
                )
            
            # 筛选该阶段的样本
            selected_stages = stage_mapping[stage]
            selected_samples = self.clinical_data[
                self.clinical_data["ajcc_pathologic_stage"].isin(selected_stages)
            ]["case_submitter_id"]
            
            stage_samples = [
                col for col in data_tumor.columns
                if col[:12] in selected_samples.values
            ]
            
            if len(stage_samples) == 0:
                continue
            
            # 合并数据（如果是combined模式）
            if stage_type == "combined" and len(normal_samples) > 0:
                data_normal = data[normal_samples]
                stage_data = pd.concat([data_normal, data[stage_samples]], axis=1)
            else:
                stage_data = data[stage_samples]
            
            # 构建网络
            network = self.build_network_mrnetb(
                stage_data,
                data_type=data_type,
                max_features=max_features,
                feature_selection=feature_selection,
                min_weight_threshold=min_weight_threshold,
                progress_callback=None
            )
            
            results[stage] = network
        
        return results
    
    def calculate_mfe(self,
                     expression_values: np.ndarray,
                     adj_matrix: np.ndarray) -> np.ndarray:
        """
        计算网络流熵 (MFE - Network Flow Entropy)
        
        参数:
            expression_values: 表达值数组
            adj_matrix: 邻接矩阵
        
        返回:
            MFE值数组
        """
        # 归一化邻接矩阵
        row_sums = np.sum(adj_matrix, axis=1, keepdim=True)
        row_sums[row_sums == 0] = 1
        p_ij = adj_matrix / row_sums
        
        # 归一化表达值
        pi = expression_values / np.sum(expression_values)
        
        # 计算流
        pi_expanded = np.expand_dims(pi, axis=1)
        flow = pi_expanded * p_ij
        
        # 计算MFE
        flow_log = np.log(flow + 1e-10)
        mfe = -np.sum(flow * flow_log, axis=1)
        
        return mfe
    
    def network_to_adj_matrix(self, network_df: pd.DataFrame, 
                             feature_names: List[str],
                             data_type: str = "gene") -> np.ndarray:
        """
        将网络边列表转换为邻接矩阵
        
        参数:
            network_df: 网络边列表
            feature_names: 特征名称列表
            data_type: 数据类型
        
        返回:
            邻接矩阵
        """
        n = len(feature_names)
        adj_matrix = np.zeros((n, n))
        feature_to_idx = {name: idx for idx, name in enumerate(feature_names)}
        
        if data_type == "gene":
            col1, col2 = "gene1", "gene2"
        else:
            col1, col2 = "mirna1", "mirna2"
        
        for _, row in network_df.iterrows():
            if row[col1] in feature_to_idx and row[col2] in feature_to_idx:
                i = feature_to_idx[row[col1]]
                j = feature_to_idx[row[col2]]
                weight = row["weight"]
                adj_matrix[i, j] = weight
                adj_matrix[j, i] = weight
        
        # 归一化
        max_val = np.max(np.abs(adj_matrix))
        if max_val > 0:
            adj_matrix = adj_matrix / max_val
        
        return adj_matrix

