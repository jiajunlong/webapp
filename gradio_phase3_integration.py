"""
Gradio App Integration Module for Phase 3: SIS-as-Biomarker-Discovery

This module provides UI components and functions for integrating Phase 3 
(Cross-Scale Parameter Propagation and SIS Biomarker Discovery) into the 
main Gradio app.

Functions:
- create_phase3_biomarker_tab(): Creates the SIS biomarker discovery subtab
- load_phase3_resources(): Initializes Phase 3 data on app startup
- run_sis_biomarker_pipeline(): Executes full Phase 3 pipeline

Status: Production Ready
Last Updated: 2026-04-15
"""

import gradio as gr
import pandas as pd
import numpy as np
import logging
from pathlib import Path
import os
import networkx as nx

# Import Phase 1-2 modules
from pathway_activity import PathwayActivityScorer
from wgcna_analysis import WGCNAAnalyzer
from disease_module_detection import DiseaseNetworkBuilder

# Import Phase 3 modules
from parameter_extraction import ParameterExtractor
from sis_network_propagation import SISNetworkPropagation
from biomarker_validation import BiomarkerValidator

logger = logging.getLogger(__name__)


class Phase3DataLoader:
    """Loads and manages Phase 1, 2, and 3 data"""
    
    def __init__(self):
        self.expr_data = None
        self.sample_metadata = None
        self.disease_modules = None
        self.ppi_network = None
        self.mirna_data = None
        self.pathway_genes = None
        self.wgcna_modules = None
        self.loaded = False
    
    def load_all(self):
        """Load all data required for Phase 3"""
        try:
            logger.info("Loading Phase 3 data (Phase 1-2 integration)...")
            
            # Load Phase 1: Gene expression
            self._load_expression_data()
            
            # Load Phase 1: Metadata
            self._load_clinical_data()
            
            # Load Phase 2: Disease modules
            self._load_disease_modules()
            
            # Load Phase 2: PPI network
            self._load_ppi_network()
            
            # Load Phase 2: WGCNA modules (optional)
            self._load_wgcna_modules()
            
            self.loaded = True
            logger.info("✓ Phase 3 data loaded successfully")
            return True
            
        except Exception as e:
            logger.error(f"Error loading Phase 3 data: {e}")
            return False
    
    def _load_expression_data(self):
        """Load gene expression from Phase 1"""
        try:
            expr_file = "data/TCGA-COAD/filtered_hiseq_data.csv"
            if os.path.exists(expr_file):
                self.expr_data = pd.read_csv(expr_file, index_col=0)
                logger.info(f"✓ Loaded expression: {self.expr_data.shape}")
        except Exception as e:
            logger.error(f"Error loading expression: {e}")
    
    def _load_clinical_data(self):
        """Load clinical metadata"""
        try:
            clinical_file = "data/TCGA-COAD/filtered_clinical.csv"
            if os.path.exists(clinical_file):
                self.sample_metadata = pd.read_csv(clinical_file, index_col=0)
                logger.info(f"✓ Loaded clinical data: {self.sample_metadata.shape}")
        except Exception as e:
            logger.error(f"Error loading clinical data: {e}")
    
    def _load_disease_modules(self):
        """Load disease modules from Phase 2"""
        try:
            modules_file = "data/phase2_disease_modules.pkl"
            if os.path.exists(modules_file):
                import pickle
                with open(modules_file, 'rb') as f:
                    self.disease_modules = pickle.load(f)
                logger.info(f"✓ Loaded {len(self.disease_modules)} disease modules")
        except Exception as e:
            logger.error(f"Error loading disease modules: {e}")
            # Fallback: create empty dict
            self.disease_modules = {}
    
    def _load_ppi_network(self):
        """Load PPI network from Phase 2"""
        try:
            ppi_file = "data/phase2_ppi_network.pkl"
            if os.path.exists(ppi_file):
                import pickle
                with open(ppi_file, 'rb') as f:
                    self.ppi_network = pickle.load(f)
                logger.info(f"✓ Loaded PPI network: {self.ppi_network.number_of_nodes()} nodes")
        except Exception as e:
            logger.error(f"Error loading PPI network: {e}")
    
    def _load_wgcna_modules(self):
        """Load WGCNA modules from Phase 2 (optional)"""
        try:
            wgcna_file = "data/phase2_wgcna_modules.pkl"
            if os.path.exists(wgcna_file):
                import pickle
                with open(wgcna_file, 'rb') as f:
                    self.wgcna_modules = pickle.load(f)
                logger.info(f"✓ Loaded {len(self.wgcna_modules)} WGCNA modules")
        except Exception as e:
            logger.debug(f"Optional WGCNA data not available: {e}")


def create_phase3_biomarker_tab():
    """
    Creates the Phase 3 SIS Biomarker Discovery subtab for Tab 7
    
    Returns:
        None (UI components added directly to parent container)
    """
    
    # Initialize data loader
    data_loader = Phase3DataLoader()
    data_loaded = data_loader.load_all()
    
    with gr.Group():
        gr.Markdown("""
        ## 🔬 SIS网络传播生物标志物发现 (SIS Network Propagation Biomarker Discovery)
        
        基于Phase 1-2多尺度整合的创新方法:
        
        ### 分析步骤
        1. **参数提取**: 从疾病模块提取SIS模型参数 (β, γ, I₀)
        2. **网络传播**: 在蛋白质相互作用网络上运行SIS动力学
        3. **持久性评分**: 计算基因在感染态下的持久性
        4. **生物标志物识别**: 选择高持久性基因作为候选生物标志物
        5. **临床验证**: 与临床结果和文献进行对标
        
        ### 生物学原理
        - 将dysregulated基因视为"感染"状态
        - 通过PPI网络的传播模型疾病的扩散动力学
        - 持久存在的基因作为关键生物标志物
        """)
        
        with gr.Row():
            with gr.Column(scale=1):
                gr.Markdown("### ⚙️ 分析参数")
                
                # Disease module selection
                disease_module_choices = list(data_loader.disease_modules.keys()) if data_loader.disease_modules else []
                
                disease_module_sel = gr.Dropdown(
                    choices=disease_module_choices,
                    value=disease_module_choices[0] if disease_module_choices else None,
                    label="🔗 疾病模块",
                    info="选择要分析的疾病模块"
                )
                
                # SIS parameters
                beta_slider = gr.Slider(
                    minimum=0.01,
                    maximum=1.0,
                    value=0.3,
                    step=0.05,
                    label="β (传播速率)",
                    info="较高值表示疾病更具传染性"
                )
                
                gamma_slider = gr.Slider(
                    minimum=0.01,
                    maximum=1.0,
                    value=0.2,
                    step=0.05,
                    label="γ (恢复速率)",
                    info="较低值表示疾病更持久"
                )
                
                n_steps_slider = gr.Slider(
                    minimum=100,
                    maximum=2000,
                    value=500,
                    step=100,
                    label="⏱️ 仿真步数",
                    info="SIS动力学仿真的时间步数"
                )
                
                n_runs_slider = gr.Slider(
                    minimum=10,
                    maximum=200,
                    value=50,
                    step=10,
                    label="🔄 随机实现数",
                    info="随机SIS轨迹的数量"
                )
                
                biomarker_percentile = gr.Slider(
                    minimum=50,
                    maximum=99,
                    value=75,
                    step=1,
                    label="📊 百分位阈值",
                    info="选择持久性得分最高的基因百分比"
                )
                
                validate_checkbox = gr.Checkbox(
                    value=True,
                    label="🔍 进行临床验证",
                    info="运行生物标志物临床相关性验证"
                )
                
                run_button = gr.Button(
                    "▶️ 开始SIS分析",
                    variant="primary",
                    size="lg"
                )
                
                status_text = gr.Textbox(
                    label="状态",
                    value="就绪" if data_loaded else "⚠️ 数据加载失败",
                    interactive=False,
                    lines=3
                )
            
            with gr.Column(scale=3):
                gr.Markdown("### 📊 分析结果")
                
                with gr.Tabs():
                    with gr.Tab("📈 生物标志物排序"):
                        biomarker_table = gr.DataFrame(
                            label="SIS发现的生物标志物",
                            headers=["基因", "持久性", "网络度", "初始感染", "综合评分"],
                            datatype=["str", "number", "number", "number", "number"],
                            row_count=15
                        )
                    
                    with gr.Tab("📉 感染动力学"):
                        dynamics_plot = gr.Plot(label="平均感染比例随时间变化")
                    
                    with gr.Tab("🧬 表达验证"):
                        expr_validation_table = gr.DataFrame(
                            label="基因表达改变验证",
                            headers=["基因", "健康平均表达", "疾病平均表达", "log2倍数变化", "Dysregulated"],
                            datatype=["str", "number", "number", "number", "bool"],
                            row_count=12
                        )
                    
                    with gr.Tab("❤️ 临床相关性"):
                        clinical_table = gr.DataFrame(
                            label="与临床结果的相关性",
                            headers=["基因", "皮尔逊相关", "p值", "Spearman相关", "显著"],
                            datatype=["str", "number", "number", "number", "bool"],
                            row_count=12
                        )
                    
                    with gr.Tab("📚 文献对标"):
                        literature_text = gr.Textbox(
                            label="已知生物标志物对标",
                            lines=6,
                            interactive=False
                        )
                    
                    with gr.Tab("📊 参数统计"):
                        params_text = gr.Textbox(
                            label="提取的SIS参数摘要",
                            lines=8,
                            interactive=False
                        )
        
        # State for storing results
        results_state = gr.State({
            'biomarkers': None,
            'persistence_scores': None,
            'dynamics': None,
            'expr_validation': None,
            'clinical_validation': None,
            'parameters': None
        })
        
        def run_phase3_analysis(
            module_name, beta, gamma, n_steps, n_runs, 
            percentile, validate, progress=gr.Progress()
        ):
            """
            Execute full Phase 3 SIS biomarker discovery pipeline
            
            Parameters:
            -----------
            module_name : str
                Disease module name
            beta : float
                Transmission rate
            gamma : float
                Recovery rate
            n_steps : int
                Number of simulation steps
            n_runs : int
                Number of stochastic runs
            percentile : int
                Biomarker selection percentile
            validate : bool
                Whether to perform validation
            
            Returns:
            --------
            tuple : (biomarker_table, dynamics_plot, expr_table, 
                    clinical_table, literature_text, params_text, 
                    status, results_state)
            """
            
            try:
                if not data_loaded or module_name is None:
                    return None, None, None, None, "无法加载数据", {}, "❌ 数据未加载"
                
                progress(0.05, desc="初始化...")
                
                # ==================== Step 1: Parameter Extraction ====================
                progress(0.15, desc="提取SIS参数...")
                
                disease_module_genes = data_loader.disease_modules.get(module_name, [])
                if not disease_module_genes:
                    return (None, None, None, None, "无法找到模块", 
                           {}, "❌ 模块不存在")
                
                # Create parameter extractor
                extractor = ParameterExtractor(
                    expression_data=data_loader.expr_data,
                    sample_metadata=data_loader.sample_metadata,
                    disease_modules={module_name: disease_module_genes},
                    ppi_network=data_loader.ppi_network if data_loader.ppi_network else None
                )
                
                params_all = extractor.extract_all_parameters()
                params = params_all[module_name]
                
                # Override with user-selected parameters if provided
                if beta > 0:
                    params['beta'] = beta
                if gamma > 0:
                    params['gamma'] = gamma
                
                # ==================== Step 2: SIS Network Propagation ====================
                progress(0.35, desc="运行SIS动力学...")
                
                propagation = SISNetworkPropagation(
                    network_adjacency=params['adjacency'],
                    gene_names=disease_module_genes,
                    beta=params['beta'],
                    gamma=params['gamma'],
                    initial_infection=params['initial_infection'],
                    random_seed=42
                )
                
                propagation.run_dynamics(n_steps=n_steps, n_runs=n_runs)
                
                # Get biomarkers
                persistence_scores = propagation.get_persistence_scores()
                biomarker_df = propagation.get_biomarker_table(
                    top_n=max(1, len(disease_module_genes) // 4)
                )
                
                # Get dynamics
                dynamics_df = propagation.get_infection_dynamics()
                
                progress(0.55, desc="准备生物标志物...")
                
                biomarkers = biomarker_df['gene'].tolist()
                
                # Convert biomarker table for display
                display_biomarker_df = pd.DataFrame({
                    '基因': biomarker_df['gene'],
                    '持久性': np.round(biomarker_df['persistence'], 3),
                    '网络度': biomarker_df['degree'].astype(int),
                    '初始感染': np.round(biomarker_df['initial_infection'], 3),
                    '综合评分': np.round(biomarker_df['weighted_biomarker_score'], 3)
                })
                
                # ==================== Step 3: Validation (Optional) ====================
                expr_validation_df = None
                clinical_validation_df = None
                literature_text = "未进行文献对标"
                
                if validate and len(biomarkers) > 0:
                    progress(0.65, desc="表达验证...")
                    
                    validator = BiomarkerValidator(
                        expression_data=data_loader.expr_data,
                        sample_metadata=data_loader.sample_metadata,
                        biomarkers=biomarkers
                    )
                    
                    expr_validation_df = validator.validate_expression_changes(
                        disease_stage_col='Age_Group' if 'Age_Group' in data_loader.sample_metadata.columns else 'Stage'
                    )
                    
                    # Format expression validation for display
                    if expr_validation_df is not None:
                        display_expr_df = pd.DataFrame({
                            '基因': expr_validation_df.index,
                            '健康平均表达': np.round(expr_validation_df['mean_expr_healthy'], 2),
                            '疾病平均表达': np.round(expr_validation_df['mean_expr_disease'], 2),
                            'log2倍数变化': np.round(expr_validation_df['fold_change'], 2),
                            'Dysregulated': expr_validation_df['is_dysregulated']
                        })
                        expr_validation_df = display_expr_df
                    
                    progress(0.75, desc="临床相关性验证...")
                    
                    # Find numeric outcome variable
                    outcome_var = None
                    for col in data_loader.sample_metadata.columns:
                        if pd.api.types.is_numeric_dtype(data_loader.sample_metadata[col]):
                            outcome_var = col
                            break
                    
                    if outcome_var:
                        clinical_validation_df = validator.validate_clinical_correlation(
                            outcome_variable=outcome_var
                        )
                        
                        # Format clinical validation for display
                        display_clinical_df = pd.DataFrame({
                            '基因': clinical_validation_df['biomarker'],
                            '皮尔逊相关': np.round(clinical_validation_df['pearson_corr'], 3),
                            'p值': clinical_validation_df['pearson_pval'].apply(lambda x: f"{x:.2e}"),
                            'Spearman相关': np.round(clinical_validation_df['spearman_corr'], 3),
                            '显著': clinical_validation_df['significant']
                        })
                        clinical_validation_df = display_clinical_df
                    
                    progress(0.85, desc="文献对标...")
                    
                    # Try to compare with known biomarkers
                    lit_comp = validator.compare_with_literature(
                        known_biomarkers=[]  # Would need literature reference data
                    )
                    
                    literature_text = f"""
                    已识别的生物标志物: {len(biomarkers)}
                    其中dysregulated: {(expr_validation_df['Dysregulated'].sum() if expr_validation_df is not None else 0)}/{len(biomarkers)}
                    临床相关的: {(clinical_validation_df is not None and len(clinical_validation_df) or 0)} 个
                    """
                
                progress(0.92, desc="生成参数摘要...")
                
                # Generate parameter summary
                param_summary = f"""
                SIS模型参数摘要
                ===============
                
                传播速率 (β): {params['beta']:.4f}
                恢复速率 (γ): {params['gamma']:.4f}
                
                模块信息
                -------
                疾病模块: {module_name}
                基因数: {len(disease_module_genes)}
                网络边数: {params['n_edges']}
                平均连接度: {params['n_edges'] * 2 / len(disease_module_genes):.2f}
                
                仿真参数
                -------
                时间步数: {n_steps}
                随机实现: {n_runs}
                总轨迹数: {n_steps * n_runs}
                
                结果统计
                -------
                已识别生物标志物: {len(biomarkers)}
                平均持久性: {np.mean(persistence_scores):.3f}
                最高持久性: {np.max(persistence_scores):.3f}
                最低持久性: {np.min(persistence_scores):.3f}
                """
                
                progress(0.98, desc="完成...")
                
                status_msg = f"""
                ✅ SIS生物标志物发现完成
                
                📊 结果摘要:
                - 识别生物标志物: {len(biomarkers)}
                - 模块大小: {len(disease_module_genes)} 基因
                - β (传播速率): {params['beta']:.3f}
                - γ (恢复速率): {params['gamma']:.3f}
                - 平均持久性: {np.mean(persistence_scores):.3f}
                """
                
                # Create dynamics visualization
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(10, 6))
                
                if len(dynamics_df) > 0:
                    ax.plot(dynamics_df['timestep'], dynamics_df['mean_infected'], 
                           linewidth=2, color='red', label='平均感染比例')
                    ax.fill_between(dynamics_df['timestep'], 
                                   dynamics_df['mean_infected'] - dynamics_df['std_infected'],
                                   dynamics_df['mean_infected'] + dynamics_df['std_infected'],
                                   alpha=0.2, color='red')
                    ax.set_xlabel('时间步')
                    ax.set_ylabel('感染基因比例')
                    ax.set_title(f'SIS动力学: {module_name}')
                    ax.legend()
                    ax.grid(True, alpha=0.3)
                
                # Store results
                results_state = {
                    'biomarkers': biomarkers,
                    'persistence_scores': persistence_scores,
                    'dynamics': dynamics_df,
                    'expr_validation': expr_validation_df,
                    'clinical_validation': clinical_validation_df,
                    'parameters': params
                }
                
                return (
                    display_biomarker_df,
                    fig,
                    expr_validation_df,
                    clinical_validation_df,
                    literature_text,
                    param_summary,
                    status_msg,
                    results_state
                )
                
            except Exception as e:
                logger.error(f"Error in Phase 3 analysis: {e}", exc_info=True)
                error_msg = f"❌ 错误: {str(e)}"
                return None, None, None, None, "分析失败", {}, error_msg, {}
        
        # Connect button callback
        run_button.click(
            run_phase3_analysis,
            inputs=[
                disease_module_sel, beta_slider, gamma_slider, 
                n_steps_slider, n_runs_slider, biomarker_percentile, 
                validate_checkbox
            ],
            outputs=[
                biomarker_table, dynamics_plot, expr_validation_table,
                clinical_table, literature_text, params_text,
                status_text, results_state
            ]
        )
    
    return results_state


def get_phase3_integration_instructions():
    """
    Returns instructions for integrating Phase 3 into app_full.py
    """
    return """
    # Phase 3 Gradio Integration Instructions
    
    ## Location in app_full.py
    Add Phase 3 as a new Tab 7 (or extend existing Tab 5/6).
    
    ## Integration Code
    
    1. Add imports at top of app_full.py:
    ```python
    from gradio_phase3_integration import (
        create_phase3_biomarker_tab, 
        Phase3DataLoader
    )
    ```
    
    2. In main Gradio interface, add new tab:
    ```python
    with gr.Tab("🔬 Phase 3: SIS生物标志物"):
        create_phase3_biomarker_tab()
    ```
    
    ## Data Requirements
    
    The following files must exist for Phase 3 to work:
    - data/TCGA-COAD/filtered_hiseq_data.csv (Phase 1)
    - data/TCGA-COAD/filtered_clinical.csv (Phase 1)
    - data/phase2_disease_modules.pkl (Phase 2)
    - data/phase2_ppi_network.pkl (Phase 2)
    
    ## Performance Notes
    
    - Parameter extraction: ~2-5 seconds
    - SIS dynamics (500 steps × 50 runs): ~10-30 seconds
    - Validation: ~5-10 seconds
    - Total: ~20-50 seconds per analysis
    
    """


if __name__ == "__main__":
    print(get_phase3_integration_instructions())
