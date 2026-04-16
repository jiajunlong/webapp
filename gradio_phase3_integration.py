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
        """Load gene expression - prefer gene symbol indexed file"""
        try:
            # Prefer gene-symbol indexed (same fix as Phase 2)
            for expr_file in ["TCGA-COAD/filtered_hiseq_data.csv", "data/TCGA-COAD/filtered_hiseq_data.csv"]:
                if os.path.exists(expr_file):
                    self.expr_data = pd.read_csv(expr_file, index_col=0)
                    if not str(self.expr_data.index[0]).startswith('ENSG'):
                        logger.info(f"✓ Loaded expression: {self.expr_data.shape}")
                        break
                    else:
                        self.expr_data = None
        except Exception as e:
            logger.error(f"Error loading expression: {e}")

    def _load_clinical_data(self):
        """Load clinical metadata"""
        try:
            for clinical_file in ["data/TCGA-COAD/filtered_clinical.csv", "TCGA-COAD/clinical.tsv"]:
                if os.path.exists(clinical_file):
                    sep = '\t' if clinical_file.endswith('.tsv') else ','
                    self.sample_metadata = pd.read_csv(clinical_file, index_col=0, sep=sep)
                    logger.info(f"✓ Loaded clinical data: {self.sample_metadata.shape}")
                    break
        except Exception as e:
            logger.error(f"Error loading clinical data: {e}")
    
    def _load_disease_modules(self):
        """Build disease modules from gene_disease.tsv (no pkl dependency)"""
        try:
            gd_file = "data/gene_disease.tsv"
            if not os.path.exists(gd_file):
                self.disease_modules = {}
                return

            gd = pd.read_csv(gd_file, sep='\t')
            # Build modules: disease -> list of gene symbols
            self.disease_modules = {}
            for disease_name, group in gd.groupby('disease_name'):
                raw = group['gene_symbol'].dropna().tolist()
                genes = list(set(
                    str(g).split(',')[0].strip()
                    for g in raw if str(g).strip() and str(g).strip() != 'NA'
                ))
                # Only include diseases with enough genes for SIS analysis
                if len(genes) >= 10 and self.expr_data is not None:
                    available = [g for g in genes if g in self.expr_data.index]
                    if len(available) >= 10:
                        self.disease_modules[disease_name] = available

            logger.info(f"✓ Built {len(self.disease_modules)} disease modules (≥10 genes in TCGA)")
        except Exception as e:
            logger.error(f"Error building disease modules: {e}")
            self.disease_modules = {}
    
    def _load_ppi_network(self):
        """Build PPI-like network from expression correlation (no pkl dependency)"""
        # Will be built on-demand per disease module in the callback
        self.ppi_network = None
        logger.info("  PPI network: will build on-demand from expression correlation")
    
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
                    return pd.DataFrame(), None, pd.DataFrame(), pd.DataFrame(), "无法加载数据", "", "❌ 数据未加载或未选择模块", {}

                progress(0.05, desc="初始化...")

                # Get disease module genes
                disease_module_genes = data_loader.disease_modules.get(module_name, [])
                if len(disease_module_genes) < 5:
                    return (pd.DataFrame(), None, pd.DataFrame(), pd.DataFrame(),
                            "模块基因不足", "", f"❌ 模块 '{module_name}' 仅有 {len(disease_module_genes)} 个基因", {})

                # ==================== Step 1: Parameter Extraction ====================
                progress(0.15, desc="提取SIS参数...")

                extractor = ParameterExtractor(
                    disease_modules={module_name: disease_module_genes},
                    expression_data=data_loader.expr_data,
                    clinical_data=data_loader.sample_metadata,
                )

                params = extractor.extract_all_parameters(module_name)

                # Override with user-selected parameters
                params['beta'] = beta
                params['gamma'] = gamma

                # ==================== Step 2: SIS Network Propagation ====================
                progress(0.35, desc="运行SIS动力学...")

                propagation = SISNetworkPropagation(
                    adjacency_matrix=params['adjacency'],
                    gene_names=disease_module_genes,
                    parameters={
                        'beta': params['beta'],
                        'gamma': params['gamma'],
                        'initial_infection': params['initial_infection'],
                    }
                )

                propagation.run_dynamics(n_steps=n_steps, n_runs=n_runs)

                # Get biomarker results
                top_n = max(1, len(disease_module_genes) // 4)
                biomarker_df = propagation.get_biomarker_table(top_n=top_n)

                progress(0.55, desc="准备生物标志物...")

                # Format for display
                display_cols = {}
                display_cols['基因'] = biomarker_df['gene'] if 'gene' in biomarker_df.columns else biomarker_df.index
                for col in ['persistence', 'degree', 'initial_infection', 'weighted_biomarker_score']:
                    if col in biomarker_df.columns:
                        display_cols[col] = np.round(biomarker_df[col], 3)
                display_biomarker_df = pd.DataFrame(display_cols)

                biomarkers = list(display_cols['基因'])

                # ==================== Step 3: Validation ====================
                expr_validation_df = pd.DataFrame()
                clinical_validation_df = pd.DataFrame()
                lit_text = "未进行文献对标"

                if validate and len(biomarkers) > 0:
                    progress(0.65, desc="表达验证...")
                    try:
                        validator = BiomarkerValidator(
                            expression_data=data_loader.expr_data,
                            sample_metadata=data_loader.sample_metadata,
                            biomarkers=biomarkers
                        )

                        stage_col = 'Age_Group' if 'Age_Group' in data_loader.sample_metadata.columns else 'Stage'
                        ev = validator.validate_expression_changes(disease_stage_col=stage_col)
                        if ev is not None and len(ev) > 0:
                            expr_validation_df = pd.DataFrame({
                                '基因': ev.index if hasattr(ev, 'index') else ev.get('gene', []),
                                '健康平均表达': np.round(ev.get('mean_expr_healthy', ev.iloc[:, 0] if len(ev.columns) > 0 else 0), 2),
                                '疾病平均表达': np.round(ev.get('mean_expr_disease', ev.iloc[:, 1] if len(ev.columns) > 1 else 0), 2),
                                'log2倍数变化': np.round(ev.get('fold_change', 0), 2),
                                'Dysregulated': ev.get('is_dysregulated', False),
                            })
                    except Exception as ve:
                        logger.warning(f"Expression validation failed: {ve}")

                    progress(0.75, desc="临床相关性...")
                    try:
                        outcome_var = None
                        for col in data_loader.sample_metadata.columns:
                            if pd.api.types.is_numeric_dtype(data_loader.sample_metadata[col]):
                                outcome_var = col
                                break
                        if outcome_var:
                            cv = validator.validate_clinical_correlation(outcome_variable=outcome_var)
                            if cv is not None and len(cv) > 0:
                                clinical_validation_df = pd.DataFrame({
                                    '基因': cv.get('biomarker', cv.index if hasattr(cv, 'index') else []),
                                    '皮尔逊相关': np.round(cv.get('pearson_corr', 0), 3),
                                    'p值': cv.get('pearson_pval', 1).apply(lambda x: f"{x:.2e}") if 'pearson_pval' in cv else 'N/A',
                                    'Spearman相关': np.round(cv.get('spearman_corr', 0), 3),
                                    '显著': cv.get('significant', False),
                                })
                    except Exception as ce:
                        logger.warning(f"Clinical validation failed: {ce}")

                    n_dys = expr_validation_df['Dysregulated'].sum() if 'Dysregulated' in expr_validation_df.columns else 0
                    lit_text = f"已识别生物标志物: {len(biomarkers)}\n其中dysregulated: {n_dys}/{len(biomarkers)}\n临床相关: {len(clinical_validation_df)} 个"

                progress(0.92, desc="生成参数摘要...")

                param_summary = f"""SIS模型参数摘要
===============
传播速率 (β): {params['beta']:.4f}
恢复速率 (γ): {params['gamma']:.4f}

模块信息
-------
疾病模块: {module_name}
基因数: {len(disease_module_genes)}
网络边数: {params.get('n_edges', 0)}

仿真参数
-------
时间步数: {n_steps}
随机实现: {n_runs}

结果统计
-------
已识别生物标志物: {len(biomarkers)}
"""

                progress(0.96, desc="生成可视化...")
                # Dynamics plot using Plotly
                import plotly.graph_objects as go
                fig = go.Figure()
                if propagation.infection_history is not None:
                    mean_infected = propagation.infection_history.mean(axis=(0, 2))  # avg over runs and nodes
                    steps = np.arange(len(mean_infected))
                    fig.add_trace(go.Scatter(x=steps, y=mean_infected,
                        mode='lines', name='平均感染比例',
                        line=dict(color='#f5576c', width=2)))
                    fig.update_layout(
                        title=dict(text=f'SIS 感染动力学 — {module_name}', x=0.5),
                        xaxis_title='时间步', yaxis_title='感染基因比例',
                        plot_bgcolor='#fafafa', paper_bgcolor='white', height=400)

                progress(0.98, desc="完成")
                status_msg = f"""✅ SIS生物标志物发现完成

📊 结果摘要:
- 识别生物标志物: {len(biomarkers)}
- 模块大小: {len(disease_module_genes)} 基因
- β={params['beta']:.3f}, γ={params['gamma']:.3f}"""

                return (
                    display_biomarker_df, fig,
                    expr_validation_df, clinical_validation_df,
                    lit_text, param_summary,
                    status_msg, {
                        'biomarkers': biomarkers,
                        'parameters': params,
                    }
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
