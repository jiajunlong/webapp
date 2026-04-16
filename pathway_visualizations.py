"""
Pathway Visualization Module - Phase 1 Module 4

Creates interactive visualizations for pathway analysis results:
- Heatmaps: Pathway activity across patient groups
- Violin plots: Distribution of pathway activity by clinical strata
- Bar charts: Hub gene rankings within pathways
- Network plots: Pathway-gene networks with hub genes highlighted
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def plot_pathway_activity_heatmap(pathway_activity: pd.DataFrame,
                                  clinical_data: pd.DataFrame,
                                  group_by: str,
                                  figsize: tuple = (1200, 600)) -> go.Figure:
    """
    Create heatmap of pathway activity across patient groups.
    
    Rows: Pathways (347 KEGG pathways)
    Columns: Patient samples (ordered by clinical variable)
    Color: Z-score normalized pathway activity (red=high, blue=low)
    
    Parameters:
    -----------
    pathway_activity : pd.DataFrame
        Pathways × samples matrix (347 × 257)
    clinical_data : pd.DataFrame
        Sample × clinical variables (must have index matching pathway_activity columns)
    group_by : str
        Clinical variable to group samples by (e.g., 'gender', 'ajcc_pathologic_stage')
    figsize : tuple
        Figure size (width, height)
        
    Returns:
    --------
    go.Figure : Interactive Plotly heatmap
    """
    logger.info(f"Creating heatmap: pathway activity grouped by '{group_by}'")
    
    # Align indices
    common_samples = list(set(pathway_activity.columns) & set(clinical_data.index))
    pathway_activity = pathway_activity[common_samples]
    clinical_data = clinical_data.loc[common_samples]
    
    # Sort by grouping variable
    sort_idx = clinical_data[group_by].argsort()
    pathway_activity = pathway_activity.iloc[:, sort_idx]
    clinical_data = clinical_data.iloc[sort_idx]
    
    # Z-score normalize each pathway
    pathway_activity_norm = (pathway_activity - pathway_activity.mean(axis=1).values.reshape(-1, 1)) / \
                           pathway_activity.std(axis=1).values.reshape(-1, 1)
    
    # Create color bar labels for groups
    group_values = clinical_data[group_by].values
    unique_groups = sorted([g for g in set(group_values) if pd.notna(g)], key=str)
    group_colors = {}
    colors = px.colors.qualitative.Plotly
    for i, group in enumerate(unique_groups):
        group_colors[group] = colors[i % len(colors)]
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=pathway_activity_norm.values,
        x=pathway_activity_norm.columns,
        y=pathway_activity_norm.index,
        colorscale='RdBu_r',
        zmid=0,
        colorbar=dict(
            title="Activity<br>(z-score)",
            thickness=20,
            len=0.7
        ),
        hoverongaps=False,
        hovertemplate='<b>%{y}</b><br>Sample: %{x}<br>Activity: %{z:.3f}<extra></extra>'
    ))
    
    # Update layout
    fig.update_layout(
        title=f"Pathway Activity Heatmap - Grouped by {group_by}",
        xaxis=dict(
            title=f"Patient Samples (ordered by {group_by})",
            side="bottom",
            tickangle=90
        ),
        yaxis=dict(
            title="Pathways (347 KEGG)",
            autorange="reversed"
        ),
        width=figsize[0],
        height=figsize[1],
        font=dict(size=10),
        margin=dict(l=300, r=150, t=80, b=150)
    )
    
    logger.info(f"✓ Created heatmap: {pathway_activity.shape[0]} pathways × {pathway_activity.shape[1]} samples")
    return fig


def plot_pathway_violin(pathway_activity: pd.DataFrame,
                       clinical_data: pd.DataFrame,
                       pathway: str,
                       group_by: str) -> go.Figure:
    """
    Create violin plot: pathway activity distribution by clinical group.
    
    Shows distribution of pathway activity within each group (e.g., by gender or stage).
    Points colored by individual samples.
    
    Parameters:
    -----------
    pathway_activity : pd.DataFrame
        Pathways × samples matrix
    clinical_data : pd.DataFrame
        Sample × clinical variables
    pathway : str
        Pathway name to visualize
    group_by : str
        Clinical variable to group by
        
    Returns:
    --------
    go.Figure : Plotly violin plot
    """
    if pathway not in pathway_activity.index:
        raise ValueError(f"Pathway '{pathway}' not in pathway_activity")
    
    logger.info(f"Creating violin plot: {pathway} by {group_by}")
    
    # Extract pathway activity for this pathway
    pathway_vals = pathway_activity.loc[pathway]
    
    # Align with clinical data
    common_samples = list(set(pathway_vals.index) & set(clinical_data.index))
    pathway_vals = pathway_vals[common_samples]
    clinical_groups = clinical_data.loc[common_samples, group_by]
    
    # Create dataframe for plotting
    plot_data = pd.DataFrame({
        'pathway_activity': pathway_vals.values,
        'group': clinical_groups.values,
        'sample': pathway_vals.index
    })
    
    # Create violin plot
    fig = px.violin(plot_data, 
                   x='group', 
                   y='pathway_activity',
                   points='all',
                   box=True,
                   title=f"{pathway} Activity by {group_by}",
                   labels={
                       'pathway_activity': 'Pathway Activity (Mean Expression)',
                       'group': group_by
                   },
                   hover_data={'sample': True, 'group': True, 'pathway_activity': ':.4f'})
    
    fig.update_layout(
        width=800,
        height=500,
        font=dict(size=12),
        showlegend=False,
        xaxis=dict(categoryorder='total ascending')
    )
    
    logger.info(f"✓ Created violin plot: {len(plot_data)} samples, {len(set(plot_data['group']))} groups")
    return fig


def plot_hub_genes_bar(hub_genes_df: pd.DataFrame,
                      top_n: int = 15) -> go.Figure:
    """
    Create bar chart of hub gene scores.
    
    Parameters:
    -----------
    hub_genes_df : pd.DataFrame
        Hub genes dataframe with columns: gene, hub_score, degree, betweenness, expr_variance
    top_n : int
        Number of top genes to display
        
    Returns:
    --------
    go.Figure : Plotly bar chart
    """
    logger.info(f"Creating hub gene bar chart: top {top_n} genes")
    
    top_genes = hub_genes_df.head(top_n).sort_values('hub_score', ascending=True)
    
    # Create bar chart with color gradient
    fig = go.Figure(data=[
        go.Bar(
            y=top_genes['gene'],
            x=top_genes['hub_score'],
            orientation='h',
            marker=dict(
                color=top_genes['hub_score'],
                colorscale='Viridis',
                showscale=True,
                colorbar=dict(title='Hub Score')
            ),
            hovertemplate='<b>%{y}</b><br>Hub Score: %{x:.4f}<extra></extra>',
            text=top_genes['hub_score'].round(4),
            textposition='outside'
        )
    ])
    
    fig.update_layout(
        title=f"Top {top_n} Hub Genes by Hub Score",
        xaxis_title="Hub Score",
        yaxis_title="Gene",
        height=400 + top_n * 15,
        width=600,
        margin=dict(l=150, r=100, t=80, b=50)
    )
    
    logger.info(f"✓ Created bar chart: {len(top_genes)} genes")
    return fig


def plot_differential_pathways(diff_results: pd.DataFrame,
                              top_n: int = 20) -> go.Figure:
    """
    Create visualization of differential pathway analysis results.
    
    Shows -log10(p-value) for top differential pathways.
    Colored by FDR significance.
    
    Parameters:
    -----------
    diff_results : pd.DataFrame
        Differential pathway analysis results with columns:
        pathway, pvalue, padj, significant
    top_n : int
        Number of top pathways to display
        
    Returns:
    --------
    go.Figure : Plotly bar/scatter plot
    """
    logger.info(f"Creating differential pathway plot: top {top_n} pathways")
    
    # Sort by p-value and take top
    top_results = diff_results.nsmallest(top_n, 'pvalue').sort_values('pvalue', ascending=True)
    
    # Calculate -log10(p-value)
    top_results['neg_log_pval'] = -np.log10(top_results['pvalue'])
    
    # Color by FDR significance
    colors = ['red' if sig else 'blue' for sig in top_results['significant']]
    
    fig = go.Figure(data=[
        go.Bar(
            y=top_results['pathway'],
            x=top_results['neg_log_pval'],
            orientation='h',
            marker=dict(color=colors),
            hovertemplate='<b>%{y}</b><br>-log10(p): %{x:.2f}<br>p-value: %{customdata:.2e}<extra></extra>',
            customdata=top_results['pvalue']
        )
    ])
    
    # Add significance threshold line
    fig.add_vline(x=-np.log10(0.05), line_dash="dash", line_color="gray",
                 annotation_text="p=0.05", annotation_position="top right")
    
    fig.update_layout(
        title=f"Top {top_n} Differential Pathways (p-value ranked)",
        xaxis_title="-log10(p-value)",
        yaxis_title="Pathway",
        height=400 + top_n * 15,
        width=900,
        margin=dict(l=300, r=50, t=80, b=50),
        legend=dict(
            x=0.7, y=0.95,
            bgcolor='rgba(255,255,255,0.8)',
            bordercolor='black',
            borderwidth=1
        )
    )
    
    # Add legend
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                            marker=dict(size=10, color='red'),
                            name='FDR < 0.05'))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                            marker=dict(size=10, color='blue'),
                            name='p < 0.05 (not FDR sig)'))
    
    logger.info(f"✓ Created differential pathway plot: {len(top_results)} pathways")
    return fig


def plot_pathway_comparison_boxplot(pathway_activity: pd.DataFrame,
                                   clinical_data: pd.DataFrame,
                                   pathways: list,
                                   group_by: str) -> go.Figure:
    """
    Create comparison plot of multiple pathways across groups.
    
    Rows: Pathways
    Columns: Groups
    Color: Activity level
    
    Parameters:
    -----------
    pathway_activity : pd.DataFrame
        Pathways × samples matrix
    clinical_data : pd.DataFrame
        Sample × clinical variables
    pathways : list
        List of pathway names to compare
    group_by : str
        Clinical variable to group by
        
    Returns:
    --------
    go.Figure : Plotly subplots with box plots
    """
    logger.info(f"Creating pathway comparison for {len(pathways)} pathways")
    
    # Get unique groups
    unique_groups = sorted([g for g in set(clinical_data[group_by]) if pd.notna(g)], key=str)
    
    # Create subplots
    n_pathways = len(pathways)
    n_cols = 3
    n_rows = (n_pathways + n_cols - 1) // n_cols
    
    fig = make_subplots(
        rows=n_rows, cols=n_cols,
        subplot_titles=[p[:30] + "..." if len(p) > 30 else p for p in pathways],
        vertical_spacing=0.12,
        horizontal_spacing=0.12
    )
    
    for idx, pathway in enumerate(pathways):
        if pathway not in pathway_activity.index:
            continue
            
        row = idx // n_cols + 1
        col = idx % n_cols + 1
        
        # Get data for this pathway
        pathway_vals = pathway_activity.loc[pathway]
        common_samples = list(set(pathway_vals.index) & set(clinical_data.index))
        
        pathway_vals = pathway_vals[common_samples]
        groups = clinical_data.loc[common_samples, group_by]
        
        # Add box plot for each group
        for group in unique_groups:
            mask = groups == group
            y_data = pathway_vals[mask].values
            
            fig.add_trace(
                go.Box(y=y_data, name=group, boxmean='sd'),
                row=row, col=col
            )
    
    fig.update_layout(
        title_text=f"Pathway Activity Comparison Across {group_by}",
        height=300 * n_rows,
        width=400 * n_cols,
        showlegend=True
    )
    
    logger.info(f"✓ Created pathway comparison plot")
    return fig


# ============================================================================
# Utility functions for creating combined visualizations
# ============================================================================

def create_pathway_summary_dashboard(pathway_activity: pd.DataFrame,
                                    clinical_data: pd.DataFrame,
                                    diff_results: pd.DataFrame,
                                    hub_genes_dict: dict,
                                    pathway: str) -> go.Figure:
    """
    Create comprehensive dashboard for single pathway analysis.
    
    Combines:
    - Violin plot of activity by group
    - Hub gene rankings
    - Differential expression results
    
    Parameters:
    -----------
    pathway : str
        Pathway to analyze
        
    Returns:
    --------
    go.Figure : Multi-panel Plotly figure
    """
    logger.info(f"Creating summary dashboard for {pathway}")
    
    # Get unique groups for grouping (assume first variable in clinical_data)
    group_by = 'ajcc_pathologic_stage' if 'ajcc_pathologic_stage' in clinical_data.columns else 'gender'
    
    # Create subplots
    from plotly.subplots import make_subplots
    
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=('Pathway Activity', 'Hub Genes', 'Differential Expression'),
        specs=[[{'type': 'box'}, {'type': 'bar'}, {'type': 'scatter'}]]
    )
    
    # 1. Violin plot
    pathway_vals = pathway_activity.loc[pathway]
    common_samples = list(set(pathway_vals.index) & set(clinical_data.index))
    pathway_vals = pathway_vals[common_samples]
    groups = clinical_data.loc[common_samples, group_by]
    
    for group in sorted([g for g in set(groups) if pd.notna(g)], key=str):
        mask = groups == group
        fig.add_trace(
            go.Box(y=pathway_vals[mask].values, name=group),
            row=1, col=1
        )
    
    # 2. Hub genes
    if pathway in hub_genes_dict:
        hub_genes = hub_genes_dict[pathway].head(10)
        fig.add_trace(
            go.Bar(y=hub_genes['gene'], x=hub_genes['hub_score'], 
                  orientation='h', name='Hub Genes'),
            row=1, col=2
        )
    
    # 3. Differential results
    if not diff_results.empty:
        pathway_diff = diff_results[diff_results['pathway'] == pathway]
        if not pathway_diff.empty:
            fig.add_trace(
                go.Scatter(y=pathway_diff['pathway'],
                          x=-np.log10(pathway_diff['pvalue']),
                          mode='markers',
                          name='p-value'),
                row=1, col=3
            )
    
    fig.update_layout(height=400, width=1400, showlegend=True)
    
    logger.info(f"✓ Created dashboard for {pathway}")
    return fig


if __name__ == "__main__":
    print("Pathway Visualization Module")
    print("=" * 60)
    print("\nUsage:")
    print("  from pathway_visualizations import plot_pathway_activity_heatmap")
    print("  fig = plot_pathway_activity_heatmap(pathway_activity, clinical_df, 'gender')")
    print("  fig.show()")
