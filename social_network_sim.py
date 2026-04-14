"""
社交网络仿真模块
Social Network Simulation Module
整合自demo目录下的SIS和community代码
"""

import networkx as nx
import numpy as np
import plotly.graph_objects as go
from typing import Tuple, List
import random

class SocialNetworkSimulator:
    """社交网络仿真器（不依赖PyTorch）"""
    
    def __init__(self):
        random.seed(0)
        np.random.seed(0)
    
    def build_community_network(self, N=100, c=5, k1=10, Z_in=8) -> Tuple[nx.Graph, List[List[int]]]:
        """
        构建社区网络
        
        参数：
        - N: 节点总数
        - c: 社区数量
        - k1: 平均度数
        - Z_in: 社区内平均连接数
        
        返回：
        - G: NetworkX图对象
        - communities: 社区划分列表
        """
        nodes_per_comm = N // c
        communities = [list(range(i * nodes_per_comm, (i + 1) * nodes_per_comm)) for i in range(c)]
        
        Z_out = k1 - Z_in
        pin = Z_in / (nodes_per_comm - 1)
        pout = Z_out / (N - nodes_per_comm)
        
        G = nx.Graph()
        G.add_nodes_from(range(N))
        
        # 添加边
        for i in range(N):
            for j in range(i + 1, N):
                same_comm = any(i in comm and j in comm for comm in communities)
                p = pin if same_comm else pout
                if random.random() < p:
                    G.add_edge(i, j)
        
        return G, communities
    
    def SIS_simulation(self, G: nx.Graph, beta=0.05, gamma=0.2, ini=0.01, max_step=100) -> Tuple[np.ndarray, np.ndarray]:
        """
        SIS传播仿真（简化版，不使用超边）
        
        参数：
        - G: NetworkX图对象
        - beta: 感染率
        - gamma: 恢复率
        - ini: 初始感染比例
        - max_step: 仿真步数
        
        返回：
        - infected_density: 每步的感染密度
        - node_state: 每步每个节点的状态矩阵
        """
        N = G.number_of_nodes()
        node_state = np.zeros((max_step, N))
        
        # 初始感染
        aa = int(N * ini)
        initial_infect = np.random.choice(N, aa, replace=False)
        node_state[0, initial_infect] = 1
        
        # 构建邻接矩阵
        A = nx.adjacency_matrix(G).todense()
        A = np.array(A)
        
        # 仿真传播过程
        for t in range(max_step - 1):
            # 恢复项
            recovery = (1 - gamma) * node_state[t]
            
            # 感染项（来自邻居）
            infection_prob = np.sum(beta * A * node_state[t], axis=1)
            
            # 总感染概率
            prob = recovery + (1 - node_state[t]) * infection_prob
            prob = np.clip(prob, 0, 1)
            
            # 随机决定下一状态
            rand_vals = np.random.rand(N)
            node_state[t+1] = (rand_vals <= prob).astype(float)
        
        # 计算感染密度
        infected_density = node_state.sum(axis=1) / N
        
        return infected_density, node_state
    
    def get_network_metrics(self, G: nx.Graph, communities: List[List[int]]) -> dict:
        """
        计算网络拓扑指标

        参数:
            G: NetworkX图对象
            communities: 社区划分列表

        返回:
            dict: 包含 density, avg_degree, clustering, modularity 等指标
        """
        n = G.number_of_nodes()
        if n == 0:
            return {"density": 0, "avg_degree": 0, "clustering": 0, "modularity": 0}

        density = nx.density(G)
        avg_degree = sum(dict(G.degree()).values()) / n
        clustering = nx.average_clustering(G) if n > 2 else 0

        # 模块度
        try:
            partition = [{node for node in comm if node in G} for comm in communities]
            partition = [s for s in partition if len(s) > 0]
            modularity = nx.community.modularity(G, partition) if partition else 0
        except Exception:
            modularity = 0

        return {
            "density": round(density, 4),
            "avg_degree": round(avg_degree, 2),
            "clustering": round(clustering, 4),
            "modularity": round(modularity, 4),
        }

    def visualize_community_network(self, G: nx.Graph, communities: List[List[int]],
                                    title="社区网络结构") -> go.Figure:
        """
        可视化社区网络
        
        参数：
        - G: NetworkX图对象
        - communities: 社区划分
        - title: 图表标题
        
        返回：
        - Plotly图表对象
        """
        # 使用spring布局
        pos = nx.spring_layout(G, seed=42, k=0.5, iterations=50)
        
        # 为每个社区分配颜色
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', 
                 '#F7DC6F', '#BB8FCE', '#85C1E2', '#F8B739', '#52B788']
        
        node_colors = []
        node_community = {}
        for comm_idx, comm in enumerate(communities):
            for node in comm:
                node_community[node] = comm_idx
                node_colors.append(colors[comm_idx % len(colors)])
        
        # 绘制边
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
        
        # 绘制节点
        node_x = [pos[node][0] for node in G.nodes()]
        node_y = [pos[node][1] for node in G.nodes()]
        
        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            hoverinfo='text',
            marker=dict(
                color=node_colors,
                size=10,
                line=dict(width=1, color='white')
            )
        )
        
        # 添加悬停信息
        node_text = []
        for node in G.nodes():
            degree = G.degree(node)
            comm = node_community[node]
            node_text.append(f'节点: {node}<br>度数: {degree}<br>社区: {comm}')
        
        node_trace.text = node_text
        
        # 创建图表
        fig = go.Figure(data=[edge_trace, node_trace])
        
        fig.update_layout(
            title=dict(
                text=title + f"<br><sub>节点数: {G.number_of_nodes()}, 边数: {G.number_of_edges()}, 社区数: {len(communities)}</sub>",
                x=0.5,
                xanchor='center'
            ),
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=60),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white',
            height=600
        )
        
        return fig
    
    def visualize_infection_spread(self, infected_density: np.ndarray, 
                                   max_step: int, title="SIS传播动态") -> go.Figure:
        """
        可视化感染传播曲线
        
        参数：
        - infected_density: 感染密度数组
        - max_step: 仿真步数
        - title: 图表标题
        
        返回：
        - Plotly图表对象
        """
        steps = np.arange(max_step)
        
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=steps,
            y=infected_density,
            mode='lines+markers',
            name='感染密度',
            line=dict(color='#FF6B6B', width=2),
            marker=dict(size=4)
        ))
        
        fig.update_layout(
            title=dict(
                text=title,
                x=0.5,
                xanchor='center'
            ),
            xaxis=dict(
                title='时间步',
                showgrid=True,
                gridwidth=1,
                gridcolor='#E0E0E0'
            ),
            yaxis=dict(
                title='感染密度',
                showgrid=True,
                gridwidth=1,
                gridcolor='#E0E0E0',
                range=[0, 1]
            ),
            hovermode='x unified',
            plot_bgcolor='white',
            height=400,
            margin=dict(l=60, r=20, t=60, b=60)
        )
        
        return fig
    
    def visualize_infection_snapshot(self, G: nx.Graph, communities: List[List[int]], 
                                     node_state: np.ndarray, step: int) -> go.Figure:
        """
        可视化某一时刻的感染状态
        
        参数：
        - G: NetworkX图对象
        - communities: 社区划分
        - node_state: 节点状态矩阵
        - step: 时间步
        
        返回：
        - Plotly图表对象
        """
        # 使用spring布局
        pos = nx.spring_layout(G, seed=42, k=0.5, iterations=50)
        
        # 绘制边
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
        
        # 绘制节点（根据感染状态着色）
        node_x = [pos[node][0] for node in G.nodes()]
        node_y = [pos[node][1] for node in G.nodes()]
        
        # 感染状态：红色=感染，蓝色=健康
        node_colors = ['#FF6B6B' if node_state[step, node] == 1 else '#4ECDC4' 
                      for node in G.nodes()]
        
        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            hoverinfo='text',
            marker=dict(
                color=node_colors,
                size=10,
                line=dict(width=1, color='white')
            )
        )
        
        # 添加悬停信息
        node_text = []
        node_community = {}
        for comm_idx, comm in enumerate(communities):
            for node in comm:
                node_community[node] = comm_idx
        
        for node in G.nodes():
            status = "感染" if node_state[step, node] == 1 else "健康"
            comm = node_community.get(node, -1)
            node_text.append(f'节点: {node}<br>状态: {status}<br>社区: {comm}')
        
        node_trace.text = node_text
        
        # 创建图表
        fig = go.Figure(data=[edge_trace, node_trace])
        
        infected_count = int(node_state[step].sum())
        total_nodes = len(node_state[step])
        
        fig.update_layout(
            title=dict(
                text=f"感染状态快照 (步数: {step})<br><sub>红色=感染 ({infected_count}/{total_nodes}), 蓝色=健康</sub>",
                x=0.5,
                xanchor='center'
            ),
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=60),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white',
            height=600
        )
        
        return fig

