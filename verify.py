"""
快速验证脚本 - 检查应用是否能正常初始化
"""

import sys

print("🔍 验证环境...")

# 检查Python版本
print(f"Python版本: {sys.version}")

# 检查依赖
try:
    import gradio as gr
    print(f"✅ Gradio {gr.__version__}")
except ImportError as e:
    print(f"❌ Gradio未安装: {e}")
    sys.exit(1)

try:
    import plotly
    print(f"✅ Plotly {plotly.__version__}")
except ImportError as e:
    print(f"❌ Plotly未安装: {e}")
    sys.exit(1)

try:
    import networkx as nx
    print(f"✅ NetworkX {nx.__version__}")
except ImportError as e:
    print(f"❌ NetworkX未安装: {e}")
    sys.exit(1)

try:
    import pandas as pd
    print(f"✅ Pandas {pd.__version__}")
except ImportError as e:
    print(f"❌ Pandas未安装: {e}")
    sys.exit(1)

try:
    import numpy as np
    print(f"✅ NumPy {np.__version__}")
except ImportError as e:
    print(f"❌ NumPy未安装: {e}")
    sys.exit(1)

print("\n🔍 验证应用代码...")

try:
    # 导入应用模块（不运行）
    import importlib.util
    spec = importlib.util.spec_from_file_location("app_full", "app_full.py")
    module = importlib.util.module_from_spec(spec)
    
    print("✅ 应用代码可以导入")
    
    # 测试数据库初始化
    print("\n🔍 测试数据库...")
    spec.loader.exec_module(module)
    
    db = module.GeneNetworkDatabase()
    print(f"✅ 疾病数: {len(db.get_all_diseases())}")
    print(f"✅ 通路数: {len(db.get_all_pathways())}")
    print(f"✅ 基因数: {len(db.genes)}")
    
    print("\n✅ 所有验证通过！")
    print("\n📝 启动命令:")
    print("   python3 app_full.py")
    print("\n🌐 访问地址:")
    print("   http://127.0.0.1:7860")
    
except Exception as e:
    print(f"❌ 验证失败: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

