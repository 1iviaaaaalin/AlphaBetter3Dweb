# AlphaBetter

蛋白质结构稳定性分析与 3D 可视化平台。

## 使用方式

部署完成后，用户只需要打开网页：

1. 选择 `.pdb` 文件
2. 点击“开始分析”
3. 等待云端计算
4. 直接查看 3D 蛋白质结构
5. 红色 = 不稳定，黄色 = 接近中性，绿色 = 稳定
6. 可下载带 B-factor 稳定性分数的彩色 PDB

本项目已经把 React 前端和 FastAPI 后端合并到同一个 Web Service 中，因此最终只需要一个网址，不需要用户安装 VS Code、Python、Node.js 或启动本地服务。

## 部署

项目使用 Docker，可直接部署到 Render。Render 支持从 Git 仓库中的 Dockerfile 构建 Web Service；服务需要监听 `0.0.0.0`，本项目已经在 Dockerfile 中处理。部署后会获得一个 `onrender.com` 公网网址。

### GitHub

将整个项目文件夹上传到一个 GitHub 仓库，根目录必须能看到：

- `Dockerfile`
- `app.py`
- `utils.py`
- `calc_stability.py`
- `data/`
- `frontend/`
- `render.yaml`

### Render

在 Render 中创建 **Web Service**，连接这个 GitHub 仓库，Runtime 选择 **Docker**。Dockerfile 在项目根目录，所以无需修改路径。

健康检查地址：`/api/health`

## 算法说明

稳定性计算仍使用项目现有的二肽/五肽参考数据库和打分逻辑。

二级结构识别优先读取 PDB 自带的 `HELIX` / `SHEET` 注释；如果 PDB 没有这些注释，则使用项目内置的 USSA 风格几何/氢键规则作为后备方案。这样可以避免当前测试结构被误判成全部 Coil，从而导致所有稳定性分数都是 0。
