# 项6：公开仓库 v0.4.1 + Zenodo 发布清单

更新时间：2026-11-05（公开导入）  
状态：**公开仓库为 https://github.com/YancyYang02/hapGWAS ，版本 0.4.1。** 版本 DOI 为 https://doi.org/10.5281/zenodo.23984394 （记录 https://zenodo.org/records/23984394 ）；概念 DOI 为 https://doi.org/10.5281/zenodo.23984393 。下方清单是导入前的本地核对记录。

## 0. 当前本地版本状态（已核对）

| 文件 | 状态 |
|---|---|
| `VERSION` | **0.4.1** |
| `pyproject.toml` → `[project].version` | **0.4.1** |
| `hapgwas/__init__.py` → `__version__` | **0.4.1** |
| `CITATION.cff` → `version` | **0.4.1**（已由 0.3.1 对齐） |
| `README.md` 标题 | **hapGWAS (v0.4.1)** |
| `CHANGELOG.md` | 已写 **[0.4.1]** 条目（`--grm` / LMM；SACJ 文档降级；hardening 指针） |
| `.gitignore` | 已补（venv / `__pycache__` / `*.egg-info` / `.DS_Store` / 大 VCF·GRM 等） |
| 功能 | CLI `--grm` / `--kinship` → spectral LMM（GCTA `.grm.bin+.id` / npy·npz / 文本方阵） |

工作树：`/workspace/hapgwas`（并同步 Desktop `hapgwas/` 对应文件）。  
**尚无**本地 `.git`（工作区未初始化远程）；**禁止**在未拿到公开 URL 前 push。

## 1. GitHub 账号现状（只读确认）

- 已连接 Cursor 的 GitHub 账号：**YancyYang02**（`https://github.com/YancyYang02`）
- 当前 **public_repos = 0**（无公开仓库）
- **Cursor / Grok Bot 不能替你创建 GitHub 仓库**（无 create-repo 授权路径；也不应绕过）
- **需要你本人**在 GitHub 网页新建一个 **空的 Public** 仓库，建议名：`hapGWAS`  
  - 勾选 Public  
  - **不要**勾选 “Add a README / .gitignore / license”（保持空仓，便于首次 push）  
  - 建好后把仓库 URL 发给助手，例如：`https://github.com/YancyYang02/hapGWAS.git`  
  - **在收到该 URL 之前，本文档与 Availability 一律不写死仓库地址**

## 2. 拿到公开仓库 URL 之后的操作清单（待你授权执行）

> 下列步骤**尚未执行**。助手在收到明确 URL 与授权前，不做 `git push`、不创建 Release、不联系 Zenodo、不编造 DOI。

### A. 本地 git 初始化并推送（示例；URL 以你提供的为准）

- [ ] `cd` 到发布用源码树（建议干净拷贝，排除 `.venv*`、大 VCF、HPC 结果）
- [ ] `git init`（若尚无 `.git`）
- [ ] `git add` 源码与文档（遵循 `.gitignore`）
- [ ] `git commit -m "Release v0.4.1: GRM/LMM --grm and App Note hardening docs"`
- [ ] `git branch -M main`
- [ ] `git remote add origin <YOUR_PUBLIC_REPO_URL>`
- [ ] `git push -u origin main`
- [ ] `git tag -a v0.4.1 -m "hapGWAS 0.4.1"`
- [ ] `git push origin v0.4.1`

### B. GitHub Release

- [ ] 在 GitHub → Releases → “Draft a new release”
- [ ] Choose tag：**v0.4.1**
- [ ] 标题建议：`hapGWAS v0.4.1`
- [ ] 正文可粘贴 `CHANGELOG.md` 中 `[0.4.1]` 段
- [ ] Publish release（Public）

### C. Zenodo ↔ GitHub

- [ ] 用同一 GitHub 账号登录 Zenodo（或 CERN Zenodo）
- [ ] 启用 GitHub 授权 / “Enable GitHub” 并选中公开仓库 `hapGWAS`（名称以实际为准）
- [ ] 对 **v0.4.1** Release 触发存档，由 Zenodo **铸造 DOI**
- [ ] **禁止编造 DOI**；仅使用 Zenodo 页面显示的正式 DOI

### D. 回填 DOI 与 Availability（DOI 到手后）

- [x] 更新 `CITATION.cff`：`doi` 为版本 DOI 10.5281/zenodo.23984394；`identifiers` 另记概念 DOI 10.5281/zenodo.23984393
- [x] 更新 `APPLICATION_NOTE_en.md` Availability：仓库 URL 与版本 DOI 已写入
- [x] 更新 README “git clone” 为 https://github.com/YancyYang02/hapGWAS
- [ ] 可选：再打一个文档回填小版本或在同一 tag 说明中注明 DOI（按期刊要求）

## 3. 明确不做 / 已停止的事项

- 不发明 / 不占位假 DOI（例如 `10.5281/zenodo.xxxxxx`）
- 不伪造仓库 URL
- 不执行 `git push`、不调用 GitHub “create repository”
- 不登录 Zenodo、不代填密码或 token
- SACJ 保持文档降级（主表 joint / ACJ / 正式 GHap allele / SNP）；代码可保留

## 4. 你现在需要回复的一句话

请创建空的公开仓库后回复，例如：

> 公开仓库 URL：`https://github.com/YancyYang02/hapGWAS.git`，授权推送 `main` + tag `v0.4.1` 并起草 GitHub Release。

收到后再继续项6后半段。
