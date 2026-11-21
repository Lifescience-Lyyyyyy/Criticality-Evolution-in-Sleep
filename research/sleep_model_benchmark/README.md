# 已有睡眠计算模型的临界距离基准测试

本目录冻结已发表模型参数，先验证原论文生理指标，再用固定 Sooter MATLAB 估计器比较计数观测的 d₂。正式规格为 40 ms、AR(5)、Yule–Walker、β=2，AR(10) 为敏感性规格。连续脑电代理单列。三个模拟种子是三次技术重复，大鼠记录按动物聚合。

- [中文科学报告](reports/benchmark.md)
- [模型审计与可行性排名](reports/model_audit.md) · [机器可读登记表](registry.json)
- [原始生理复现](results/original_validation.json) · [模拟 d₂ 和 QC](results/matlab_benchmark.json) · [实测等时 mPFC](results/empirical_240s.json)
- [连续信号验证](results/costa_continuous.json) · [分支模型 RQA 失败记录](results/branching_validation.json) · [Krishnan 原始代码冒烟](results/krishnan_smoke.json)
- [生理验证图](figures/physiological_validation.pdf) · [观测敏感性图](figures/observation_sensitivity.pdf) · [多尺度图](figures/multiscale_sensitivity.pdf)

## 已执行的复现命令

在 WSL 的 Git 检出根目录执行。下列 `$SLEEP_BENCH` 指向仓库外的缓存或已忽略的 `logs/` 子目录。需要 Python 3.12、GNU g++，以及 NumPy、SciPy、h5py 和 remfile。Actions 固定 `numpy==2.2.6 scipy==1.15.3 h5py==3.13.0 remfile==0.1.16`。生成图另需 Matplotlib。本地验证用 Python 3.12.3、GNU g++ 13.3.0、NumPy 2.5.5、SciPy 1.18.1、Matplotlib 3.11.2。

~~~bash
export PYTHONPATH=research/sleep_model_benchmark
mkdir -p "$SLEEP_BENCH"
python -m unittest discover -s research/sleep_model_benchmark/tests
python "$PYTHONPATH"/validate_original.py --out "$SLEEP_BENCH/levenstein"
python "$PYTHONPATH"/prepare_fits.py \
  --simulation "$SLEEP_BENCH/levenstein" --out "$SLEEP_BENCH/counts"
python "$PYTHONPATH"/validate_branching.py --out "$SLEEP_BENCH/branching"

git clone https://github.com/ModelDBRepository/226474.git "$SLEEP_BENCH/Costa"
git -C "$SLEEP_BENCH/Costa" checkout eb9fad3bf9e19e3e7f3018ec471c93cfb0f31ebe
python "$PYTHONPATH"/validate_costa.py --source "$SLEEP_BENCH/Costa" \
  --out "$SLEEP_BENCH/continuous" --seconds 20
python "$PYTHONPATH"/validate_costa.py --source "$SLEEP_BENCH/Costa" \
  --out "$SLEEP_BENCH/continuous-180s" --seconds 180

git show 5fe4ab5efbb4888dd63f7b2972d4195a628c83b0:data/fixtures/sleep_dandi_input.mat \
  > "$SLEEP_BENCH/archived-counts.mat"
python "$PYTHONPATH"/prepare_empirical.py \
  --archive "$SLEEP_BENCH/archived-counts.mat" --out "$SLEEP_BENCH/empirical"
python "$PYTHONPATH"/plot_results.py --preview "$SLEEP_BENCH/figure-previews"
~~~

Krishnan 冒烟在源码提交 `f672f122a07a57a8c8d54fca81d8309431e7b848` 执行以下命令。预置 `one_state.h=4` 保持原样。提交中的 C++ 改动仅为本仓库的独立模拟器及 Costa 输出接口，第三方源码保存在外部缓存。

~~~bash
make -j2 mindcurrent
make network network_config=network.cfg
mkdir -p "$SLEEP_BENCH/krishnan/recording"
sed -E 's/^tmax .*/tmax 1000/; s/^num_mp_threads .*/num_mp_threads 3/' \
  params.txt > "$SLEEP_BENCH/krishnan/params.txt"
/usr/bin/time -v timeout 120 ./mindcurrent "$SLEEP_BENCH/krishnan/params.txt" \
  "$SLEEP_BENCH/krishnan/recording" connection_info2
~~~

正式拟合由 GitHub Actions 触发，分别选择 `model` 和 `empirical`。工作流现场生成计数输入，验证摘要后用 MATLAB R2024b 运行 `fit_matlab`。原始轨迹、计数数组、第三方构建产物留在缓存。工作流保留输出诊断和校验，Git 仅保存审核后的摘要、报告、图及代码。校验结果时将下载的原始 artifact ZIP 保留为 `actions.zip` 并解包至同一目录，再使用 `summarize.py` 或 `summarize_empirical.py`。具体 ID、哈希、模型代码提交和拟合提交保存于正式结果 JSON。[模拟运行](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44895142233) · [实测等时运行](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44896254896)
