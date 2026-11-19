# 第 1 阶段合成验证

## 首次运行 · FAIL

[运行 44553697477](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44553697477) 在 MATLAB R2024b Update 10 上运行了固定版本的作者示例。种子为 2026。两次原始 AR(3) 计算和仅记录优化器输出的复跑得到同一 d₂。SQP 求解器返回正的退出标志。

额外的 AR(1) 解析对照验证了方向。生成系数从 0.50 增至 0.95 时，估计的 d₂ 从 0.2106864104 降至 0.01892511488 bits/sample。远离临界的例子与拟合系数对应解析值的绝对差为 2.2357 × 10⁻⁸ bits/sample，超过测试最初设置的 10⁻⁸ 阈值。MATLAB 返回非零退出码，真实数据分析等待复验。

复验采用作者函数积分使用的精度准则。绝对误差 10⁻¹⁰ nats 换算为 bits，相对误差使用 MATLAB `integral` 的默认值 10⁻⁶。两者按积分器的组合规则计算允许误差。验证对拟合 AR(1) 系数使用 $d_2=\frac12\log_2[2/(1+\hat\phi)]$。另用独立的绝对阈值 0.015 bits/sample 检查有限采样估计与生成系数对应的解析值。[MathWorks 积分误差准则](https://www.mathworks.com/help/matlab/ref/integral.html)

## 复验 · PASS

[运行 44554797385](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44554797385) 在项目提交 `eccaec9d3183e1c21148e6f53e3b4a332486ba18` 执行。未修改的作者函数两次均得到 `0.03465842415362974` bits/s。原始示例的 AR(3) 拟合系数为 `[0.7952113, 0.1972218, -0.0964218]`，估计残差标准差为 0.4465526。只记录诊断的临时副本复算一致，SQP 退出标志为 1，迭代次数为 3，约束违反量为 0。整项验证耗时 4.60 秒。

| 生成 AR(1) 系数 | 拟合系数 | MATLAB d₂（bits/sample） | 生成过程解析 d₂（bits/sample） |
| --- | --- | --- | --- |
| 0.50 | 0.49342743 | 0.21068641 | 0.20751875 |
| 0.95 | 0.94821077 | 0.01892511 | 0.01826294 |

两组有限样本各包含 10,000 个点，共用种子 2026 生成的一条高斯输入噪声，并丢弃前 2,000 点。每组为一条合成时间序列，图中显示单次估计和解析值。两组都通过拟合系数解析值的积分精度准则。可执行入口为 `scripts/ci/run_synthetic.m`。

[合成验证图 PDF](../figures/pnas/synthetic_d2.pdf) · [数值与诊断](../results/pnas/synthetic.json)
