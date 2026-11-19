# 第 1 阶段合成验证

## 首次运行 · FAIL

[运行 44553697477](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44553697477) 在 MATLAB R2024b Update 10 上运行了固定版本的作者示例。种子为 2026。两次原始 AR(3) 计算和仅记录优化器输出的复跑得到同一 (d_2)。SQP 求解器返回正的退出标志。

额外的 AR(1) 解析对照验证了方向。生成系数从 0.50 增至 0.95 时，估计的 (d_2) 从 0.2106864104 降至 0.01892511488 bits/sample。远离临界的例子与拟合系数对应解析值的绝对差为 (2.2357 \times 10^{-8}) bits/sample，超过测试最初设置的 (10^{-8}) 阈值。MATLAB 返回非零退出码，真实数据分析等待复验。

复验采用作者函数积分使用的精度准则。绝对误差 (10^{-10}) nats 换算为 bits，相对误差使用 MATLAB `integral` 的默认值 (10^{-6})。两者按积分器的组合规则计算允许误差。验证对拟合 AR(1) 系数使用 (d_2=\frac12\log_2[2/(1+\hat\phi)])。另用独立的绝对阈值 0.015 bits/sample 检查有限采样估计与生成系数对应的解析值。[MathWorks 积分误差准则](https://www.mathworks.com/help/matlab/ref/integral.html)
