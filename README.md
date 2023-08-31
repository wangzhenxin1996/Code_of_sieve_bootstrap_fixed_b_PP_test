# Sieve Bootstrap for Fixed-$b$ Phillips-Perron Unit Root Test

This repository contains code and scripts related to the manuscript "the Sieve Bootstrap for the fixed-$b$ Phillips-Perron unit root test". The primary focus is on the `ztest.m` function for conducting the PP^b(fb) test, along with a series of simulation scripts and empirical test scripts.

## Table of Contents

- [PP^b(fb)test](#PP^b(fb)test)
- [Simulation](#simulation)
- [Empirical Test](#empirical-test)


## PP^b(fb)test


The `ztest.m` function performs the Phillips-Perron (PP) test, PP(fb) test proposed by Vogelsang and Wagner (2013) and PP^b(fb) test on a given time series data.

### Function Signature:

```matlab
[zalpha, p_zalpha, zt, p_zt] = ztest(y, dt, q, m, ic, B, b, type);
```
### Input Parameters:

- `y`: Original series for testing.
- `dt`: Detrending method:
  - `1` for one-step detrending,
  - `2` for two-step detrending,
  - `3` for GLS detrending.
- `q`: Highest power of deterministic linear time trend to be removed:
  - `-1` for no-detrending,
  - `0` for demeaning,
  - `1` for linear trend detrending.
- `m`: Method for PP test:
  - `1` for traditional PP test using Newey-West for long-term variance estimation,
  - `2` for fixed-b PP unit root test (PP(fb)) proposed by Vogelsang and Wagner (2013), using fixed-b HAR long-run variance estimation,
  - `4` for obtaining critical values of PP(fb) test using the Sieve Bootstrap method (PP^b(fb)).
- `ic`: Information criterion method for selecting the number of autocorrelation lags used in the Sieve Bootstrap method:
  - `"aic"` for AIC criterion,
  - `"bic"` for BIC criterion,
  - `"maic"` for Modified AIC criterion,
  - `"mbic"` for Modified BIC criterion.
- `B`: Number of times bootstrap samples are constructed in the Bootstrap method.
- `b`: Truncation parameter for fixed smoothing parameter HAR long-term variance estimation:
  - When `type` is `"bt"`, `"pr"`, or `"qs"`, `b` represents the truncation parameter.
  - When `type` is `"os"`, it represents the length of the orthogonal sequence used.
- `type`: Type of kernel or non-parametric method used when estimating long-term variance-covariance matrix:
  - `"bt"` for Bartlett kernel,
  - `"pr"` for Parzen kernel,
  - `"qs"` for Quadratic Spectral kernel,
  - `"os"` for Orthonormal Series method.

### Output Parameters:

- `zalpha`: Statistic obtained from the PP test.
- `zt`: Statistic obtained from the PP test.
- `p_zalpha`: Corresponding p-value for the `zalpha` statistic.
- `p_zt`: Corresponding p-value for the `zt` statistic.

### Reference:

Vogelsang, T. J. and Wagner, M. (2013). A Fixed-b Perspective on the Phillips-Perron Unit Root Tests. Econometric Theory, 29(3):609-628.

### Example Usage of `ztest.m`

To illustrate how to use the `ztest.m` function for conducting a Phillips-Perron unit root test, consider the following example:

```matlab
% Set a fixed random seed for consistent results
rng(0);

% Create a synthetic unit root dataset
T = 200;
unit_root_data = cumsum(randn(T, 1));

% Parameters for the test
detrending_method = 3;     % GLS detrending
trend_order = 1;           % linear detrending
test_method = 3;           % PP^b(fb) unit root test
information_criterion = 'aic'; % AIC for selecting autocorrelation lags
bootstrap_samples = 499;
truncation_param = 0.02;   % Truncation parameter for HAR
kernel_type = 'bt';        % Bartlett kernel

% Perform the PP unit root test
[zalpha, p_zalpha, zt, p_zt] = ztest(unit_root_data, detrending_method, trend_order, test_method, information_criterion, bootstrap_samples, truncation_param, kernel_type);

% Display the results
disp(['zalpha Statistic: ', num2str(zalpha)]);
disp(['p-value for zalpha Statistic: ', num2str(p_zalpha)]);
disp(['zt Statistic: ', num2str(zt)]);
disp(['p-value for zt Statistic: ', num2str(p_zt)]);
```

## Simulation

This section focuses on simulating the performance of three unit root tests using the `ztest.m` function in different scenarios. The simulations aim to analyze the size and power of three tests under various conditions.


### Size Table for Different Bootstrap Replications

The script `SizeTableforT_100withDifferentBootsam` simulates the PP^b(fb) test's size under different number of bootstrap replications. The results are stored in the file `SizeTableforT_100withDifferentBootsam`. This corresponds to Table 1 in the supplementary materials.

### Power Table for Different Bootstrap Replications

The script `PowerTableforT_100withDifferentBootsam` simulates the PP^b(fb) test's power under different number of bootstrap replications. The results are stored in the file `PowerTableforT_100withDifferentBootsam`. This corresponds to Table 2 in the supplementary materials.

### Size Curves for Different Sample Sizes

The scripts `sizecurveforT_50withdifferentb`, `sizecurveforT_100withdifferentb`, `sizecurveforT_200withdifferentb`, and `sizecurveforT_400withdifferentb` simulate size curves for various unit root tests under different sample sizes (T=50, 100, 200, 400) and varying truncation parameter `b` values. The simulation results are saved in the following files:
- `SizeCurveforT_50withDifferentB`
- `SizeCurveforT_100withDifferentB`
- `SizeCurveforT_200withDifferentB`
- `SizeCurveforT_400withDifferentB`

These size curves demonstrate how the size of the tests changes with different truncation parameter `b`. To visualize the results, the following scripts create plots for the size curves:
- `plotsizecurveforT_50withdifferentb`
- `plotsizecurveforT_100withdifferentb`
- `plotsizecurveforT_200withdifferentb`
- `plotsizecurveforT_400withdifferentb`

These plots correspond to Figures 1 to 8 in the supplementary materials. Additionally, Figures 2 to 7 in the main text are excerpts from Figures 1 to 8 in the supplementary materials.

Table 1 in the main text and Table 3 in the supplementary materials present the selected results for `b = 0.02` from the aforementioned size curves.

### Power Curve for Sample Sizes

The script `powercurveforsamplesize` simulates the power of the PP^b(fb) test under various sample sizes. The results are stored in the file `PowerCurveforSampleSizewithSizeAdjusted`. The script `plotpowercurveforsamplesize` visualizes the results as graphs.

These simulation scripts aim to provide insights into the behavior of the PP^b(fb) test's size and power under different conditions, sample sizes, and bootstrap iterations.

## Empirical Test

This section focuses on empirical tests conducted using the `ztest.m` function on real-world datasets. Two specific cases are considered: one related to U.S. inflation (`empcpi`) and the other related to Chinese stock indices (`empstock`).
