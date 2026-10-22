# Archived qPCR plotting scripts

These three historical scripts were moved here from `qpcr/` and preserved unchanged. They are not a raw Ct analysis workflow.

`deepseek_r_20260827_605096.r` contains a demonstration ΔCt plot that adds random noise with `rnorm()`. The spread in that plot is synthetic; do not use it as experimental replicate data or as an uncertainty estimate.

The current `final_qpcr.r` also plots fixed relative-quantity values. Recalculate ΔCt and ΔΔCt from the original Ct measurements before treating qPCR summaries as a reproducible statistical result.
