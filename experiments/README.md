# Reproducing the EAlign experiments

Run every command below from the EAlign repository root. The commands use the
parameters reported in the paper.

## 0. Setup

Extract each released archive into the output directory shown below:

| Archive | Contents | Output directory |
| --- | --- | --- |
| `Manual_curated_benchmarks.zip` | MALISAM, MALIDUP, HOMSTRAD, and the two pair-list CSVs | `dataset/` |
| `SABmark.zip` | SABmark structures and SupFP/TwiFP JSONL records | `dataset/` |
| `Pairwise_Quality_Runtime_Benchmark.zip` | Saved per-method pairwise runs and `run_summary.csv` | `test_results/pairwise_runs/` |
| `Detailed_experiments_results.zip` | Saved outputs from Experiments 01, 02, and 03 | `test_results/` |

Download: [Zenodo](https://zenodo.org/uploads/23476316)

```bash
mkdir -p dataset test_results/pairwise_runs

unzip Manual_curated_benchmarks.zip -d dataset
unzip SABmark.zip -d dataset
unzip Pairwise_Quality_Runtime_Benchmark.zip \
  -d test_results/pairwise_runs
unzip Detailed_experiments_results.zip -d test_results

make -j
python3 -m pip install -r experiments/requirements.txt
```

Expected executables:

```text
binary/ealign
binary/baselines/DeepAlign/DeepAlign
binary/baselines/TMalign/TMalign
binary/baselines/DaliLite.v5/bin/dali.pl
binary/baselines/gtalign_alpha/Linux_installer_GPU/bin/gtalign
binary/baselines/gtalign_alpha/patched_cpu/bin/gtalign
binary/baselines/foldseek/bin/foldseek
```

The paper used DeepAlign 1.4, TM-align 20220412, DaliLite v5, GTalign
1.0.1, and Foldseek revision `463739e0014a`. Apply the supplied GTalign CPU
finalization patch before running the CPU experiments:

## 1. Pairwise alignment

Datasets: 130 MALISAM, 241 MALIDUP, and 1,032 HOMSTRAD pairs. Runtime
experiments use one outer worker.

Prepare HOMSTRAD:

```bash
python3 experiments/01_pairwise_alignment/split_HOMSTRAD.py \
  --homstrad-root dataset/HOMSTRAD \
  --output-root test_results/HOMSTRAD_split_structures
```

Run EAlign, DeepAlign, TM-align, Foldseek, and GTalign-CUDA:

```bash
for method in ealign deepalign tmalign foldseek gtalign_cuda; do
  python3 experiments/01_pairwise_alignment/run_pairwise_alignment.py \
    --method "$method" --workers 1 \
    --output-root test_results/pairwise_runs
done
```

Run GTalign-CPU with one and 20 internal threads:

```bash
python3 experiments/01_pairwise_alignment/run_pairwise_alignment.py \
  --method gtalign_patched_cpu --workers 1 --cpu-threads 1 \
  --output-root test_results/pairwise_runs

python3 experiments/01_pairwise_alignment/run_pairwise_alignment.py \
  --method gtalign_patched_cpu --workers 1 --cpu-threads 20 \
  --output-root test_results/pairwise_runs
```

Run DALI with one job:

```bash
for benchmark in malisam malidup homstrad; do
  python3 experiments/01_pairwise_alignment/run_homstrad_dali_direct.py \
    --benchmark "$benchmark" --jobs 1 \
    --output-root "test_results/pairwise_runs/dali/$benchmark"
done
```
We also provided runtime results obtained from a CentOS Linux 8 server with two
Intel Xeon Gold 5218 CPUs (32 physical cores, 64 logical CPUs):
these raw runs and their combined `run_summary.csv` are provided in
`Pairwise_Quality_Runtime_Benchmark.zip`.

Generate all paper statistics and figures:

```bash
python3 experiments/01_pairwise_alignment/evaluate_pairwise_alignment.py
python3 experiments/01_pairwise_alignment/analyze_pairwise_evaluation.py

python3 experiments/01_pairwise_alignment/summarize_pairwise_runtime.py
python3 experiments/01_pairwise_alignment/summarize_pairwise_runtime.py \
  --methods gtalign_cuda gtalign_patched_cpu \
            gtalign_patched_cpu_t20 foldseek \
  --output test_results/01_pairwise_alignment/runtime_summary_supplement.csv

python3 experiments/01_pairwise_alignment/summarize_geometry_superposition.py
python3 experiments/01_pairwise_alignment/plot_homstrad_runtime_deciles.py
python3 experiments/01_pairwise_alignment/plot_homstrad_runtime_deciles.py \
  --all-methods
```

Outputs are written to `test_results/01_pairwise_alignment/`. Reference
residue pairs are uppercase-uppercase manual-alignment columns. DALI failures
and no-hit cases contribute empty predictions to the all-pair accuracy macro.
Bootstrap intervals use 10,000 resamples and seed 42.

## 2. EAlign score ablation

The paper compares Full EAlign with four ablations:

| Setting | Removed component |
| --- | --- |
| `baseline` | none |
| `without-aa` | amino-acid evidence |
| `without-3di` | 3Di evidence |
| `without-vector` | orientation term |
| `without-distance` | distance term |

Build and run all five settings with one benchmark worker:

```bash
python3 experiments/02_score_ablation/ablation_ealign_score.py \
  --baseline --all --jobs 8 --workers 1
```

Reproduce the paper table and paired-bootstrap intervals from the five released
per-pair CSVs in `test_results/02_score_ablation/ablation_inputs/`:

```bash
python3 experiments/02_score_ablation/analyze_ealign_score_ablation.py \
  --bootstrap 10000 --seed 42
```

Analyze the MALISAM rescued failures and generate their correspondence plots:

```bash
python3 experiments/02_score_ablation/analyze_ealign_score_mechanisms.py
python3 experiments/02_score_ablation/plot_ealign_score_maximum_cases.py
```

Outputs are written to `test_results/02_score_ablation/`. To analyze a new
full rerun instead of the released CSVs, pass the generated files through
`--full-csv` and repeatable `--variant LABEL=CSV` arguments.

## 3. SABmark retrieval

Extracting `SABmark.zip` into `dataset/` produces the following layout:

```text
dataset/SABmark/structures/
dataset/SABmark/supfp/sup_fp.jsonl
dataset/SABmark/twifp/twi_fp.jsonl
```

The runners and evaluator derive TT/TF records directly from the two JSONL
files; no precomputed manifest is required. FF records are skipped. The paper
uses 16 outer workers for all five methods. GTalign-CPU additionally uses two
internal CPU threads per process; Foldseek uses one internal thread per outer
worker.

Run SupFP and TwiFP:

```bash
for subset in supfp twifp; do
  if [ "$subset" = "supfp" ]; then
    input="dataset/SABmark/supfp/sup_fp.jsonl"
  else
    input="dataset/SABmark/twifp/twi_fp.jsonl"
  fi

  structures="dataset/SABmark/structures"
  results="test_results/03_sabmark_retrieval/results/$subset"

  python3 experiments/03_sabmark_retrieval/run_twi_sup_fp_ealign.py \
    "$subset" \
    --input-jsonl "$input" \
    --structures "$structures" \
    --binary binary/ealign \
    --workers 16 \
    --output-root "$results/ealign"

  python3 experiments/03_sabmark_retrieval/run_twi_sup_fp_deepalign.py \
    "$subset" \
    --input-jsonl "$input" \
    --structures "$structures" \
    --binary binary/baselines/DeepAlign/DeepAlign \
    --workers 16 \
    --output-root "$results/deepalign"

  python3 experiments/03_sabmark_retrieval/run_twi_sup_fp_tmalign.py \
    "$subset" \
    --input-jsonl "$input" \
    --structures "$structures" \
    --binary binary/baselines/TMalign/TMalign \
    --workers 16 \
    --output-root "$results/tmalign"

  python3 experiments/03_sabmark_retrieval/run_twi_sup_fp_foldseek.py \
    "$subset" \
    --input-jsonl "$input" \
    --structures "$structures" \
    --binary binary/baselines/foldseek/bin/foldseek \
    --workers 16 \
    --output-root "$results/foldseek"

  python3 experiments/03_sabmark_retrieval/run_twi_sup_fp_gtalign.py \
    "$subset" \
    --input-jsonl "$input" \
    --structures "$structures" \
    --backend cpu \
    --binary binary/baselines/gtalign_alpha/patched_cpu/bin/gtalign \
    --workers 16 \
    --cpu-threads 2 \
    --output-root "$results/gtalign_cpu_patched_t2_w16"
done
```

These commands write the summary files expected by the evaluator under
`test_results/03_sabmark_retrieval/results/{supfp,twifp}/`.

Compute the query-macro AUROC, mean average precision, and ROC1:

```bash
python3 experiments/03_sabmark_retrieval/evaluate_record_metrics.py supfp
python3 experiments/03_sabmark_retrieval/evaluate_record_metrics.py twifp
```

Outputs are written to `test_results/03_sabmark_retrieval/evaluation/`.
EAlign and DeepAlign scores are divided by the shorter structure length;
TM-align and GTalign use shorter-length-normalized TM-score. Foldseek is
evaluated with both bit score and probability. The seven single-class SupFP
queries are excluded consistently from all metrics.

## Useful options

- `--limit N`: run a small pairwise or ablation smoke test.
- `--overwrite`: rerun jobs instead of reusing compatible saved results.
- `--binary PATH`: override a default executable.
- `--allow-incomplete`: preview a partial runtime result; do not use it for the
  reported paper table.
