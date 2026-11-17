# Nextflow pipeline: methfast by sample

Runs `methfast <methylation_bed(.gz)> <target_bed> -o <output_file>` independently for each WGBS sample in a metadata TSV.

## Inputs

The TSV must have a header called `wgbs_directory_name` (other columns are ignored):

```tsv
wgbs_directory_name	group
273919-dx-wgbs	IDHmut
514066-dx-wgbs	IDHwt
```

Each ID resolves to:

```
/storage2/fs1/dspencer/Active/spencerlab/data/wgbs/<id>/<id>.meth.bed.gz
```

## Run on SLURM

Run from the directory containing `main.nf` and `nextflow.config`:

```bash
nextflow run main.nf \
    --samples /path/to/samples.tsv \
    --target_bed /path/to/IDH_DMRs.bed \
    --outdir /path/to/results/DMR_methylation \
    -profile slurm \
    -resume
```

Optional parameters:

```bash
    --postfix custom_label \
    --wgbs_root /another/wgbs/root
```

Default postfix: the target BED basename with `.bed` or `.bed.gz` removed. For `IDH_DMRs.bed`, the output is:

```
results/DMR_methylation/
├── 273919-dx-wgbs/
│   └── 273919-dx-wgbs.IDH_DMRs.bed
└── 514066-dx-wgbs/
    └── 514066-dx-wgbs.IDH_DMRs.bed
```

`--samples`, `--target_bed`, and `--outdir` are required. Duplicate IDs are run once. Missing methylation inputs, empty IDs, or invalid postfixes stop the run with an error.

## Runtime and resources

- Requires Nextflow, Java, and `methfast` executable on the SLURM compute nodes.
- Defaults per sample: 4 CPUs, 16 GB RAM, 3 hours (adjust if `methfast` is single-threaded). Edit `nextflow.config` for your cluster and job sizes.
- The `slurm` profile dispatches independent jobs to `condo-dspencer` under account `compute2-dspencer`. The `local` profile can be used for testing.
- Up to 15 tasks may be queued concurrently, with submissions throttled to one per second.
- `report.html`, `timeline.html`, and `trace.txt` are generated in the launch directory; subsequent runs overwrite them.
- The Nextflow launcher must remain running (e.g., in `tmux` or a persistent SLURM allocation). Nextflow submits and monitors the compute jobs.
- Use `-resume` to reuse successful tasks; retain the `work/` directory and `.nextflow/` cache from the original launch directory.
- The original `.meth.bed.gz` files are staged as links into task work directories rather than copied.

Micromamba support is enabled, but **this alone does not create or activate an environment**. The process currently assumes `methfast` is accessible on compute nodes. If it is installed in a specific environment, add a `conda` directive for `METHFAST` pointing to that environment or an environment YAML. `errorStrategy = 'terminate'` is intentional: a failed methylation job should not be silently reported as overall success.

For debugging, run with `-with-trace trace.tsv` and inspect `work/<hash>/.command.err` for failed tasks. To see task commands without submitting SLURM jobs, use `-profile local` with a small test sample (this **executes** the command; it is not a dry run).
