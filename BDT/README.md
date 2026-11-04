## How to run

Start from the python env:

```bash
cd /home/users/garciaja/fullRun3/CMSSW_15_0_2/src/run3_scouting
source setup.sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
EV=BDT/PNN_output/parametrized/dqcd_mratio/eval_testonly_prevlim
DATA="/ceph/cms/store/group/Run3Scouting/tuples_data_nhfix/tuples_Data*_2024_*.root"
```

Train the models. Skip this if you already have them (they are in PNN_output/parametrized/dqcd_mratio/models/binned_lxy):

```bash
python3 BDT/train_discriminant.py --model pnn --signal dqcd --conditional
```

Then the eval. This is the slow one:

```bash
python3 BDT/evaluate_workingpoint.py --combine-jobs 32 --model dqcd --model-dir BDT/PNN_output/parametrized/dqcd_mratio --conditional --cond-var mratio --sig-xsec-hepdata previous_limits --eval-suffix prevlim
```

Score the data at both working points:

```bash
python3 BDT/data_at_wp.py $EV --wp tight --rescore --write-scored --data "$DATA"
python3 BDT/data_at_wp.py $EV --wp loose --rescore --write-scored --data "$DATA"
```

Score the signal. You only need this once, loose reuses the same files:

```bash
python3 BDT/make_fit_inputs.py $EV --stage score
```

Now switch to cmsenv and write the RooDataSets:

```bash
deactivate
cd /home/users/garciaja/fullRun3/CMSSW_15_0_2/src
cmsenv
cd run3_scouting

python3 BDT/make_fit_inputs.py $EV --stage write --wp tight
python3 BDT/make_fit_inputs.py $EV --stage write --wp loose
```

The RooDataSets end up in `fitInputs/2024_pnn/dqcd_mratio_WP_tight/` and `fitInputs/2024_pnn/dqcd_mratio_WP_loose/`, with one folder per mass point.

## Remaking the plots

```bash
python3 BDT/plot_eval.py $EV
```

With no argument it goes through every eval folder under `BDT/*_output/`.

## Good to know

- `$DATA` has to be passed explicitly. The default in `data_at_wp.py` is `tuples_data/`, not the `tuples_data_nhfix/` tuples.
- Right now the data is eras C, D, E, H and I, 43.43 fb⁻¹. To add an era, update `ERA_LUMI_FB` in both `data_at_wp.py` and `make_fit_inputs.py`. `data_at_wp.py` picks up new tuples by itself, but the signal normalization only changes if you edit those two tables.
- The loose `--stage write` fails for 9 short ctau points, because there is no signal in a displaced bin to set a threshold with. Those files just don't get written.
- The signal in the RooDataSets is normalized to 1 pb.
- The resonance veto is always on, and only applies to background.
- `combine` runs from the python env through a small wrapper in `~/bin` that sets up its own CMSSW. If you're running this on another account you'll need one too.
