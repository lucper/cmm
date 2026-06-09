## Accuracy experiment

Steps:
- Make sure the directory `../string_db_cleaned` (path relative to the current directory) exists with files `<dataset>.s700.cleaned.fa` and `<dataset>.s700.cleaned.int` therein. They should have been gone through the `clean_dataset.py` script.
- Make sure SLIDER is installed: the jar file should be available at `../competitors/SliderLight/dist/SliderLight.jar` (path relative to the current directory).
- Create the directory `cmm_out` wih outputs from `cmm`. Files should be names `<dataset>.cmm.out`.
- Run `make` to compile `cmm_eval` and `cmm_dedup`.
- Run `run.sh`: it will run SLIDER and produce the directory `slider_out`.
- Run `make_plots.sh`: it will use the data in `slider_out` and produce the directory `plots`.
