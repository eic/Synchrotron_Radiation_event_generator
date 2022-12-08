# Synchrotron-radiation event generator from Synrad+ table

authors: Reynier Cruz-Torres^, Benjamen Sterwerf^^

-^ Lawrence Berkeley National Laboratory

-^^ University of California, Berkeley

-----

1. Download csv file stored [here](https://drive.google.com/file/d/1XX78_qeuoMK8xhuOB5QgbUyye7Lv_xPg/view?usp=sharing). You can get this file following one of the two methods below:

```bash
wget -O combined_data.csv 'https://drive.google.com/uc?export=download&id=1XX78_qeuoMK8xhuOB5QgbUyye7Lv_xPg&confirm=no_antivirus'
```

or

```bash
curl -L 'https://drive.google.com/uc?export=download&id=1XX78_qeuoMK8xhuOB5QgbUyye7Lv_xPg&confirm=no_antivirus' > combined_data.csv
```

2. Create a yaml configuration file (e.g. ```config.yaml```) with the following information:

	- ```input_single_photons```: path to csv file downloaded in step 1.
	- ```n_events```: number of events to be generated.
	- ```integration_window```: time window that will define one event.
	- ```seed```: random seed for reproducibility. Set to 0 to leave the seed unconstrained.

3. Run the generator as: 

```bash
python3 sr_generator.py --configFile config.yaml
```

-----

Libraries needed to run this notebook: ```pyyaml```, ```numpy```, ```matplotlib```, ```pandas```, ```pyhepmc```, ```ROOT```, ```os```, ```argparse``` and ```time```.

If you have trouble getting pyroot to work, [this](https://root-forum.cern.ch/t/cannot-import-root-6-22-in-python-3-9-1-could-not-load-cppyy-backend-library/43764/8) may be useful.













