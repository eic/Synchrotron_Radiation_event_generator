# Synchrotron-radiation event generator from Synrad+ table

authors: Reynier Cruz-Torres^, Benjamen Sterwerf^^

-^ Lawrence Berkeley National Laboratory

-^^ University of California, Berkeley

-----

1. Download csv file stored [here](https://drive.google.com/file/d/1XX78_qeuoMK8xhuOB5QgbUyye7Lv_xPg/view?usp=sharing).

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

Libraries needed to run this notebook: ```yaml```, ```numpy```, ```matplotlib.pyplot```, ```pandas```, ```pyhepmc```, ```ROOT```, ```os```, ```argparse``` and ```time```.













