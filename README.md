# Synchrotron-radiation event generator from Synrad+ table

authors: Reynier Cruz-Torres^, Benjamen Sterwerf^^

-^ Lawrence Berkeley National Laboratory

-^^ University of California, Berkeley

-----

1. Download csv file stored [here](https://drive.google.com/file/d/1XX78_qeuoMK8xhuOB5QgbUyye7Lv_xPg/view?usp=sharing).

2. Create a yaml configuration file with the following information:

	- ```input_single_photons```: path to csv file downloaded in step 1.
	- ```n_events```: number of events to be generated.
	- ```integration_window```: time window that will define one event.
	- ```seed```: random seed for reproducibility. Set to 0 to leave the seed unconstrained.

3. Open the Jupyter notebook ```sr_event_generator.ipynb```, edit the first cell to point to the yaml configuration file created in step 2, and run the notebook.

-----

Libraries needed to run this notebook: ```yaml```, ```numpy```, ```matplotlib.pyplot```, ```pandas```, ```pyhepmc_ng```, ```ROOT```, ```os``` and ```time```.













