## Boltzmann Reweighting in Machine Learning for Tailoring Pheromone Molecule Design to Insect Odorant Receptors

# How to Reweighting

First make sure you have all dependencies installed by running `pip install -r requirements.txt`.

Then, you must prepare a Excel include smiles and docking score like example.xlsx. You can easily get the reweighting dataset (reweighting.csv) as follow.

```shell
python reweight.py --file example.xlsx
```

# How to train
We thank the [previous work by the SeqGAN team and ORGAN team](https://github.com/LantaoYu/SeqGAN, https://github.com/gablg1/ORGAN). ORGAN is been utilized to train our model.

We provide a working example that can be run with :

```shell
python organ.py
```


