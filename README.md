# BinAlloy 

This module allows one to calculate temperature profiles of the long-range order, energy and short-range order in a binary alloy. The back-end is fully written in C for optimal speed.

The current version will only let one execute the algortihm for CuZn.
In the next version, the user will be able to pick metals from a list and even supply custom bond energies.

# Installation
The project is not yet listed on pipy
To install the module to your computer, download the .whl file from the ```build``` folder and execute
```
pip install $PATH/BinAlloy-0.1-cp39-cp39-macosx_10_9_x86_64.whl
```

You can now use the module like you would any other python package.
To understand better what the functionality offers, have a look at the following example:

https://github.com/Thomas-Jaeken/BinAlloy/blob/4e402cfcd7568829ad3ce771b6328bb2bb7b1bf6/example.py#L1-L39
