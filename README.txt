cfg - Folder containing configuration files
documentation - Folder containing files used to construct documentation
output - Folder containing raw output for all experiments performed during testing

If you do not have access to a built version of the documentation, we suggest
building one by going to the documentation/ folder and doing:
make html
This will build an html documentation setup accessible in documentation/_build/html/index.html

To run experiments, use main.py
To analyze and plot experimental data, use plotter.py
To recreate the predicted amounts of wasted evaluations plot from the paper, run wasteplot.py
