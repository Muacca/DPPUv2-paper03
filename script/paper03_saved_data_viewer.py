"""Notebook widget for saved Paper03 geometric slices. No recomputation."""
from pathlib import Path
import numpy as np

def show(data_dir=None):
    import ipywidgets as widgets
    import matplotlib.pyplot as plt
    from IPython.display import display,clear_output
    from paper_results.data import load_data
    source=Path(data_dir) if data_dir else Path(__file__).resolve().parents[1]/"data/revised"
    data,_=load_data(source)
    options={"Radius: MX/AX/VT":("radius.csv","r",["MX","AX","VT"]),
             "Shear: exponential paths":("shear.csv","epsilon",["H0","H2","H1"]),
             "Twist and rotated frame":("twist.csv","omega1",["U","U_rotated_frame"]),
             "Constrained shape minimum":("shape_minima.csv","epsilon",["F_axis","F_min","F_global"]),
             "Quartic twist term":("twist_kinetic.csv","rho",["static","rotating_derivative"]),
             "Velocity cross coefficients":("cross_kinetic.csv","w",["G_rw","G_sw","G_rs"]),
             "Coupled periodic action":("coupled_action.csv","n",["delta_static","derivative","delta_total"]),
             "Tensor expansion":("invariants.csv","q1",["delta_U","quartic"])}
    menu=widgets.Dropdown(options=list(options),description="Slice:")
    output=widgets.Output()
    def draw(change=None):
        name,x,ys=options[menu.value];d=data[name]
        with output:
            clear_output(wait=True)
            fig,ax=plt.subplots(figsize=(7,4))
            for y in ys:ax.plot(d[x],d[y],label=y)
            ax.set_xlabel(x);ax.set_ylabel("Geometric quantity")
            ax.legend();fig.tight_layout();display(fig);plt.close(fig)
    menu.observe(draw,names="value")
    display(widgets.VBox([widgets.HTML("<b>Paper03 geometric sectors — saved data</b>"),menu,output]))
    draw()
    return menu,output
