import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import plotly.express as px
import matplotlib.colors as mcolors
import matplotlib
import numpy as np
import pandas as pd
import math

def create_plots_dir(plots_path):
    """ Creates plots directory if it doesn't exist.
    """
    if not os.path.exists(plots_path):
        os.makedirs(plots_path)

def add(x, y):
    return x + y

def square(x):
    return x ** 2