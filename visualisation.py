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

class Visualisation():
    """ Class used to generate all the plots to analyse the results: the factors
    matrix heat map which displays all stable factors, the factors composition
    pie charts which display the variance caotured within each modality by each
    factor, and the detailed factors plots which displays the top features in
    each factor in each relevant modality.

    Parameters
    ----------
    plots_path : str
        Path to directory where plots should be saved.
    """
    def __init__(self, plots_path: str):
        self.plots_path = plots_path
        self.create_plots_dir()


    def create_plots_dir(self):
        """ Creates plots directory if it doesn't exist.
        """
        if not os.path.exists(self.plots_path):
            os.makedirs(self.plots_path)
    

    def get_colors(self, values: np.ndarray, cmap_name: str ="yellow") -> dict:
        """ Generated a color mapping for a list of values that are displayed in
        a plot.

        Parameters
        ----------
        values : np.ndarray
            List of all values.
        cmap_name : str, optional
            Name of a type of color map, by default "yellow".

        Returns
        -------
        dict
            Color mapping for each value in the list.
        """
        # Generate colors where negative values use a pink/yellow gradient and positive values use a green/blue gradient.
        if cmap_name == "yellow":
            cmap = sns.diverging_palette(40, 240, s=100, l=50, as_cmap=True)
        else:
            cmap = sns.color_palette(cmap_name, as_cmap=True)
        
        # Both positive & negative values (Use entire colormap)
        if values.min() < 0 and values.max() > 0:
            norm = mcolors.TwoSlopeNorm(vmin=values.min(), vcenter=0, vmax=values.max())
            colors = {v: cmap(norm(v)) for v in values}

        # Only positive values 
        elif values.min() >= 0:
            green_cmap = cmap(np.linspace(0.5, 1, 256)) 
            norm = mcolors.Normalize(vmin=values.min(), vmax=values.max())  
            colors = {v: green_cmap[int(norm(v) * 255)] for v in values}

        # Only negative values
        elif values.max() <= 0:
            pink_cmap = cmap(np.linspace(0, 0.5, 256))
            norm = mcolors.Normalize(vmin=values.min(), vmax=values.max())
            colors = {v: pink_cmap[int(norm(v) * 255)] for v in values}

        return colors


    def generate_color_map(self, Ws: list[np.ndarray]) -> tuple:
        """ Based on the values in a list of factor weights, generated a gradient 
        color map that will be used in a plot. 

        Parameters
        ----------
        Ws : list[np.ndarray]
            List of factor weights that will be displayed in plot. 

        Returns
        -------
        tuple
            Color map information needed for plot.
        """
        # Generate color for each value
        sorted_Ws = np.sort(np.concatenate(Ws).flatten())
        color_mapping = self.get_colors(sorted_Ws)
        
        # Extract values and colors
        values = list(color_mapping.keys())
        colors = list(color_mapping.values())

        # Create color map object
        cmap = mcolors.ListedColormap(colors)
        norm = mcolors.BoundaryNorm(values + [max(values) + 1], cmap.N)

        return cmap, norm, color_mapping


    def get_factors_matrix_plot_params(self, factors_Ws: list[np.ndarray], d: list[int], n_features: int) -> dict:
        """ Compiles and returns all parameters needed to plot the factors 
        matrix.

        Parameters
        ----------
        factors_Ws : list[np.ndarray]
            List of all stable factors weights identified by all models. 
            Dimension should be n_stable_factors * n_models * n_features, where
            n_stable_factors is the number of stable factors,
            n_models is the number of models that identified the stable factor,
            and n_features is the total number of features (n_features = sum(d)).
        d : list[int]
            List of dimensionalities of each modality
        n_features : int
            Total number of features in analysis.

        Returns
        -------
        dict
            Structure with factors matrix plot parameters.
        """
        # print(factors_Ws)
        print('LALA')
        # a = min(arr.min() for arr in factors_Ws)
        a = np.array(factors_Ws).min()
        print('LALa2')
        b = np.array(factors_Ws).max()
        print('LALa3')
        c = [sum(d[:i]) for i in range(1, len(d))]
        plot_params = {
            "vmin": np.array(factors_Ws).min(),
            "vmax": np.array(factors_Ws).max(),
            "axes_sep": [sum(d[:i]) for i in range(1, len(d))],
            "label_size": 10,
            "top_margin": 0.04,
            "bottom_margin": 0.04
        }
        print('LALa4')
        dpi = 72.27
        matrix_height_pt = plot_params["label_size"] * n_features
        matrix_height_in = matrix_height_pt / dpi
        plot_params["figure_height"] = matrix_height_in / (1 - plot_params["top_margin"] - plot_params["bottom_margin"])

        print('LALa5')
        return plot_params


    def plot_stable_factors_matrix(self, factors_Ws: list[np.ndarray], d: list[int], feature_labels: list[str], save: bool = False):
        """ Function that plots the stable factors matrix - a heatmap
        for each stable factor, where the weights are displayed for all
        "versions" of the stable factor, i.e. from all models that identified
        the factor. Shows the similarities between the different versions of
        the stable factor, the differences between the factors, and the most
        relevant features and modalities within each factor.

        Parameters
        ----------
        factors_Ws : list[np.ndarray]
            List of all stable factors weights identified by all models. 
            Dimension should be n_stable_factors * n_models * n_features, where
            n_stable_factors is the number of stable factors,
            n_models is the number of models that identified the stable factor,
            and n_features is the total number of features (n_features = sum(d)).
        d : list[int]
            List of dimensionalities of each modality
        feature_labels : list[str]
            List of feature labels to be displayed in plot.
        save : bool, optional
            Bool value that indicates whether to save the plot or not, by default 
            False.
        """
        # Generate color map
        cmap, norm, _ = self.generate_color_map(factors_Ws)

        # Define plot parameters
        plot_params = self.get_factors_matrix_plot_params(factors_Ws, d, len(feature_labels))
        print('Ola')
        plt.rcParams["ytick.labelsize"] = plot_params["label_size"]

        print('Ola2')
        # Plot matrix
        print(f"len(factors_Ws): {len(factors_Ws)}")
        print(f"Figure height: {plot_params['figure_height']}")
        print(f"Top margin: {plot_params['top_margin']}")
        print(f"Bottom margin: {plot_params['bottom_margin']}")

        plt.figure()
        plt.plot([0, 1], [0, 1])
        plt.savefig("test_plot.jpg")
        plt.close()
        print("Saved test figure.")

        fig, axes = plt.subplots(ncols=len(factors_Ws),
                                            figsize=(len(factors_Ws) * 4, plot_params["figure_height"]),
                                            gridspec_kw=dict(top=1 - plot_params["top_margin"], bottom=plot_params["bottom_margin"]))

        print("LALA KK1")
        for i, factor in enumerate(factors_Ws):
            print("LALA KK2", i)
            if i == 0:
                # First factor, display feature labels
                ax = sns.heatmap(factor.T, cmap=cmap, norm=norm,
                                 fmt="", ax=axes[i], cbar=False,
                                 vmin=plot_params["vmin"], vmax=plot_params["vmax"],
                                 yticklabels=feature_labels)
                ax.set_xticks(ticks=[], labels=[])
            elif i == len(factors_Ws) - 1:
                # Last factor, display color map
                ax = sns.heatmap(factor.T, cmap=cmap, norm=norm,
                                 fmt="", ax=axes[i], cbar=True,
                                 vmin=plot_params["vmin"], vmax=plot_params["vmax"])
                ax.set_yticks(ticks=[], labels=[])
                ax.set_xticks(ticks=[], labels=[])
            else:
                ax = sns.heatmap(factor.T, cmap=cmap, norm=norm,
                                 fmt="", ax=axes[i], cbar=False,
                                 vmin=plot_params["vmin"], vmax=plot_params["vmax"])
                ax.set_yticks(ticks=[], labels=[])
                ax.set_xticks(ticks=[], labels=[])

            # Plot horizontal line to separate modalities
            ax.hlines(plot_params["axes_sep"], *ax.get_xlim())

        print("LALA KK _ alnost")
        # Add plot and axes labels
        fig.supxlabel("Stable factor", fontsize=plot_params["label_size"] * 3)
        fig.supylabel("Feature", fontsize=plot_params["label_size"] * 3)
        fig.suptitle(f"Stable factors", fontsize=plot_params["label_size"] * 3)

        print("LALA KK")
        # Save plot
        plt.savefig("test_debug.jpg")
        if save:
            plt.savefig(f"{self.plots_path}/stable_factors_matrix.svg", format="svg", bbox_inches="tight")
            plt.savefig(f"{self.plots_path}/stable_factors_matrix.jpg", format="jpg", bbox_inches="tight")
            plt.close()
        else:
            plt.show()
    

    def plot_factors_composition(self, n_factors: int, variances: list[np.ndarray], modalities: list[str], cutoff: float = 0, save: bool = False):
        """ Plots pie chart for each factor which shows the amount of variance
        captured by the factor in each modality. Highlights the most significant
        modalities within each factor.

        Parameters
        ----------
        n_factors : int
            Number of stable factors.
        variances : list[np.ndarray]
            Factors variances, shape is n_factors * n_modalities.
        modalities : list[str]
            List of modalities display names.
        cutoff : float, optional
            Value which determines which of the modalities are pruned from the 
            factor plot, by default 0. All modalities whose variance is less or
            equal to the cutoff are removed from plot. 
        save : bool, optional
            Bool value that indicates whether to save the plot or not, by default 
            False.
        """
        colors = px.colors.qualitative.Set2
        modality_colors = dict(zip(modalities, colors))

        for i in range(n_factors):
            factor_df = pd.DataFrame(data={
                "Modality": modalities,
                "Variance": variances[i]
            })
            # Filter out modalities with variance <= cutoff
            factor_df = factor_df[factor_df["Variance"] > cutoff]

            # Plot pie chart
            fig = px.pie(factor_df, names="Modality", values="Variance",
                         color="Modality", color_discrete_map=modality_colors,
                         hole=0.3)
            
            # Add annotations to the chart, format "Modality [variance]"
            fig.update_traces(textinfo="none", texttemplate="%{label} [%{customdata:.2f}%]", 
                              customdata=factor_df["Variance"], textposition="outside",
                              pull=[0.05] * len(factor_df), textfont_size=18)
            
            # Add title
            fig.update_layout(title=dict(
                text=f"Stable factor {i + 1} composition",
                font=dict(size=22),
                y=0.98 ),
                legend_font_size=18, width=1000, height=1050)
            
            # Save plot
            if not save:
                fig.show()
            else:
                fig.write_image(f"{self.plots_path}/factor_{i+1}_composition.jpg", format="jpg", scale=3)
                fig.write_image(f"{self.plots_path}/factor_{i+1}_composition.svg", format="svg")
    
    
    def get_factors_features_plot_params(self, n_modalities: int, all_weights: list[np.ndarray]) -> dict:
        """ Compiles and returns all parameters needed to plot the detailed 
        factor.

        Parameters
        ----------
        n_modalities : list[np.ndarray]
            Number of modalities to be displayed in detailed plot of factors,
            which determines the number of subplots.
        all_weights : list[np.ndarray]
            List of factor weights that will be displayed in plot. 

        Returns
        -------
        dict
            Structure with detailed factors plot parameters.
        """
        params = {}
        params["cmap"], params["norm"], params["color_mapping"] = self.generate_color_map(all_weights)
        params["n_subplots"] = n_modalities
        params["n_cols"] = math.ceil(math.sqrt(params["n_subplots"]))
        params["n_rows"] = math.ceil(params["n_subplots"] / params["n_cols"])
        params["fig_width"] = 4 * params["n_cols"]
        params["fig_height"] = 8 * params["n_rows"]
        params["subplots_wspace"]=1.5
        params["subplots_hspace"]=0.2
        return params

    
    def filter_factor_weights(self, Ws: list[np.ndarray], variances: list[np.ndarray], 
                              modalities: list[str], top_n: int, cutoff: float, 
                              split_indeces: list[int], labels_mapping: dict, 
                              i: int) -> tuple:
        """ Filter factor modalities and weights which will be displayed in plot 
        based on cutoff and and the number of top features that have to be 
        displayed.

        Parameters
        ----------
        Ws : list[np.ndarray]
            Factors weights, shape is n_factors * n_features.
        variances : list[np.ndarray]
            Factors variances, shape is n_factors * n_modalities.
        modalities : list[str]
            List of modalities display names.
        top_n : int
            Number of top features to be added to the plot.
        cutoff : float, optional
            Value which determines which of the modalities are pruned from the 
            factor plot. All modalities whose variance is less or equal to the 
            cutoff are removed from plot. 
        split_indeces : list[int]
            Indeces in the factor weights array where each modality begins.
        labels_mapping : dict
            Mappings of each modality to its corresponding feature labels list.
        i : int
            Index of stable latent factor.

        Returns
        -------
        tuple
            The filtered weights of the factor split for each modality and a 
            list of all the weghts which will be displayed in plot (needed
            for color mapping).
        """
        factor_df = pd.DataFrame(data={
                "Modality": modalities,
                "Variance": variances[i]
            })
        # Filter out modalities with variance <= cutoff
        factor_df = factor_df[factor_df["Variance"] > cutoff]
        plot_modalities = factor_df["Modality"].to_list()

        # Get selected modality weights
        split_W = np.split(Ws[i], split_indeces)
        split_W = {modality: split_W[j] for j, modality in 
                   enumerate(modalities) if modality in plot_modalities}
            
        # Select top modality weights
        all_weights = []
        for j, modality in enumerate(split_W):
            modality_weights = split_W[modality]
            modality_labels = labels_mapping[modality]
            top_features = sorted(zip(modality_weights, modality_labels), 
                                  key=lambda x: abs(x[0]), reverse=True)[:top_n]
            top_weights, top_labels = zip(*sorted(top_features, key=lambda x: x[0], reverse=True))
            all_weights.append(top_weights)
            split_W[modality] = {
                    "weights": np.array(list(top_weights)),
                    "labels": np.array(list(top_labels))
                }
            
        return split_W, all_weights  


    def plot_factors_features(self, n_factors: int, Ws: list[np.ndarray], 
                              variances: list[np.ndarray], modalities: list[str], 
                              d: list[int], feature_labels: list[str], 
                              top_n: int, cutoff: float = 0, save: bool = False):
        """ For each stable factor, displays the top_n features from each
        relevant modality of each stable factor. 

        Parameters
        ----------
        n_factors : int
            Number of stable factors.
        Ws : list[np.ndarray]
            Factors weights, shape is n_factors * n_features.
        variances : list[np.ndarray]
            Factors variances, shape is n_factors * n_modalities.
        modalities : list[str]
            List of modalities display names.
        d : list[int]
            List of dimensionality of each modality.
        feature_labels : list[str]
            List of feature names which will be displayed in plots.
        top_n : int
            Number of top features to be added to the plot.
        cutoff : float, optional
            Value which determines which of the modalities are pruned from the 
            factor plot, by default 0. All modalities whose variance is less or
            equal to the cutoff are removed from plot. 
        save : bool, optional
            Bool value that indicates whether to save the plot or not, by default 
            False.
        """
        split_indeces = [sum(d[:i+1]) for i in range(len(d) - 1)]
        labels_mapping = np.split(feature_labels, split_indeces)
        labels_mapping = {modalities[i]: labels for i, labels in enumerate(labels_mapping)}
        
        for i in range(n_factors):
            split_W, all_weights = self.filter_factor_weights(Ws, variances, modalities, top_n, cutoff, split_indeces, labels_mapping, i)
            plot_params = self.get_factors_features_plot_params(len(split_W), all_weights)

            fig = plt.figure(figsize=(plot_params["fig_width"], plot_params["fig_height"]))
            gs = gridspec.GridSpec(plot_params["n_rows"], plot_params["n_cols"], figure=fig, wspace=plot_params["subplots_wspace"], hspace=plot_params["subplots_hspace"])  # Adjust spacing

            for j, modality in enumerate(split_W):
                # Bar plot
                row = j // plot_params["n_cols"]
                col = j % plot_params["n_cols"]
                ax = fig.add_subplot(gs[row, col])
                sns.barplot(x=split_W[modality]["weights"], y=split_W[modality]["labels"],
                            ax=ax, orient="h", errorbar=("ci", False), dodge=False,
                            palette=plot_params["color_mapping"], hue=split_W[modality]["weights"])
                ax.set_title(f"{modality}")
                ax.spines[["right", "top"]].set_visible(False)
                ax.legend_.remove()

                # Add x ticks for min and max values
                min_value = min(split_W[modality]["weights"]) if min(split_W[modality]["weights"]) < 0 else 0
                max_value = max(split_W[modality]["weights"]) if max(split_W[modality]["weights"]) > 0 else 0
                t = [f"{min_value:.2e}" if min_value != 0 else "0", f"{max_value:.2e}" if max_value != 0 else "0"]
                ax.set_xticks([min_value, max_value], t)
                ax.axvline(0, color="black", linewidth=0.8)

            # Add plot titple
            suptitle = fig.suptitle(f"Stable factor {i + 1} top features")
            fdct = {"fontsize": 20}
            suptitle.set(**fdct)

            # Add plot color bar
            sm = matplotlib.cm.ScalarMappable(cmap=plot_params["cmap"], norm=plot_params["norm"])
            fig.colorbar(sm, ax=fig.axes, orientation="horizontal", fraction=0.02, pad=0.08)

            # Save plot
            if not save:
                plt.show()
            else:
                plt.savefig(f"{self.plots_path}/factor_{i+1}_top_features.svg", format="svg", bbox_inches="tight")
                plt.savefig(f"{self.plots_path}/factor_{i+1}_top_features.jpg", format="jpg", bbox_inches="tight")
                plt.close()
