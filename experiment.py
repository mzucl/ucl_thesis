from .dataset import Dataset
from .gfa.models import GFA_DiagonalNoiseModel
from .gfa.utils import GFAtools
import os
import numpy as np
import time
import pickle
from sklearn.metrics import roc_curve, balanced_accuracy_score, roc_auc_score, precision_score, recall_score, f1_score, confusion_matrix
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from sklearn.cluster import DBSCAN
from collections import Counter
from sklearn.metrics.pairwise import cosine_similarity


class Experiment():
    """ Class that stores all information needed for experiments, perform 
    experiments and classification with GFA, computes classification results and 
    factors information resulted from analysis.

    Parameters
    ----------
    n_folds : int
        Number of folds.
    n_reps : int
        Number of training repetitions to be performed for each train-test split.
    dataset : Dataset
        Dataset object with the data used in experiment.
    models_path : str
        Path to folder where to store and read the trained GFA models.
    factors_path : str
        Path to folder where to store and read the GFA factors.
    taus_path : str
        Path to folder where to store and read the GFA taus.
    
    Attributes
    ----------
    n_folds : int
        Number of folds.
    n_reps : int
        Number of training repetitions to be performed for each train-test split.
    dataset : Dataset
        Dataset object with the data used in experiment.
    models_path : str
        Path to folder where to store and read the trained GFA models.
    factors_path : str
        Path to folder where to store and read the GFA factors.
    taus_path : str
        Path to folder where to store and read the GFA taus.
    Ks : list[int]
        Number of factors identified by each model.
    Ws : list[np.ndarray]
        List of factors matrices from each model.
    taus : list
        List of taus from each model.
    models : list[GFA_DiagonalNoiseModel]
        List of all trained models.
    elbos : list[float]
        List of ELBO from last training iteration from each model.
    y_pred : list[dict]
        List of prediction results in continuous format for each train-test 
        split.
    y_pred_binary : list[dict]
        List of prediction results in binary format for each train-test split.
    results : pandas.DataFrame
        Dataframe with performance for all models.
    total_vars : list[np.float64]
        List of total variance captured by each model.
    factors_vars : list[np.float64]
        List of total variance captured by the factors of each model.
    vars_within : list[np.ndarray]
        List of variances captured by each factor within each modality for all 
        models.
    rel_vars_within : list[np.ndarray]
        List of relative varainces  captured by each factor within each modality
        for all models.
    factors : list[np.ndarray]
        List of all factors identified by all models.
    factors_cluster_labels : np.ndarray
        List of cluster labels of all factors.
    stable_factors : list[dict]
        List of information for each stable factor identified.
    """
    def __init__(self, n_folds: int, n_reps: int, dataset: Dataset, models_path: str, factors_path: str, taus_path: str):
        self.n_folds = n_folds
        self.n_reps = n_reps
        self.dataset = dataset
        self.models_path = models_path
        self.factors_path = factors_path
        self.taus_path = taus_path
        self.models_load_check = False
        self.factors_load_check = False
        self.taus_load_check = False
        self.prediction_check = False
        self.factors_variances_check = False
        self.factors_clustered_check = False


    def save_factors(self, model: GFA_DiagonalNoiseModel, fold: int):
        """ Extracts and saves the factors inferred by the GFA model.

        Parameters
        ----------
        model : GFA_DiagonalNoiseModel
            Trained GFA model from which the factors are extracted.
        fold : int
            The number of the fold on which the model was run, used to identify
            the file in which to save the factors.
        """
        # Created factors matrix
        W = np.zeros((np.sum(model.d), model.k))
        d = 0
        # Fill factors matrix
        for m in range(len(model.d)):
            Dm = model.d[m]
            W[d:(d + Dm), :] = model.means_w[m]
            d += Dm
        
        # Check factors directory exists
        if not os.path.exists(f"{self.factors_path}/"):
            os.makedirs(f"{self.factors_path}/")
        
        # Save factors in numpy file
        with open(f"{self.factors_path}/W_{fold}.npy", "wb") as f:
            np.save(f, W)
    

    def load_factors(self):
        """ Loads factors inferred by GFA for each train-test split from files. 
        Stores the factors in object (Ws) and the number of factors identified 
        by each model (Ks).
        """
        self.Ks = []
        self.Ws = []

        for fold in range(self.n_folds):
            W = np.load(f"{self.factors_path}/W_{fold}.npy")
            self.Ws.append(W)
            self.Ks.append(W.shape[1])
        
        self.factors_load_check = True

        
    def save_taus(self, model: GFA_DiagonalNoiseModel, fold: int):
        """ Extracts and saves the taus inferred by the GFA model.

        Parameters
        ----------
        model : GFA_DiagonalNoiseModel
            Trained GFA model from which the taus are extracted.
        fold : int
            The number of the fold on which the model was run, used to identify
            the file in which to save the taus.
        """
        # Extact taus and reshape
        E_taus = [tau.reshape(-1) for tau in model.E_tau]
        E_taus_array = np.hstack(E_taus)

        # Check taus directory exists
        if not os.path.exists(f"{self.taus_path}/"):
            os.makedirs(f"{self.taus_path}/")
        
        # Save taus to numpy file
        with open(f"{self.taus_path}/tau_{fold}.npy", "wb") as f:
            np.save(f, E_taus_array)
    

    def load_taus(self):
        """ Loads taus inferred by GFA for each train-test split from files,
        reshapes the taus to original format, and stores them in object (taus).
        """
        self.taus = []
        split_indeces = [sum(self.dataset.d[:i]) for i, _ in enumerate(self.dataset.d) if i > 0]
        for fold in range(self.n_folds):
            # Read taus
            tau = np.load(f"{self.taus_path}/tau_{fold}.npy")
            # Reshape taus to original shape
            tau = np.split(tau, split_indeces)
            tau = [t.reshape(1, -1) for t in tau]
            self.taus.append(tau)
        
        self.taus_load_check = True
    

    def save_model(self, model: GFA_DiagonalNoiseModel, fold: int):
        """ Saves the trained GFA model to file.

        Parameters
        ----------
        model : GFA_DiagonalNoiseModel
            Trained GFA model.
        fold : int
            The number of the fold on which the model was run, used to identify
            the file in which to save the model.
        """
        if not os.path.exists(f"{self.models_path}/"):
            os.makedirs(f"{self.models_path}/")
        
        # Save model to file
        file = open(f"{self.models_path}/fold_{fold}.txt", "wb")
        pickle.dump(model, file)
        file.close()

    
    def load_models(self):
        """ Loads trained GFA models and stores them within object. Extracts
        the ELBOS from last iteration for each model to rank them.
        """
        self.models = []
        self.elbos = []

        for fold in range(self.n_folds):
            file = open(f"{self.models_path}/fold_{fold}.txt", "rb")
            model = pickle.load(file)
            self.models.append(model)
            self.elbos.append(model.L[-1])
            file.close()
        
        self.models_load_check = True


    def fit_model(self, x: list[np.ndarray], model_params: dict) -> GFA_DiagonalNoiseModel:
        """ Fits one GFA instance on data and returns the model.

        Parameters
        ----------
        x : list[np.ndarray]
            Data to fit the model on, list of np.ndarrays representing data for 
            each modality.
        model_params : dict
            Dictionary of model parameters required by GFA_DiagonalNoiseModel 
            class.

        Returns
        -------
        GFA_DiagonalNoiseModel
            GFA model trained on data.
        """
        model = GFA_DiagonalNoiseModel(x, model_params)
        time_start = time.process_time()
        model.fit(x)
        model.time_elapsed = time.process_time() - time_start
        print(f"Computational time: {float('{:.2f}'.format(model.time_elapsed))}s")
        return model


    def run_fold_experiment(self, fold: int, model_params: dict, x: list[np.ndarray]):
        """ Trains n_reps different instances of GFA in parallel on data from
        one train-test split. Uses at most 5 parallel workers. Saves the best 
        model, its inferred factors and taus.

        Parameters
        ----------
        fold : int
            The number of the fold on which the model is trained, used to 
            identify which files to save the models, factors and taus in.
        model_params : dict
            Dictionary of model parameters required by GFA_DiagonalNoiseModel 
            class.
        x : list[np.ndarray]
            Data to fit the model on, list of np.ndarrays representing data for 
            each modality.
        """
        with ProcessPoolExecutor(max_workers=5) as executor:
            futures = [executor.submit(self.fit_model, x, model_params) for _ in range(self.n_reps)]
            models = [f.result() for f in futures]
        
        # Save best model
        best_model_index = np.argmax([model.L[-1] for model in models])
        self.save_model(models[best_model_index], fold)
        self.save_factors(models[best_model_index], fold)
        self.save_taus(models[best_model_index], fold)


    def run_experiment(self, K: int):
        """ Runs experiment for all train-test splits.

        Parameters
        ----------
        K : int
            Number of factors to initialise GFA model with.

        Raises
        ------
        ValueError
            If data was not fully prepared for experiemnt (read -> split -> 
            scaled).
        """
        if not self.dataset.data_scaled_check:
            raise ValueError("Data is not ready for experiments.")
        
        model_params = {
            "num_groups": len(self.dataset.d),
            "K": K,
            "scenario": "incomplete"
        }

        for fold in range(self.n_folds):
            # Get data for fold
            train_data = self.dataset.scaled_data[fold]["train"]
            test_data = self.dataset.scaled_data[fold]["test"]
            x = [np.vstack((train_data[i], test_data[i])) for i in range(len(self.dataset.d))]

            # Run experiment for fold
            self.run_fold_experiment(fold, model_params, x)


    def get_predictions(self):
        """ Generates predictions for test data for each train-test split, 
        determines the binary threshold on train data, binarises predictions and
        stores them within object in the continuous (y_pred) and binary format
        (y_pred_binary).

        Raises
        ------
        ValueError
            If trained models were not loaded.
        """
        if not self.models_load_check:
            raise ValueError("Trained models were not loaded.")

        info_miss = {
            "ds": [1],
            "type": ["random"],
            "perc": [100]
        }
        
        self.y_pred = []
        self.y_pred_binary = []
    
        for fold in range(self.n_folds):
            # Get data for fold
            train_data = self.dataset.scaled_data[fold]["train"]
            test_data = self.dataset.scaled_data[fold]["test"]
            x = [np.vstack((train_data[i], test_data[i])) for i in range(len(self.dataset.d))]

            # Mask all labels
            x[0] = np.full(x[0].shape, np.nan)
            
            # Predict all labels
            y_pred = GFAtools(x, self.models[fold]).PredictMissing(info_miss)[0]
            y_pred = {
                "train": y_pred[:train_data[0].shape[0]],
                "test": y_pred[train_data[0].shape[0]:]
            }
            
            # Continuous values to sigmoid
            y_pred = {
                "train": 1 / (1 + np.exp(-y_pred["train"])),
                "test": 1 / (1 + np.exp(-y_pred["test"]))
            }

            # Get binary threshold from train data
            fpr, tpr, thresholds = roc_curve(self.dataset.folds_data[fold]["train"][f"{self.dataset.label_modality}_{self.dataset.label_modality}"], y_pred["train"])
            gmeans = np.sqrt(tpr * (1 - fpr))
            train_best_threshold = thresholds[np.argmax(gmeans)]

            # Get binary labels for train and test data
            y_pred_binary = {
                "train": np.where(y_pred["train"] >= train_best_threshold, 1, 0),
                "test": np.where(y_pred["test"] >= train_best_threshold, 1, 0),
            }
            
            self.y_pred.append(y_pred)
            self.y_pred_binary.append(y_pred_binary)
        
        self.prediction_check = True


    def get_performance(self):
        """ Computes classification performance metrics, and stores results
        within object (results).

        Raises
        ------
        ValueError
            If predictions were not made.
        """
        if not self.prediction_check:
            raise ValueError("Predictions haven't been made.")

        results = {
            "Model": list(range(1, self.n_folds + 1)),
            "ELBO": self.elbos,
            "bAcc": [],
            "ROC-AUC": [],
            "Precision": [],
            "Recall": [],
            "F1": [],
            "PPV": [],
            "NPV": [],
            "Sensitivity": [],
            "Specificity": []
        }

        for fold in range(self.n_folds):
            y_true = self.dataset.folds_data[fold]["test"][f"{self.dataset.label_modality}_{self.dataset.label_modality}"]
            y_pred = self.y_pred[fold]["test"]
            y_pred_binary = self.y_pred_binary[fold]["test"]

            tn, fp, fn, tp = confusion_matrix(y_true, y_pred_binary).ravel()

            results["bAcc"].append(balanced_accuracy_score(y_true, y_pred_binary))
            results["ROC-AUC"].append(roc_auc_score(y_true, y_pred))
            results["Precision"].append(precision_score(y_true, y_pred_binary))
            results["Recall"].append(recall_score(y_true, y_pred_binary))
            results["F1"].append(f1_score(y_true, y_pred_binary))
            results["PPV"].append(tp / (tp + fp))
            results["NPV"].append(tn / (tn + fn))
            results["Sensitivity"].append(tp / (tp + fn))
            results["Specificity"].append(tn / (tn + fp))
        
        self.results = pd.DataFrame(data=results)


    def compute_model_variance(self, W: np.ndarray, tau: list[np.ndarray]) -> tuple:
        """ Based on the models weigths/factors and the noise parameters, 
        computes the total variance in data and variance captured by all factors.

        Parameters
        ----------
        W : np.ndarray
            Model factors.
        tau : list[np.ndarray]
            Model noise.

        Returns
        -------
        tuple
            Total variance in data and variance captured by factors.
        """
        # Split W into modalities
        split_indeces = [sum(self.dataset.d[:j+1]) for j in range(len(self.dataset.d) - 1)]
        split_W = np.split(W, split_indeces, axis=0)
        # For each modality
        total_var = 0
        factors_var = 0
        for m in range(len(self.dataset.d)):
            w = split_W[m]
            T = np.diag(1 / tau[m][0, :])
            # Compute total variance
            total_var += np.trace(np.dot(w, w.T) + T)
            # Compute variance captured by factors
            factors_var += np.trace(np.dot(w, w.T))

        return total_var, factors_var


    def compute_factors_variance(self, K: int, W: np.ndarray, total_var) -> tuple:
        """ Based on the models weights/factors, computes the variance within 
        each modality and factors, and the realtive variance captured by each
        factor within a modality.

        Parameters
        ----------
        K : int
            Number of factors.
        W : np.ndarray
            Model factors.
        total_var : _type_
            Total variance.

        Returns
        -------
        tuple
            Variances and relative variances within each modality and factor.
        """
        # Split W into modalities
        split_indeces = [sum(self.dataset.d[:j+1]) for j in range(len(self.dataset.d) - 1)]
        split_W = np.split(W, split_indeces, axis=0)
        
        # Compute the variance within each modality an factors
        var_within = np.zeros((len(self.dataset.d), K))
        non_rel_var_within = np.zeros((len(self.dataset.d), K))
        for m in range(len(self.dataset.d)):
            for k in range(K):
                non_rel_var_within[m][k] = np.sum(split_W[m][:, k] ** 2)
                var_within[m][k] = non_rel_var_within[m][k] / total_var * 100
        
        # Compute realtive variance captured by each factor within a modality
        rel_var_within = np.zeros((len(self.dataset.d), K))
        for m in range(len(self.dataset.d)):
            for k in range(K):
                rel_var_within[m][k] = var_within[m][k] / np.sum(var_within[m][:]) * 100
        
        return var_within, rel_var_within


    def compute_variances(self):
        """ Computes model and factors variances for trained models, and stores 
        them within object (total_vars, factors_vars, vars_within, 
        rel_vars_within).

        Raises
        ------
        ValueError
            If factors or taus were not loaded.
        """
        if not self.factors_load_check or not self.taus_load_check:
            raise ValueError("Factors and/or taus were not loaded.")
        
        self.total_vars = []
        self.factors_vars = []
        for i in range(self.n_folds):
            total_var, factors_var = self.compute_model_variance(self.Ws[i], self.taus[i])
            self.total_vars.append(total_var)
            self.factors_vars.append(factors_var)
        
        self.vars_within = []
        self.rel_vars_within = []
        for i in range(self.n_folds):
            var_within, rel_var_within = \
                self.compute_factors_variance(self.Ks[i], self.Ws[i], self.total_vars[i])
            self.vars_within.append(var_within)
            self.rel_vars_within.append(rel_var_within)
        
        self.factors_variances_check = True
    

    def sort_factors(self):
        """ Sorts factors based on variances in descending order, and stores
        then within object (factors) for clustering.

        Raises
        ------
        ValueError
            If factors variances were not computed.
        """
        if not self.factors_variances_check:
            raise ValueError("Factors variances were not computed.")
        
        self.factors = []
        for i in range(self.n_folds):
            W = self.Ws[i]
            var_within_factors = np.sum(self.vars_within[i], axis=0)
            sorted_order = np.argsort(-var_within_factors)
            W = W[:, sorted_order]
            self.Ws[i] = W
            self.vars_within[i] = self.vars_within[i][:, sorted_order]
            self.rel_vars_within[i] = self.rel_vars_within[i][:, sorted_order]
            
            split_factors = np.split(W, W.shape[1], axis=1)
            self.factors = self.factors + split_factors


    def cluster_factors(self, cosine_similarity_threshold: float):
        """ Clusters all factors from all models to identify factors that were
        inferred by multiple models based on absolute cosine distance. Stores
        clustering labels within object (factors_cluster_labels) used to 
        determine stable factors.

        Parameters
        ----------
        cosine_similarity_threshold : float
            Threshold for cosine similarity above which factors are cosidered 
            the same
        """
        factors_matrix = np.squeeze(np.array(self.factors))
        clustering = DBSCAN(eps=1 - cosine_similarity_threshold,
                            min_samples=1,
                            metric=lambda x, y: 1 - abs(cosine_similarity(x.reshape(1, -1), y.reshape(1, -1))),
                            n_jobs=5)
        clustering = clustering.fit(factors_matrix)
        self.factors_cluster_labels = clustering.labels_

        self.factors_clustered_check = True
    

    def get_stable_factors(self, stability_threshold: float):
        """ Based on the clustering results and the stability thresholds, gets
        the stable factors and stores their information within the object
        (stable_factors), including the weights of the factors from the 
        best model, the weights of the factors across all models that identified
        it, and the variances of the factor.

        Parameters
        ----------
        stability_threshold : float
            Threshold for percentage of models a factor has to be identified by
            to be considered stable.
        
        Raises
        ------
        ValueError
            If factors were not clustered.
        """
        if not self.factors_clustered_check:
            raise ValueError("Factors were not clustered.")

        # Auxiliary lists to identify the models the factors were generated by
        factor_model_label = []
        for i, k in enumerate(self.Ks):
            factor_model_label = factor_model_label + [i] * k
        factor_model_label = np.array(factor_model_label)
        factor_labels_per_model = np.split(self.factors_cluster_labels, [sum(self.Ks[:i]) for i in range(1, len(self.Ks))])
        unique_factors_labels_per_model = [set(factors) for factors in factor_labels_per_model]

        # Get stable factors labels
        factors_prevalence = Counter(self.factors_cluster_labels)
        min_count = stability_threshold * self.n_folds
        stable_factors_labels = [label for label, count in factors_prevalence.items() if count >= min_count]

        # Save factors weights and variances for visualisations
        self.stable_factors = []
        for label in stable_factors_labels:
            # Get stable factor information for detailed plots
            factor_info = self.get_factor_info(factor_labels_per_model, unique_factors_labels_per_model, label)
            self.stable_factors.append(factor_info)


    def get_factor_info(self, factor_labels_per_model: list, unique_factors_labels_per_model: list, label: int) -> dict:
        """ Compiles information about a stable latent factor and returns a 
        structure with it, including the factor weights from the best model, 
        all factor weights from all models that identified it, and the factors 
        variances. All information is needed for visualisations.

        Parameters
        ----------
        factor_labels_per_model : list
            List of clustering labels for each factor identified by each model.
        unique_factors_labels_per_model : list
            List of unique clustering factor labels identified by each model.
        label : int
            Label of the factor for which the information is compiled.

        Returns
        -------
        dict
            Structure with factor information, including its weights from the 
            best model, all its weights from all models that identified it, and 
            its variances.
        """
        factor_info = {}
        
        # Get factor weight
        models = [i for i in range(self.n_folds) if label in unique_factors_labels_per_model[i]]
        elbos = [self.elbos[i] for i in models]
        best_model_idx = models[np.argmax(elbos)]
        factor_idx = np.where(factor_labels_per_model[best_model_idx] == label)[0][0]
        factor_info["W"] = self.Ws[best_model_idx][:, factor_idx]
            
        # Get factor variances
        factor_info["Vars within"] = self.vars_within[best_model_idx][:, factor_idx]
        total_var = sum(factor_info["Vars within"])
        factor_info["Vars percentage"] = np.floor((factor_info["Vars within"] / total_var * 100) * 100) / 100
        factor_info["Vars relative"] = self.rel_vars_within[best_model_idx][:, factor_idx]
        
        # Get weights from all models
        factor_info["all_Ws"] = np.array([self.Ws[model][:, np.where(factor_labels_per_model[model] == label)[0][0]] for model in models])

        return factor_info

        