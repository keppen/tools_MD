import math
import os
import re
import sys
import time

import numpy as np
import pandas as pd

import graph_mol
from libmath import (
    calcAngle,
    calcCentroid,
    calcDistance_form_Vector,
    calcLinearRegression_PowerIteration,
    calcMiddlePoint,
    calcPlane,
    calculate_dihedral,
    calculate_kde,
    calcVector,
)
from libread import read_pdb
from libwrite import write_log


class Visualize:
    def __init__(
        self,
        structure_pdb,
        cluster_pdb,
        log_file=None,
        options="all",
        truncate=None,
        step=None,
        no_plot=True,
        to_csv=True,
    ):
        """
        Initialize the Visualize class.

        Parameters:
            -s, --structure | structure_pdb (str): Path to the PDB file
                    containing the structure data.
            -c, --cluster   | cluster_pdb (str): Path to the PDB file
                    containing cluster data.
            -l, --log-file  | log_file (str): Path to the log file for
                    storing analysis information.
            -o, --options   | options (str, optional): Options for data
                    analysis. Default is "all".
                    d - angle statistics:
                        1 - 1d dihedral angles distribtion analysis
                        2 - 2d dihedral angles distribtion analysis
                        3 - 3d dihedral angles distribtion analysis
                        Ex. "d1" does 1d Ramachandran 1d,
                        "d23" to do 2d and 3d Ramachandran.
                    g - spatial group arrangament analysis
                    h - hydrogen bonds matrix
                    a - helix axis derivation
                    c - contact map
                    It works just like tar options. "d23ch" will run 2d and 3d
                    dihedral distribtions, hbond matrix and contact map. Order
                    does not matter.
            -t, --truncate  | truncate (int, optional): Model index to truncate
                    data collection.
            --step          | step (int, optional): Step size for collecting
                    data.
            --no-plot       | no_plot (bool, optional): If True, suppress
                    plotting. Default is True.
            --no-csv        | to_csv (bool, optional): If True, store data in
                    CSV files. Default is True.
        """
        self.structure_pdb = structure_pdb
        self.cluster_pdb = cluster_pdb
        self.log_file = log_file
        self.MOL = graph_mol.Molecule()
        self.truncate = int(truncate) if truncate else truncate
        self.options = "dghac" if options == "all" else options
        self.step = int(step) if step else 1

        self.no_plot = no_plot
        self.to_csv = to_csv
        self.no_limit = True

        self.DF_pdb: pd.DataFrame
        self.DF_dihedrals: pd.DataFrame
        self.DF_geometry: pd.DataFrame
        self.DF_hbond: pd.DataFrame
        self.DF_axis: pd.DataFrame
        self.DF_distance: pd.DataFrame

    def _init_functions(self):
        self.runs_function = [
            self._get_dih,
            self._get_geometry,
            self._get_hbond,
            self._get_axis,
            self._get_distance,
        ]
        self.calculation_function = [
            self._do_ramachandran,
            self._do_geometry,
            self._do_hbond,
            self._do_axis,
            self._do_contact,
        ]
        self.dataframe_funtion = [
            self._dihedral2DF,
            self._geometry2DF,
            self._hbond2DF,
            self._axis2DF,
            self._distance2DF,
        ]
        self.rama_function = [
                self._ramachandran_1d,
                self._ramachandran_2d,
                self._ramachandran_3d,
                ]
        self.runs_bool = [
            True if c in self.options else False for c in "dghac"
            ]
        self.rama_bool = [
                True if c in self.options else False for c in "123"
                ]

    def _init_dicts(self):
        self.models_dihedral = []
        self.resids_dihedral = []
        self.psi = []
        self.phi = []
        self.omega1 = []
        self.omega2 = []

        self.dict_dihedral = {
            "model": self.models_dihedral,
            "residue index": self.resids_dihedral,
            "Phi": self.phi,
            "Psi": self.psi,
            "Xi": self.omega1,
            "Chi": self.omega2,
        }

        self.models_geometry = []
        self.resids_geometry = []
        self.alpha = []
        self.theta1 = []
        self.theta2 = []

        self.dict_geometry = {
            "model": self.models_geometry,
            "residue index": self.resids_geometry,
            "Alpha": self.alpha,
            "Theta1": self.theta1,
            "Theta2": self.theta2,
        }

        self.models_hbond = []
        self.resid_acceptor = []
        self.resid_donor = []
        self.angle = []
        self.lenght = []

        self.dict_hbond = {
            "model": self.models_hbond,
            "residue acceptor": self.resid_acceptor,
            "residue donor": self.resid_donor,
            "angle": self.angle,
            "lenght": self.lenght,
        }

        self.models_axis = []
        self.resids_axis = []
        self.index = []
        self.distance_axis = []
        self.dict_axis = {
            "model": self.models_axis,
            "index": self.index,
            "distance": self.distance_axis,
            "resid": self.resids_axis,
        }

        self.models_distance = []
        self.index_distance_1 = []
        self.index_distance_2 = []
        self.distance_distance = []

        self.dict_distance = {
            "models": self.models_distance,
            "index 1": self.index_distance_1,
            "index 2": self.index_distance_2,
            "distance": self.distance_distance,
        }

    def _init_log(self):
        if not self.log_file:
            self.log_file = os.getcwd() + "/log_file.log"
            file = open(self.log_file, 'w')
            print(self.log_file)
            file.close()
            return 0
        else:
            file = open(self.log_file, 'r')
            content = file.readlines()
            file.close()
            return content[0].split(';')[0]

    def _update_log(self):
        run_options = [
            self.run_id,
            self.structure_pdb,
            self.cluster_pdb,
            self.truncate,
            self.step,
            self.plot_type,
            self.plot_name,
            self.plot_labels,
            self.plot_limits,
            self.plot_resolution,
            self.data_npz,
        ]
        print(run_options)
        msg = ';'.join([str(i) for i in run_options])
        write_log(self.log_file, msg)

        self.run_id += 1

    def _init_run(self):
        """Run the calculation and generation of plots"""
        self._prepare_mol()

        data = open(self.cluster_pdb, "r")
        self.pdb_gen = read_pdb(data)

        self._init_dicts()
        self._init_functions()
        self.run_id = int(self._init_log()) + 1
        # print(self.MOL.atoms)

        if self.step or self.truncate:
            self._run_truncated()
        else:
            self._run_all()
        self._map_runs(self.dataframe_funtion)

        self._map_runs(self.calculation_function)

    def _run_all(self):
        while self.no_limit:
            self._get_data()
            self._map_runs(self.runs_function, self.runs_bool)
            print("MODEL: ", self.model, end="\r")
            # self._printProgress(self.model)

    def _run_truncated(self):
        while self.no_limit:
            self._get_data()

            if self.model == self.truncate:
                break
            if self.model % self.step != 0:
                continue

            # self._printProgress(self.model)
            print("MODEL: ", self.model, end="\r")
            self._map_runs(self.runs_function, self.runs_bool)

    def _do_ramachandran(self):
        self._map_runs(self.rama_function, self.rama_bool)

    def _map_runs(self, function_list, bool_list):
        for condition, function in zip(bool_list, function_list):
            if condition:
                function()

    def _split_name(self, file):
        """Split name of file and get the information about analyzed file.

        Parameters:
            file(str): a file to be analyzed. e.g. 10ns-a1.pdb001.pdb

        Returns:
            name(str): a na of a system
            cluster(str): a number of a cluster
        """
        file = re.split("\-|\_|\.", file)
        cluster = file[-2]
        nums = [str(_) for _ in range(1, 10)]
        for letter in cluster:
            if letter in nums:
                break
            cluster = cluster[1:]

        name = file[1]
        print(f"CHECK THE FILE NAME AND PICK CAREFULLY: {file}")
        print(f"indexes picked: \n\tcluser: {cluster}\n\tname: {name}")
        return name, cluster

    def _read_file(self, file):
        """Helper function"""
        with open(file, "r") as f:
            return f.readlines()

    def _prep_line(self, line: str):
        """Helper function"""
        return line.strip().split()

    def _prepare_mol(self):
        """Helper function"""
        structure = self._read_file(self.structure_pdb)
        for line in structure:
            line = self._prep_line(line)
            if "ATOM" in line or "HETATM" in line:
                self.MOL.add_atom(line[1], line[11])
            if line[0] == "CONECT":
                self.MOL.add_bond(line[1:])
        self.MOL.run()
        self.cluster_pdb = self.MOL.newPDB(self.cluster_pdb)

    def _printProgress(self, model):
        if model < 1000 and model % 200 == 0:
            print(f"Collecting model: {model}", end="\r")
        if model >= 1000 and model % 1000 == 0:
            print(f"Collecting model: {model}", end="\r")

    def _center_points(self):
        points = []
        # for index, row in self.DF_pdb.iterrows():
        #     points.append(np.array([row.x, row.y, row.z]))
        for index in self.MOL.backbone:
            p = self._search_by_index(self.DF_pdb, int(index))
            points.append(np.array(p))
        geom_center = calcCentroid(points)

        self.points_backbone = []
        for index, row in self.DF_pdb.iterrows():
            p = np.array([row.x, row.y, row.z]) - geom_center
            self.DF_pdb.iloc[index, 7] = p[0]
            self.DF_pdb.iloc[index, 8] = p[1]
            self.DF_pdb.iloc[index, 9] = p[2]
            self.points_backbone.append(
                p if str(index) in self.MOL.backbone else np.zeros(1)
            )
        self.points_backbone = [p for p in points if p.all()]

    def _get_data(self):
        try:
            self.model, pdb_data = next(self.pdb_gen)
            # print(pdb_data)
            self.DF_pdb = pd.DataFrame.from_dict(
                pdb_data, orient="index").transpose()
            self._center_points()
        except StopIteration:
            self.no_limit = False

    def _get_dih(self):
        """Calculate dihedrals with data collected from pdb file.
        Dihedrals are stored in DataFrame object

        It iterates over all frames and (min + 1; max -1) residues indexes range
        It follows the same algorithm as _collectData method. Helper sunftions searches
        for coordianates and calculates diherdal angles.
        Truncating DataFrame by redundant data is not speeding up searching process
        """

        max_resid = self.DF_pdb.max()["residue_seqnumber"]

        for resid in range(2, max_resid):
            self.models_dihedral.append(self.model)
            self.resids_dihedral.append(resid)
            self.phi.append(self._calculate_phi(self.DF_pdb, resid))
            self.psi.append(self._calculate_psi(self.DF_pdb, resid))
            self.omega1.append(self._calculate_xi(self.DF_pdb, resid))
            self.omega2.append(self._calculate_chi(self.DF_pdb, resid))

    def _get_geometry(self):
        """Calculate dihedrals with data collected from pdb file.
        Dihedrals are stored in DataFrame object

        It iterates over all frames and (min + 1; max -1) residues indexes range
        It follows the same algorithm as _collectData method. Helper sunftions searches
        for coordianates and calculates diherdal angles.
        Truncating DataFrame by redundant data is not speeding up searching process
        THERE IS SOMETHING WRONG
        """
        max_resid = self.DF_pdb.max()["residue_seqnumber"]

        # print(self.model)

        res0 = None
        res1 = None
        for resid in range(2, max_resid):
            res0 = res1 or self._calculate_plane_and_geomcentre(
                self.DF_pdb, resid)
            res1 = self._calculate_plane_and_geomcentre(self.DF_pdb, resid + 1)
            a = calcAngle(res0[0], res1[0])
            a = np.degrees(math.acos(a))
            line = calcVector(res0[1], res1[1])
            t1 = calcAngle(res0[0], line)
            t2 = calcAngle(res1[0], line)
            t1 = np.degrees(math.asin(abs(t1)))
            t2 = np.degrees(math.asin(abs(t2)))
            self.models_geometry.append(self.model)
            self.resids_geometry.append(resid)
            self.alpha.append(a)
            self.theta1.append(t1)
            self.theta2.append(t2)
            print(
                f"plane1: {res0}\nplane2: {res1}\nalpha: {a}\ntheta1: {t1}\ntheta2: {t2}"
            )

            # libplot.debug_geometry(res0[0], res1[0], res0[1], res1[1])

    def _get_hbond(self):
        max_resid = self.DF_pdb.max()["residue_seqnumber"]

        for acceptor in range(2, max_resid + 1):
            for donor in range(1, max_resid):
                angle, lenght = self._calculate_hbond(
                    self.DF_pdb, acceptor, donor)
                angle = math.acos(angle) / math.pi * 180
                self.models_hbond.append(self.model)
                self.resid_acceptor.append(acceptor)
                self.resid_donor.append(donor)
                self.angle.append(angle)
                self.lenght.append(lenght)

    def _get_axis(self):
        self.axis = calcLinearRegression_PowerIteration(self.points_backbone)
        for index in self.MOL.backbone:
            self.models_axis.append(self.model)
            point = self._search_by_index(self.DF_pdb, int(index))
            self.distance_axis.append(
                calcDistance_form_Vector(point, self.axis))
            self.index.append(self.MOL.backbone.index(index))
            self.resids_axis.append(self.MOL._findbyID(index)["resid"])

    def _get_distance(self):
        atoms_indexes = []
        for index in self.MOL.backbone:
            atom = self.DF_pdb.loc[
                # (df_model['model'] == model) &
                (self.DF_pdb["atom_serialnumber"] == int(index))
            ]
            if atom.iloc[0, 1] in ["C", "CB"]:
                atoms_indexes.append(atom.iloc[0, 0])

        for index_1, atom_id_1 in enumerate(atoms_indexes):
            for index_2, atom_id_2 in enumerate(atoms_indexes):
                atom_1 = self._search_by_index(self.DF_pdb, int(atom_id_1))
                atom_2 = self._search_by_index(self.DF_pdb, int(atom_id_2))
                distance = np.linalg.norm(np.array(atom_1) - np.array(atom_2))
                self.models_distance.append(self.model)
                self.index_distance_1.append(index_1)
                self.index_distance_2.append(index_2)
                self.distance_distance.append(distance)

    def _dihedral2DF(self):
        self.DF_dihedrals = pd.DataFrame.from_dict(
            self.dict_dihedral, orient="index"
        ).transpose()

    def _hbond2DF(self):
        self.DF_hbond = pd.DataFrame.from_dict(
            self.dict_hbond, orient="index"
        ).transpose()
        self.DF_hbond["Hbond presence"] = np.where(
            (self.DF_hbond["angle"] <= 30) & (
                self.DF_hbond["lenght"] <= 3.5), 1, 0
        )
        self.pivot_hbond = pd.pivot_table(
            self.DF_hbond,
            values="Hbond presence",
            index="residue acceptor",
            columns="residue donor",
            aggfunc=np.mean,
        )

    def _axis2DF(self):
        self.DF_axis = pd.DataFrame.from_dict(
            self.dict_axis, orient="index"
        ).transpose()

    def _distance2DF(self):
        self.DF_distance = pd.DataFrame.from_dict(
            self.dict_distance, orient="index"
        ).transpose()
        self.pivot_distance = pd.pivot_table(
            self.DF_distance,
            values="distance",
            index="index 1",
            columns="index 2",
            aggfunc=np.mean,
        )

    def _search_by_type(self, df_model, resid, atom_type):
        """Helper function"""
        p = df_model.loc[
            # (df_model['model'] == model) &
            (df_model["residue_seqnumber"] == resid)
            & (df_model["atom_name"] == atom_type)
        ]
        # for a, b in p.to_numpy():
        #     print(f'{a}: {b}')
        return [float(p.iat[0, 7]), float(p.iat[0, 8]), float(p.iat[0, 9])]

    def _search_by_index(self, df_model, index):
        """Helper function"""
        p = df_model.loc[
            # (df_model['model'] == model) &
            (df_model["atom_serialnumber"] == index)
        ]
        # df1 = p.stack().reset_index(level=0, drop=True).rename_axis('a').reset_index(name='b')
        # for a, b in df1.to_numpy():
        #     print(f'{a}: {b}')
        return [float(p.iat[0, 7]), float(p.iat[0, 8]), float(p.iat[0, 9])]

    def _calculate_phi(self, df, resid):
        """Helper function"""
        p1 = self._search_by_type(df, resid - 1, "C")
        p2 = self._search_by_type(df, resid, "N")
        p3 = self._search_by_type(df, resid, "CG")
        p4 = self._search_by_type(df, resid, "CB")
        return calculate_dihedral(p1, p2, p3, p4)

    def _calculate_psi(self, df, resid):
        """Helper function"""
        p1 = self._search_by_type(df, resid, "CB")
        p2 = self._search_by_type(df, resid, "OA")
        p3 = self._search_by_type(df, resid, "C")
        p4 = self._search_by_type(df, resid + 1, "N")
        return calculate_dihedral(p1, p2, p3, p4)

    def _calculate_xi(self, df, resid):
        """Helper function"""
        p1 = self._search_by_type(df, resid, "N")
        p2 = self._search_by_type(df, resid, "CG")
        p3 = self._search_by_type(df, resid, "CB")
        p4 = self._search_by_type(df, resid, "OA")
        return calculate_dihedral(p1, p2, p3, p4)

    def _calculate_chi(self, df, resid):
        """Helper function"""
        p1 = self._search_by_type(df, resid, "CG")
        p2 = self._search_by_type(df, resid, "CB")
        p3 = self._search_by_type(df, resid, "OA")
        p4 = self._search_by_type(df, resid, "C")
        return calculate_dihedral(p1, p2, p3, p4)

    def _calculate_plane_and_geomcentre(self, df, resid):
        p0 = self._search_by_type(df, resid - 1, "O")
        p1 = self._search_by_type(df, resid - 1, "OA")
        p2 = self._search_by_type(df, resid, "N")
        p3 = self._search_by_type(df, resid - 1, "C")
        # print(p0, p1, p2, p3)
        return calcPlane(p0, p2, p1), p3

    def _calculate_hbond(self, df, acceptor, donor):
        p0 = self._search_by_type(df, acceptor, "N")
        p1 = self._search_by_type(df, acceptor - 1, "C")
        p2 = self._search_by_type(df, acceptor, "CG")
        p3 = self._search_by_type(df, donor, "O")
        p4 = calcMiddlePoint(p1, p2)
        line1 = calcVector(p4, p0)
        line2 = calcVector(p0, p3)
        return calcAngle(line1, line2), np.linalg.norm(line2)

    def _do_hbond(self):
        import matplotlib.pyplot as plt
        import seaborn as sb

        name, cluster = self._split_name(self.cluster_pdb)

        sb.set_context("paper", font_scale=1.35, rc={"lines.linewidth": 0.85})

        print(self.pivot_hbond)
        f, ax = plt.subplots()
        sb.heatmap(
            self.pivot_hbond,
            cmap="crest",
            vmin=0,
            vmax=1,
        )
        f.savefig(
            f"hbond_{name}_{cluster}.png",
            dpi=100,
        )

        # plot.savefig(f"{'test1.png' if i == 0 else 'test2'}", dpi=1000)
        plt.close()

    def _do_contact(self):
        import matplotlib.pyplot as plt
        import seaborn as sb
        name, cluster = self._split_name(self.cluster_pdb)

        sb.set_context("paper", font_scale=1.35, rc={"lines.linewidth": 0.85})

        print(self.pivot_distance)
        f, ax = plt.subplots()
        sb.heatmap(
            self.pivot_distance,
            cmap="crest",
            vmin=0,
            vmax=8,
        )
        f.savefig(
            f"contact_{name}_{cluster}.png",
            dpi=100,
        )

        # plot.savefig(f"{'test1.png' if i == 0 else 'test2'}", dpi=1000)
        plt.close()

    def _ramachandran_3d(self):
        from intertools import chain
        from libplot import mavi_contour

        DF_angles = self.DF_dihedrals[["Phi", "Xi", "Chi"]]
        np_angles = DF_angles.to_numpy(dtype=np.float64)

        grid = np.array([
            [-180, 180],
            [-180, 180],
            [-180, 180]
        ])
        resolution = 180
        name = "test"
        labels = DF_angles.columns.tolist()

        density, coordinates = calculate_kde(np_angles, grid, resolution)

        if self.to_csv:

            self.data_npz = f"{self.run_id}_data.npz"

            np.savez_compressed(
                    self.data_npz,
                    data=density,
                    coords=coordinates
                    )

            self.plot_type = "rama_3d"
            self.plot_name = name
            self.plot_labels = " ".join(labels)
            self.plot_limits = " ".join([str(i) for i in chain(*grid)])
            self.plot_resolution = resolution

            self._update_log()

        if not self.no_plot:
            mavi_contour(density, coordinates, limits=grid,
                         name=name, labels=labels)

    def _ramachandran_2d(self):
        from itertools import chain, combinations

        from libplot import ramachandran_plot

        DF_angles = self.DF_dihedrals[["Phi", "Xi", "Chi"]]
        DF_combinations = combinations(DF_angles, 2)

        grid = np.array(
            [[-180, 180],
             [-180, 180]]
        )

        resolution = 180

        for comb in DF_combinations:
            DF_data = (DF_angles[list(comb)])
            np_angles = DF_data.to_numpy(dtype=np.float64)

            name = f"test-{''.join(comb)}"
            labels = DF_angles.columns.tolist()

            density, coordinates = calculate_kde(np_angles, grid, resolution)

            if self.to_csv:

                self.data_npz = f"{self.run_id}_data.npz"

                np.savez_compressed(
                        self.data_npz,
                        data=density,
                        coords=coordinates
                        )

                self.plot_type = "rama_2d"
                self.plot_name = name
                self.plot_labels = " ".join(labels)
                self.plot_limits = " ".join([str(i) for i in chain(*grid)])
                self.plot_resolution = resolution

                self._update_log()

            if not self.no_plot:
                ramachandran_plot(density, coordinates,
                                  name=name, limits=grid, labels=labels)

    def _ramachandran_1d(self):
        from libplot import distribution_plot

        DF_angles = self.DF_dihedrals[["Phi", "Xi", "Chi"]]
        np_angles = DF_angles.to_numpy(dtype=np.float64)

        grid = np.array(
            [[-180, 180]]
        )
        resolution = 180

        distribution_plot(np_angles, DF_angles.columns.tolist(),
                          "test", grid, resolution)

    def _do_geometry(self):
        import matplotlib.pyplot as plt
        import seaborn as sb

        name, cluster = self._split_name(self.cluster_pdb)

        letters = {
            "greek": ["\u03b1", "\u03b8\u2081", "\u03b8\u2082"],
            "latin": ["alpha", "theta1", "theta2"],
        }

        sb.set_context("paper", font_scale=1.35, rc={"lines.linewidth": 0.85})
        colors = ["black", "red", "blue"]

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7, 5), sharex=True)
        x1 = self.DF_geometry.Alpha
        x2 = self.DF_geometry.Theta1
        x3 = self.DF_geometry.Theta2

        sb.kdeplot(
            data=self.DF_geometry, x=x1, hue="residue index", palette=colors, ax=ax1
        )
        ax1.axhline(0, color="k", clip_on=False)
        ax1.set_ylabel(letters["greek"][0])
        ax1.legend(labels=["2", "3", "4"], title="No. Mer", frameon=False)

        sb.kdeplot(
            data=self.DF_geometry,
            x=x2,
            hue="residue index",
            palette=colors,
            ax=ax2,
            legend=False,
        )
        ax2.axhline(0, color="k", clip_on=False)
        ax2.set_ylabel(letters["greek"][1])

        sb.kdeplot(
            data=self.DF_geometry,
            x=x3,
            hue="residue index",
            palette=colors,
            ax=ax3,
            legend=False,
        )
        ax3.axhline(0, color="k", clip_on=False)
        ax3.set_ylabel(letters["greek"][2])

        ax3.set_xlim(0, 180)
        ax3.set_xticks([_ for _ in range(0, 181, 30)])
        ax3.set_xlabel("Degrees")

        fig.savefig(f"geom_{name}_{cluster}.png", dpi=500)
        plt.close()

    def _do_axis(self):
        import matplotlib.pyplot as plt
        import seaborn as sb

        name, cluster = self._split_name(self.cluster_pdb)

        plot = sb.lineplot(data=self.DF_axis, x="index",
                           y="distance", marker="o")
        plot.figure.savefig(f"axis_{name}_{cluster}.png", dpi=500)
        plt.close()


if __name__ == "__main__":
    print(sys.argv)

    truncate = None
    options = None
    step = None
    no_plot = None
    to_csv = None
    log_file = None
    structure = None
    cluster = None

    start_time = time.time()
    for index, arg in enumerate(sys.argv):
        if arg in ["-t", "--truncate"]:
            truncate = sys.argv[index + 1]
        if arg in ["-o", "--options"]:
            options = sys.argv[index + 1]
        if arg in ["--step"]:
            step = sys.argv[index + 1]
        if arg in ["--no-plot"]:
            no_plot = True
        if arg in ["-l", "--log-file"]:
            log_file = sys.argv[index + 1]
        if arg in ["-s", "--structure"]:
            structure = sys.argv[index + 1]
        if arg in ["-c", "--cluster"]:
            cluster = sys.argv[index + 1]
        if arg in ["--no-csv"]:
            to_csv = False

    R = Visualize(
        structure,
        cluster,
        truncate=truncate,
        step=step,
        options=options,
        no_plot=no_plot,
        to_csv=to_csv,
        log_file=log_file
    )
    R._init_run()

    # R._collectDih()
    # R._rama_plot()
    #
    # R._collectGeometry()
    # R._geom_plot()
    #
    # R._collectHbond()
    # R._hbond_plot()
    #
    # R._collectAxis()
    # R._axis_plot()
    print("--- %s seconds ---" % (time.time() - start_time))
