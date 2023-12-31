import math
import os
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
from utils import split_name
from liblogg import timing


class Visualize:
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

    def __init__(
        self,
        structure_pdb,
        cluster_pdb,
        log_file=None,
        options="all",
        truncate=None,
        step=None,
        no_plot=False,
        to_csv=True,
    ):
        """
        Initialize the Visualize class.
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
        self.DF_contact: pd.DataFrame

    def _init_functions(self):
        self.runs_function = [
            self._get_dih,
            self._get_geometry,
            self._get_hbond,
            self._get_axis,
            self._get_contact,
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
            self._contact2DF,
        ]
        self.rama_function = [
            self._do_rama_1d,
            self._do_rama_2d,
            self._do_rama_3d,
        ]
        self.geom_function = [
            self._do_geom_1d,
            self._do_geom_2d,
            self._do_geom_3d,
        ]

        self.runs_bool = [
            True if c in self.options else False for c in "dghac"
        ]
        self.dimmentions_bool = [
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

        self.dict_contact = {
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
            no_entry = 0
        else:
            file = open(self.log_file, 'r')
            content = file.readlines()
            file.close()
            no_entry = content[-1].split(';')[0]

        self.data_id = int(no_entry) + 1

    def _update_log(self):
        run_options = [
            self.data_id,
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
        # print(run_options)
        msg = ';'.join([str(i) for i in run_options])
        write_log(self.log_file, msg)

        self.data_id += 1

    def _init_run(self):
        """Run the calculation and generation of plots"""
        self._prepare_mol()

        data = open(self.cluster_pdb, "r")
        self.pdb_gen = read_pdb(data)

        self._init_dicts()
        self._init_functions()
        self._init_log()
        # print(self.MOL.atoms)

        if self.step or self.truncate:
            self._run_truncated()
        else:
            self._run_all()
        self._map_runs(self.dataframe_funtion, self.runs_bool)

        self._map_runs(self.calculation_function, self.runs_bool)

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
        self._map_runs(self.rama_function, self.dimmentions_bool)

    def _do_geometry(self):
        self._map_runs(self.geom_function, self.dimmentions_bool)

    def _map_runs(self, function_list, bool_list):
        for condition, function in zip(bool_list, function_list):
            if condition:
                function()

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
            # print(
            #     f"plane1: {res0}\nplane2: {res1}\nalpha: {a}\ntheta1: {t1}\ntheta2: {t2}"
            # )

            # libplot.debug_geometry(res0[0], res1[0], res0[1], res1[1])

    def _get_hbond(self):
        max_resid = self.DF_pdb.max()["residue_seqnumber"]

        for donor in range(2, max_resid + 1):
            for acceptor in range(1, max_resid):
                angle, lenght = self._calculate_hbond(
                    self.DF_pdb, donor, acceptor)
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

            self.index.append(index)

            # self.index.append(self.MOL.backbone.index(index))
            self.resids_axis.append(self.MOL._findbyID(index)["resid"])

    def _get_contact(self):
        atoms_indexes = []
        for index in self.MOL.backbone:
            atom = self.DF_pdb.loc[
                # (df_model['model'] == model) &
                (self.DF_pdb["atom_serialnumber"] == int(index))
            ]
            if atom.iloc[0, 1] in ["C", "CB"]:
                atoms_indexes.append(atom.iloc[0, 0])

        for atom_id_1 in atoms_indexes:
            for atom_id_2 in atoms_indexes:
                atom_1 = self._search_by_index(self.DF_pdb, int(atom_id_1))
                atom_2 = self._search_by_index(self.DF_pdb, int(atom_id_2))
                distance = np.linalg.norm(np.array(atom_1) - np.array(atom_2))
                self.models_distance.append(self.model)
                self.index_distance_1.append(atom_id_1)
                self.index_distance_2.append(atom_id_2)
                self.distance_distance.append(distance)

    def _dihedral2DF(self):
        self.DF_dihedrals = pd.DataFrame.from_dict(
            self.dict_dihedral, orient="index"
        ).transpose()

        if self.to_csv:
            np_data = self.DF_dihedrals.to_numpy(dtype=np.float64)
            self.data_npz = f"{self.data_id}_rawdata.npz"

            np.savez_compressed(
                self.data_npz,
                rawdata=np_data,
            )

            self.plot_type = "rawdata"
            self.plot_name = "dihedral"
            self.plot_labels = " "
            self.plot_limits = " "
            self.plot_resolution = " "

            self._update_log()

    def _geometry2DF(self):
        self.DF_geometry = pd.DataFrame.from_dict(
            self.dict_geometry, orient="index"
        ).transpose()

        if self.to_csv:
            np_data = self.DF_geometry.to_numpy(dtype=np.float64)
            self.data_npz = f"{self.data_id}_rawdata.npz"

            np.savez_compressed(
                self.data_npz,
                rawdata=np_data,
            )

            self.plot_type = "rawdata"
            self.plot_name = "geometry"
            self.plot_labels = " "
            self.plot_limits = " "
            self.plot_resolution = " "

            self._update_log()

    def _hbond2DF(self):
        self.DF_hbond = pd.DataFrame.from_dict(
            self.dict_hbond, orient="index"
        ).transpose()

        if self.to_csv:
            np_data = self.DF_hbond.to_numpy(dtype=np.float64)
            self.data_npz = f"{self.data_id}_rawdata.npz"

            np.savez_compressed(
                self.data_npz,
                rawdata=np_data,
            )

            self.plot_type = "rawdata"
            self.plot_name = "hbonds"
            self.plot_labels = " "
            self.plot_limits = " "
            self.plot_resolution = " "

            self._update_log()

    def _axis2DF(self):
        self.DF_axis = pd.DataFrame.from_dict(
            self.dict_axis, orient="index"
        ).transpose()

        if self.to_csv:
            np_data = self.DF_axis.to_numpy(dtype=np.float64)
            self.data_npz = f"{self.data_id}_rawdata.npz"

            np.savez_compressed(
                self.data_npz,
                rawdata=np_data,
            )

            self.plot_type = "rawdata"
            self.plot_name = "axis"
            self.plot_labels = " "
            self.plot_limits = " "
            self.plot_resolution = " "

            self._update_log()

    def _contact2DF(self):
        self.DF_contact = pd.DataFrame.from_dict(
            self.dict_contact, orient="index"
        ).transpose()

        if self.to_csv:
            np_data = self.DF_contact.to_numpy(dtype=np.float64)
            self.data_npz = f"{self.data_id}_rawdata.npz"

            np.savez_compressed(
                self.data_npz,
                rawdata=np_data,
            )

            self.plot_type = "rawdata"
            self.plot_name = "contact"
            self.plot_labels = " "
            self.plot_limits = " "
            self.plot_resolution = " "

            self._update_log()

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

    @timing
    def _do_axis(self, data=None, name=None, limits=None, resolution=None):
        from libplot import plt_lineplot

        if data is None:
            data = self.DF_axis

        name = name or "axis-" + "-".join(split_name(self.cluster_pdb))

        data = data[["index", "distance"]]
        data.set_index("index", inplace=True)
        data = data.groupby(data.index).mean()

        labels = ['Backbone atoms', 'Distance from axis [A]']
        limits = limits or [0, data.index.max() + 1]

        if not self.no_plot:
            plt_lineplot(data, labels=labels,
                         limits=limits, name=name)

    @timing
    def _do_hbond(self, data=None, name=None, limits=None, resolution=None):
        from libplot import plt_heatmap
        import numpy as np

        if data is None:
            data = self.DF_hbond

        data["Hbond presence"] = np.where(
            (data["angle"] <= 30) & (
                data["lenght"] <= 3.5), 1, 0
        )

        data = pd.pivot_table(
            data,
            values="Hbond presence",
            index="residue acceptor",
            columns="residue donor",
            aggfunc=np.mean
        )
        print(data)

        name = name or "hbond-" + "-".join(split_name(self.cluster_pdb))

        limits = None
        labels = [
                [f"Mer {int(x)}" for x in data.columns],
                [f"Mer {int(y)}" for y in data.index]
                ]
        ax_labels = ["Donor", "Acceptor"]

        if not self.no_plot:
            plt_heatmap(
                    data,
                    labels=labels,
                    limits=limits,
                    name=name,
                    ax_labels=ax_labels
                    )

    @timing
    def _do_contact(self, data=None, name=None, limits=None, resolution=None):
        from libplot import plt_heatmap

        if data is None:
            data = self.DF_contact

        data = pd.pivot_table(
            data,
            values="distance",
            index="index 1",
            columns="index 2",
            aggfunc=np.mean,
        )

        name = name or "contact-" + "-".join(split_name(self.cluster_pdb))

        limits = None
        labels = [
                [f"Atom {int(x)}" for x in data.columns],
                [f"Atom {int(y)}" for y in data.index]
                ]

        palette = "Greens"

        # ax_labels = ["Donor", "Acceptor"]

        if not self.no_plot:
            plt_heatmap(
                    data,
                    labels=labels,
                    limits=limits,
                    name=name,
                    palette=palette
                        )

    @timing
    def _do_rama_3d(self, data=None, name=None, limits=None, resolution=180):
        from itertools import chain
        from libplot import mavi_contour

        if data is None:
            data = self.DF_dihedrals

        DF_angles = data[["Phi", "Xi", "Chi"]]
        np_angles = DF_angles.to_numpy(dtype=np.float64)

        if not limits:
            limits = np.array([
                [-180, 180],
                [-180, 180],
                [-180, 180]
            ])
        name = name or "rama3d-" + "-".join(split_name(self.cluster_pdb))
        labels = DF_angles.columns.tolist()

        density, coordinates = calculate_kde(np_angles, limits, resolution)

        if self.to_csv:

            self.data_npz = f"{self.data_id}_data.npz"

            np.savez_compressed(
                self.data_npz,
                data=density,
                coords=coordinates
            )

            self.plot_type = "rama_3d"
            self.plot_name = name
            self.plot_labels = " ".join(labels)
            self.plot_limits = " ".join([str(i) for i in chain(*limits)])
            self.plot_resolution = resolution

            self._update_log()

        if not self.no_plot:
            mavi_contour(density, coordinates, limits=limits,
                         name=name, labels=labels)

    @timing
    def _do_rama_2d(self, data=None, name=None, limits=None, resolution=180):
        from itertools import chain, combinations
        from libplot import plt_ramachandran

        if data is None:
            data = self.DF_dihedrals

        DF_angles = data[["Phi", "Xi", "Chi"]]
        DF_combinations = combinations(DF_angles, 2)

        if not limits:
            limits = np.array(
                [[-180, 180],
                 [-180, 180]]
            )

        for comb in DF_combinations:
            DF_data = (DF_angles[list(comb)])
            np_angles = DF_data.to_numpy(dtype=np.float64)

            if name:
                name_sub = f"{name}-{''.join(comb)}"
            else:
                name_sub = f"rama2d-{'-'.join(split_name(self.cluster_pdb))}-{''.join(comb)}"

            labels = DF_data.columns.tolist()

            density, coordinates = calculate_kde(np_angles, limits, resolution)

            if self.to_csv:

                self.data_npz = f"{self.data_id}_data.npz"

                np.savez_compressed(
                    self.data_npz,
                    data=density,
                    coords=coordinates
                )

                self.plot_type = "rama_2d"
                self.plot_name = name_sub
                self.plot_labels = " ".join(list(comb))
                self.plot_limits = " ".join([str(i) for i in chain(*limits)])
                self.plot_resolution = resolution

                self._update_log()

            if not self.no_plot:
                plt_ramachandran(density, coordinates,
                                 name=name_sub, limits=limits, labels=labels)

    def _split_DF(self, data, column_name):
        """
    Split DataFrame into DF with respcect to angle_names and res_ids.
    Return a dictionary where each entry contains a produced Dataframes
    converted into numpy arrays.
        """

        # lists of items to split by
        res_ids = data[column_name].unique()
        angle_names = data.columns.tolist()[1:]

        data_dict = {}
        # filter by res id
        for res_id in res_ids:
            condition = (data[column_name] == res_id)
            by_res_id = data[condition].copy()
            by_res_id.drop(by_res_id.columns[0], axis=1, inplace=True)
            # filter by angle name
            for angle_name in angle_names:
                by_angle = by_res_id[angle_name]
                vector = by_angle.to_numpy(dtype=np.float64).reshape(-1, 1)
                key = (angle_name, int(res_id))
                data_dict[key] = vector

        return data_dict

    def _multiple_kde(self, data, limits, resolution):
        """
    Calculate a kde from data produced by _split_DF.
    Return a two dictionaries containing densites and coordinates.
        """

        density_dict = {}
        coordinate_dict = {}
        for key, item in data.items():
            angle, res_id = key
            density, coordinates = calculate_kde(item, limits, resolution)
            density_dict[f"{angle} {res_id} dens"] = np.array(density)
            coordinate_dict[f"{angle} {res_id} coords"] = np.array(coordinates)

        return density_dict, coordinate_dict

    @timing
    def _do_rama_1d(self, data=None, name=None, limits=None, resolution=180):
        from libplot import plt_distribution
        from itertools import chain

        if data is None:
            data = self.DF_dihedrals

        if not limits:
            limits = np.array([
                [-180, 180],
                [-180, 180],
                [-180, 180]
            ])

        DF_angles = data[["residue index", "Phi", "Xi", "Chi"]]
        data = self._split_DF(DF_angles, "residue index")
        density_dict, coordinates_dict = self._multiple_kde(data, limits[0].reshape(1,-1), resolution)


        name = name or "rama1d-" + "-".join(split_name(self.cluster_pdb))

        labels = DF_angles.columns.tolist()[1:]

        if self.to_csv:

            self.data_npz = f"{self.data_id}_data.npz"

            np.savez_compressed(
                self.data_npz,
                **density_dict,
                **coordinates_dict
            )

            self.plot_type = "rama_1d"
            self.plot_name = name
            self.plot_labels = " ".join(list(labels))
            self.plot_limits = " ".join([str(i) for i in chain(*limits)])
            self.plot_resolution = resolution

            self._update_log()

        if not self.no_plot:
            plt_distribution(density_dict, coordinates_dict,
                             labels=labels,
                             name=name,
                             limits=limits,
                             )

    @timing
    def _do_geom_3d(self, data=None, name=None, limits=None, resolution=90):
        from itertools import chain
        from libplot import mavi_contour

        if data is None:
            data = self.DF_geometry

        DF_angles = data[["Alpha", "Theta1", "Theta2"]]
        np_angles = DF_angles.to_numpy(dtype=np.float64)

        if not limits:
            limits = np.array([
                [0, 180],
                [0, 90],
                [0, 90]
            ])
        name = name or "geom3d-" + "-".join(split_name(self.cluster_pdb))
        labels = DF_angles.columns.tolist()

        density, coordinates = calculate_kde(np_angles, limits, resolution)

        if self.to_csv:

            self.data_npz = f"{self.data_id}_data.npz"

            np.savez_compressed(
                self.data_npz,
                data=density,
                coords=coordinates
            )

            self.plot_type = "geom_3d"
            self.plot_name = name
            self.plot_labels = " ".join(labels)
            self.plot_limits = " ".join([str(i) for i in chain(*limits)])
            self.plot_resolution = resolution

            self._update_log()

        if not self.no_plot:
            mavi_contour(density, coordinates, limits=limits,
                         name=name, labels=labels)

    @timing
    def _do_geom_2d(self, data=None, name=None, limits=None, resolution=90):
        from itertools import chain, combinations
        from libplot import plt_ramachandran

        if data is None:
            data = self.DF_geometry

        DF_angles = data[["Alpha", "Theta1", "Theta2"]]
        DF_combinations = combinations(DF_angles, 2)

        if not limits:
            limits = np.array(
                [[0, 180],
                 [0, 180]]
            )

        for comb in DF_combinations:
            DF_data = (DF_angles[list(comb)])
            np_angles = DF_data.to_numpy(dtype=np.float64)

            if name:
                name_sub = f"{name}-{''.join(comb)}"
            else:
                name_sub = f"geom-2-{'-'.join(split_name(self.cluster_pdb))}-{''.join(comb)}"

            labels = DF_data.columns.tolist()

            density, coordinates = calculate_kde(np_angles, limits, resolution)

            if self.to_csv:

                self.data_npz = f"{self.data_id}_data.npz"

                np.savez_compressed(
                    self.data_npz,
                    data=density,
                    coords=coordinates
                )

                self.plot_type = "geom_2d"
                self.plot_name = name_sub
                self.plot_labels = " ".join(list(comb))
                self.plot_limits = " ".join([str(i) for i in chain(*limits)])
                self.plot_resolution = resolution

                self._update_log()

            if not self.no_plot:
                plt_ramachandran(density, coordinates,
                                 name=name_sub, limits=limits, labels=labels)

    @timing
    def _do_geom_1d(self, data=None, name=None, limits=None, resolution=90):
        from libplot import plt_distribution
        from itertools import chain

        if data is None:
            data = self.DF_geometry

        if not limits:
            limits = np.array(
                    [[0, 180],
                     [0, 90],
                     [0, 90]]
            )

        DF_angles = data[["residue index", "Alpha", "Theta1", "Theta2"]]
        data = self._split_DF(DF_angles, "residue index")
        density_dict, coordinates_dict = self._multiple_kde(data, limits[0].reshape(1,-1), resolution)

        name = name or "geom1d-" + "-".join(split_name(self.cluster_pdb))

        labels = DF_angles.columns.tolist()[1:]

        if self.to_csv:

            self.data_npz = f"{self.data_id}_data.npz"

            np.savez_compressed(
                self.data_npz,
                **density_dict,
                **coordinates_dict
            )

            self.plot_type = "geom_1d"
            self.plot_name = name
            self.plot_labels = " ".join(list(labels))
            self.plot_limits = " ".join([str(i) for i in chain(*limits)])
            self.plot_resolution = resolution

            self._update_log()

        if not self.no_plot:
            plt_distribution(density_dict, coordinates_dict,
                             labels=labels,
                             name=name,
                             limits=limits,
                             )


if __name__ == "__main__":
    print(sys.argv)

    truncate = None
    options = None
    step = None
    no_plot = True
    to_csv = True
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
