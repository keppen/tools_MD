#!/home/mszatko/.conda/envs/code/bin/python3
import time
import math
import pandas as pd
import sys
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
import graph_mol
import re
from libread import read_pdb
from libmath import (
    calcDihedral,
    calcLinearRegression_PowerIteration,
    calcDistance_form_Vector,
    calcCentroid,
    calcPlane,
    calcAngle,
    calcVector,
    calcMiddlePoint,
)


class Visualize:
    def __init__(
        self, structure_pdb, cluster_pdb, options="all", truncate=None, step=None
    ):
        self.structure_pdb = structure_pdb
        self.cluster_pdb = cluster_pdb
        self.MOL = graph_mol.Molecule()
        self.truncate = int(truncate) if truncate else truncate
        self.options = "dhac" if options == "all" else options
        self.step = int(step) if step else step

        self.limit = True

        self.DF_pdb: pd.DataFrame
        self.DF_dihedrals: pd.DataFrame
        self.DF_hbond: pd.DataFrame
        self.DF_axis: pd.DataFrame
        self.DF_distance: pd.DataFrame

        self.runs_function = [
            self._collectDih,
            self._collectHbond,
            self._collectAxis,
            self._collectDistance,
        ]
        self.plot_function = [
            self._rama_plot,
            self._hbond_plot,
            self._axis_plot,
            self._contact_plot,
        ]
        self.dataframe_funtion = [
            self._dataframeDihedral,
            self._dataframeHbond,
            self._dataframeAxis,
            self._dataframeDistance,
        ]
        self.runs_bool = [True if c in self.options else False for c in "dhac"]

    def _initDicts(self):
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
            "OmegaPrim": self.omega1,
            "OmegaBis": self.omega2,
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

    def _initRun(self):
        """Run the calculation and generation of plots"""
        self._prepareMol()

        data = open(self.cluster_pdb, "r")
        self.pdb_gen = read_pdb(data)

        self._initDicts()
        # print(self.MOL.atoms)

        if self.step or self.truncate:
            self._runTruncated()
        else:
            self._runAll()
        self._mapRuns(self.dataframe_funtion)
        # print(self.DF_dihedrals, self.DF_hbond,
        #       self.DF_axis, self.DF_distance, sep="\n")
        self._mapRuns(self.plot_function)

    def _runAll(self):
        while self.limit:
            self._collectData()
            self._mapRuns(self.runs_function)

    def _runTruncated(self):
        while self.limit:
            self._collectData()

            if self.model == self.truncate:
                break
            if self.model % self.step != 0:
                continue

            print("MODEL: ", self.model, end="\r")
            self._mapRuns(self.runs_function)

    def _mapRuns(self, function_list):
        for condition, function in zip(self.runs_bool, function_list):
            if condition:
                function()

    def _splitName(self, file):
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
        print(f"CHECK THE FILE N?AME AND PICK CAREFULLY: {file}")
        print(f"indexes picked: \n\tcluser: {cluster}\n\tname: {name}")
        return name, cluster

    def _read_file(self, file):
        """Helper function"""
        with open(file, "r") as f:
            return f.readlines()

    def _prep_line(self, line: str):
        """Helper function"""
        return line.strip().split()

    def _prepareMol(self):
        """Helper function"""
        structure = self._read_file(self.structure_pdb)
        for line in structure:
            line = self._prep_line(line)
            if "ATOM" in line or "HETATM" in line:
                self.MOL.add_atom(line[1], line[11])
            if line[0] == "CONECT":
                self.MOL.add_bond(line[1:])
        self.MOL.run()
        self.cluster_pdb = self.MOL.updatePDB(self.cluster_pdb)

    def _printProgress(self, model):
        if model < 1000 and model % 200 == 0:
            print(f"Collecting model: {model}", end="\r")
        if model >= 1000 and model % 1000 == 0:
            print(f"Collecting model: {model}", end="\r")

    def _centerPoints(self):
        points = []
        # for index, row in self.DF_pdb.iterrows():
        #     points.append(np.array([row.x, row.y, row.z]))
        for index in self.MOL.backbone:
            p = self._searchByIndex(self.DF_pdb, int(index))
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

    def _collectData(self):
        try:
            self.model, pdb_data = next(self.pdb_gen)
            self.DF_pdb = pd.DataFrame.from_dict(pdb_data, orient="index").transpose()
            self._centerPoints()
        except StopIteration:
            self.limit = False

    def _collectDih(self):
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
            self.phi.append(self._calcPhi(self.DF_pdb, resid))
            self.psi.append(self._calcPsi(self.DF_pdb, resid))
            self.omega1.append(self._calcOmegaPrim(self.DF_pdb, resid))
            self.omega2.append(self._calcOmegaBis(self.DF_pdb, resid))

    def _collectGeometry(self):
        """Calculate dihedrals with data collected from pdb file.
        Dihedrals are stored in DataFrame object

        It iterates over all frames and (min + 1; max -1) residues indexes range
        It follows the same algorithm as _collectData method. Helper sunftions searches
        for coordianates and calculates diherdal angles.
        Truncating DataFrame by redundant data is not speeding up searching process
        THERE IS SOMETHING WRONG
        """
        max_resid = self.DF_pdb.max()["residue_seqnumber"]

        res0 = None
        res1 = None
        for resid in range(2, max_resid):
            res0 = res1 or self._calcPlane_and_Geomcentre(self.DF_pdb, resid)
            res1 = self._calcPlane_and_Geomcentre(self.DF_pdb, resid + 1)
            a = calcAngle(res0[0], res1[0])
            a = math.acos(a) / math.pi * 180
            line = calcVector(res0[1], res1[1])
            t1 = calcAngle(res0[0], line)
            t2 = calcAngle(res1[0], line)
            t1 = math.asin(abs(t1)) / math.pi * 180
            t2 = math.asin(abs(t2)) / math.pi * 180
            # t1 = t1 if t1 >= 0 else t1 + 180
            # t2 = t2 if t2 >= 0 else t2 + 180
            self.models_geometry.append(self.model)
            self.resids_geometry.append(resid)
            self.alpha.append(a)
            self.theta1.append(t1)
            self.theta2.append(t2)
            print(
                f"plane1: {res0[0]}\nplane2: {res1[0]}\nalpha: {a}\ntheta1: {t1}\ntheta2: {t2}"
            )

    def _collectHbond(self):
        NH_bondlenght = 0.99
        max_resid = self.DF_pdb.max()["residue_seqnumber"]

        for acceptor in range(2, max_resid + 1):
            for donor in range(1, max_resid):
                angle, lenght = self._calcHbond(self.DF_pdb, acceptor, donor)
                angle = math.acos(angle) / math.pi * 180
                self.models_hbond.append(self.model)
                self.resid_acceptor.append(acceptor)
                self.resid_donor.append(donor)
                self.angle.append(angle)
                self.lenght.append(lenght - NH_bondlenght)

    def _collectAxis(self):
        self.axis = calcLinearRegression_PowerIteration(self.points_backbone)
        for index in self.MOL.backbone:
            self.models_axis.append(self.model)
            point = self._searchByIndex(self.DF_pdb, int(index))
            self.distance_axis.append(calcDistance_form_Vector(point, self.axis))
            self.index.append(self.MOL.backbone.index(index))
            self.resids_axis.append(self.MOL._findbyID(index)["resid"])

    def _collectDistance(self):
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
                atom_1 = self._searchByIndex(self.DF_pdb, int(atom_id_1))
                atom_2 = self._searchByIndex(self.DF_pdb, int(atom_id_2))
                distance = np.linalg.norm(np.array(atom_1) - np.array(atom_2))
                self.models_distance.append(self.model)
                self.index_distance_1.append(index_1)
                self.index_distance_2.append(index_2)
                self.distance_distance.append(distance)

    def _dataframeDihedral(self):
        self.DF_dihedrals = pd.DataFrame.from_dict(
            self.dict_dihedral, orient="index"
        ).transpose()

    def _dataframeHbond(self):
        self.DF_hbond = pd.DataFrame.from_dict(
            self.dict_hbond, orient="index"
        ).transpose()
        self.pivot_hbond = pd.pivot_table(
            self.DF_hbond,
            values="angle",
            index="residue acceptor",
            columns="residue donor",
            aggfunc=np.mean,
        )

    def _dataframeAxis(self):
        self.DF_axis = pd.DataFrame.from_dict(
            self.dict_axis, orient="index"
        ).transpose()

    def _dataframeDistance(self):
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

    def _searchByType(self, df_model, resid, atom_type):
        """Helper function"""
        p = df_model.loc[
            # (df_model['model'] == model) &
            (df_model["residue_seqnumber"] == resid)
            & (df_model["atom_name"] == atom_type)
        ]
        return [float(p.iat[0, 7]), float(p.iat[0, 8]), float(p.iat[0, 9])]

    def _searchByIndex(self, df_model, index):
        """Helper function"""
        p = df_model.loc[
            # (df_model['model'] == model) &
            (df_model["atom_serialnumber"] == index)
        ]
        return [float(p.iat[0, 7]), float(p.iat[0, 8]), float(p.iat[0, 9])]

    def _calcPhi(self, df, resid):
        """Helper function"""
        p1 = self._searchByType(df, resid - 1, "C")
        p2 = self._searchByType(df, resid, "N")
        p3 = self._searchByType(df, resid, "CG")
        p4 = self._searchByType(df, resid, "CB")
        return calcDihedral(p1, p2, p3, p4)

    def _calcPsi(self, df, resid):
        """Helper function"""
        p1 = self._searchByType(df, resid, "CB")
        p2 = self._searchByType(df, resid, "OA")
        p3 = self._searchByType(df, resid, "C")
        p4 = self._searchByType(df, resid + 1, "N")
        return calcDihedral(p1, p2, p3, p4)

    def _calcOmegaPrim(self, df, resid):
        """Helper function"""
        p1 = self._searchByType(df, resid, "N")
        p2 = self._searchByType(df, resid, "CG")
        p3 = self._searchByType(df, resid, "CB")
        p4 = self._searchByType(df, resid, "OA")
        return calcDihedral(p1, p2, p3, p4)

    def _calcOmegaBis(self, df, resid):
        """Helper function"""
        p1 = self._searchByType(df, resid, "CG")
        p2 = self._searchByType(df, resid, "CB")
        p3 = self._searchByType(df, resid, "OA")
        p4 = self._searchByType(df, resid, "C")
        return calcDihedral(p1, p2, p3, p4)

    def _calcPlane_and_Geomcentre(self, df, resid):
        p0 = self._searchByType(df, resid - 1, "O")
        p1 = self._searchByType(df, resid - 1, "OA")
        p2 = self._searchByType(df, resid, "N")
        p3 = self._searchByType(df, resid - 1, "C")
        return calcPlane(p0, p1, p2), p3

    def _calcHbond(self, df, acceptor, donor):
        p0 = self._searchByType(df, acceptor, "N")
        p1 = self._searchByType(df, acceptor - 1, "C")
        p2 = self._searchByType(df, acceptor, "CG")
        p3 = self._searchByType(df, donor, "O")
        p4 = calcMiddlePoint(p1, p2)
        line1 = calcVector(p4, p0)
        line2 = calcVector(p0, p3)
        return calcAngle(line1, line2), np.linalg.norm(line2)

    def _hbond_plot(self):
        name, cluster = self._splitName(self.cluster_pdb)

        sb.set_context("paper", font_scale=1.35, rc={"lines.linewidth": 0.85})

        print(self.pivot_hbond)
        f, ax = plt.pyplot.subplots()
        sb.heatmap(
            self.pivot_hbond,
            cmap="crest",
        )
        f.savefig(
            f"hbond_{name}_{cluster}.png",
            dpi=100,
        )

        # plot.savefig(f"{'test1.png' if i == 0 else 'test2'}", dpi=1000)
        plt.pyplot.close()

    def _contact_plot(self):
        name, cluster = self._splitName(self.cluster_pdb)

        sb.set_context("paper", font_scale=1.35, rc={"lines.linewidth": 0.85})

        print(self.pivot_distance)
        f, ax = plt.pyplot.subplots()
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
        plt.pyplot.close()

    def _rama_plot(self, chirality=None):
        """Plots a Ramachandram plot from seaborn.JointGrid provided with pandas.DataFrame.
        3 fieds of grid:
                center: Ramachandran plot as kernel density plot of angle vs angle 2d plot
                without regard to resdue indexes.
                top and bottom: 1d kernel density plot of corresponding angle axis with
                regard to residue indexes.
        Chirality should be provided, it is presumed "S" isomers otherwise.
        """
        print("INFO: Starting plotting")
        name, cluster = self._splitName(self.cluster_pdb)
        if chirality is None:
            n = int(self.DF_dihedrals.max()["residue index"]) - 1
            chirality = ["S"] * n
            # print(chirality)

        letters = {
            "greek": ["\u03c9'", "\u03c9''", "\u03c8", "\u03c6"],
            "latin": ["omega1", "omega2", "psi", "phi"],
        }

        x = self.DF_dihedrals["Phi"]
        y = [
            self.DF_dihedrals["OmegaPrim"],
            self.DF_dihedrals["OmegaBis"],
            self.DF_dihedrals["Psi"],
        ]

        hue_order = [4, 3, 2]
        sb.set_context("paper", font_scale=1.35, rc={"lines.linewidth": 0.85})
        colors = ["black", "red", "blue"]

        for i in range(3):
            plot = sb.JointGrid(data=self.DF_dihedrals, space=0)
            joint = sb.kdeplot(
                data=self.DF_dihedrals,
                x=x,
                y=y[i],
                fill=True,
                bw_adjust=0.05,
                cmap="Greys",
                clip=[-180, 180],
                thresh=0,
                ax=plot.ax_joint,
                levels=100,
            )
            top = sb.kdeplot(
                data=self.DF_dihedrals,
                x=x,
                hue="residue index",
                hue_order=hue_order,
                fill=False,
                ax=plot.ax_marg_x,
                bw_adjust=0.2,
                legend=True,
                palette=colors,
                clip=[-180, 180],
            )
            sb.kdeplot(
                data=self.DF_dihedrals,
                y=y[i],
                hue="residue index",
                fill=False,
                ax=plot.ax_marg_y,
                bw_adjust=0.1,
                legend=False,
                palette=colors,
                clip=[-180, 180],
            )
            plot.ax_marg_x.legend(
                labels=[
                    f"2 ({chirality[0]})",
                    f"3 ({chirality[1]})",
                    f"4 ({chirality[2]})",
                ],
                title="No. Mer",
                frameon=False,
            )
            joint.text(x=-240, y=262, s=f"{name.upper()}-{cluster}")
            sb.move_legend(top, loc="upper right", bbox_to_anchor=(1.28, 1.5))
            plot.ax_marg_x.set_xlim(-180, 180)
            plot.ax_marg_y.set_ylim(-180, 180)
            plot.ax_joint.set_xticks([-120, -60, 0, 60, 120])
            plot.ax_joint.set_yticks([-120, -60, 0, 60, 120])
            joint.set(xlabel=letters["greek"][3], ylabel=letters["greek"][i])
            plot.savefig(
                f'{letters["latin"][i]}_{name}_{cluster}.png',
                dpi=500,
            )

            # plot.savefig(f"{'test1.png' if i == 0 else 'test2'}", dpi=1000)
            plt.pyplot.close()
            print(
                f'INFO: {letters["latin"][i]}_{name}_{cluster}.png has been generated'
            )

    def _geom_plot(self):
        name, cluster = self._splitName(self.cluster_pdb)

        letters = {
            "greek": ["\u03b1", "\u03b8\u2081", "\u03b8\u2082"],
            "latin": ["alpha", "theta1", "theta2"],
        }

        sb.set_context("paper", font_scale=1.35, rc={"lines.linewidth": 0.85})
        colors = ["black", "red", "blue"]

        fig, (ax1, ax2, ax3) = plt.pyplot.subplots(3, 1, figsize=(7, 5), sharex=True)
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

    def _axis_plot(self):
        name, cluster = self._splitName(self.cluster_pdb)

        plot = sb.lineplot(data=self.DF_axis, x="index", y="distance", marker="o")
        plot.figure.savefig(f"axis_{name}_{cluster}.png", dpi=500)
        plt.close()


if __name__ == "__main__":
    print(sys.argv)
    start_time = time.time()
    try:
        truncate = sys.argv[3]
        options = sys.argv[4]
        step = sys.argv[5]
    except IndexError:
        truncate = None
        options = None
        step = None
    R = Visualize(
        sys.argv[1], sys.argv[2], truncate=truncate, step=step, options=options
    )
    R._initRun()

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
