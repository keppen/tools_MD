    def _ramachandran_depreciated(self, chirality=None):
        """Plots a Ramachandram plot from seaborn.JointGrid provided with pandas.DataFrame.
        3 fieds of grid:
                center: Ramachandran plot as kernel density plot of angle vs angle 2d plot
                without regard to resdue indexes.
                top and bottom: 1d kernel density plot of corresponding angle axis with
                regard to residue indexes.
        Chirality should be provided, it is presumed "S" isomers otherwise.
        """
        import seaborn as sb
        import matplotlib.pyplot as plt
        print("INFO: Starting plotting")
        name, cluster = self._split_name(self.cluster_pdb)
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
            self.DF_dihedrals["Xi"],
            self.DF_dihedrals["Chi"],
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
            plt.close()
            print(
                f'INFO: {letters["latin"][i]}_{name}_{cluster}.png has been generated'
            )


