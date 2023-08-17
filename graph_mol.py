#!/home/keppen/.conda/envs/code/bin/python3
import sys
from libread import read_pdb, dict_2_list
from libwrite import write_pdb


def prep_line(line):
    return line.strip().split()


class Molecule:
    def __init__(self, name="mol"):
        self.name = name
        self.atoms = []
        self.CT = None
        self.NT = None
        self.Ns = []
        self.backbone = []

    def add_atom(self, index, atom):
        """Describing an atom with properties derived from pdb molecule file.
        Setting atoms as nodes for grapth representation.

        Parameters:
            index (int): index of atom in pdb file
            atom (str): an type of element, NOT an atomtype

        Returns: None
        """
        d = {
            "id": int(index),
            "element": atom,
            "type": "",
            "bonds": set(),
            "visited": False,
            "resid": None,
        }
        self.atoms.append(d)

    def add_bond(self, bond):
        """Building the molecule by connecting atoms (nodes) with bonds (edges)
        Setting up a molecule as a graph.

        Parameters:
            bond (list): a list of integers, list[0] an atom to which list[1:]
            are connected to.

        Returns:: None
        """
        for i in bond[1:]:
            for atom in self.atoms:
                if atom["id"] == int(bond[0]):
                    atom["bonds"].add(int(i))

    def run(self):
        """Run the methods of the class"""

        self._findN()
        self._findCT()
        self._findNT()
        self._parse()

    def _parse(self):
        """Parsing backbone and updating its parameters.
        Updating: type, resid
        """
        resid = 1
        for backatom in self.backbone:
            backatom = self._findbyID(backatom)
            if backatom["type"] == "N":
                resid += 1
            # print(backatom)
            self._update_type(backatom)
            backatom["resid"] = resid
            search = self._BFS(backatom, border=self.backbone)
            try:
                while True:
                    srch = next(search)[0]
                    # print(backatom['id'], srch['id'])
                    self._update_type(srch)
                    srch["resid"] = resid
            except StopIteration:
                # print('Iteration stopped')
                pass
            finally:
                del search

    def updateDataPDB(self, line):
        line = list(line)
        if line[3] == "BPA":
            return None
        atom = self._findbyID(line[0])
        line[1] = atom["type"]
        line[3] = "UNK"
        line[5] = atom["resid"]
        return line

    def newPDB(self, cluster_pdb):
        cluster_pdb_content = open(cluster_pdb, "r")
        generator = read_pdb(cluster_pdb_content)

        # pdb_file = [r[1].values() for r in generator]
        pdb_file = []
        for record in generator:
            record = record[1]
            pdb_file.append(dict_2_list(record))

        print("Updating data")

        for i1, record in enumerate(pdb_file):
            for i2, line in enumerate(record):
                pdb_file[i1][i2] = self.updateDataPDB(line)

        file_name = f"NEW_{cluster_pdb}"
        print(f"CREATING: {file_name}")

        # print(*pdb_file, sep=3*'\n')

        write_pdb(pdb_file, file_name, model=True)

        print(f"NEW_{cluster_pdb} has been created. Check it for any problems")
        return f"NEW_{cluster_pdb}"

    def printout(self):
        """A debug method"""
        for atom in self.atoms:
            if atom["element"]:
                print(
                    atom["id"], atom["element"], atom["type"], atom["resid"], sep="\t"
                )

    def _update_type(self, atom):
        """Higher level function to disciminate an atomtype of atom.

        Parameters:
            atom (dict): a dictionary of parameters of an atom
                            set up in add_atom method
        Returns: None
        """
        # print(self.CT, self.NT, self.backbone, sep='\n')
        # atom = self._findbyID(atom)
        # print(atom)
        if atom["element"] == "C":
            # print('in "C" type')
            self._typeC(atom)
        if atom["element"] == "N":
            # print('in "N" type')
            self._typeN(atom)
        if atom["element"] == "O":
            # print('in "O" type')
            self._typeO(atom)
        if atom["element"] == "H":
            atom["type"] = "H"

    def _findN(self):
        """Find N elements. Those are starting points for BFS walker
        in backbone finding task"""
        for atom in self.atoms:
            if atom["element"] == "N":
                self.Ns.append(atom["id"])

    def _findCT(self):
        """Find CT atomtype. This is C-terminus of macromolecule.
        This is an engind point of BFS walker in backbone finding task

        1. Check if O element
        2. Check if sp3
        3. Starts BFS instance. Sets N atoms as borders.
            Count O elements within borders. C-terminus should be
            equal to 2.
        It presumes no other OH groups in side chains at the termini monomers
        """
        for atom in self.atoms:
            if atom["element"] != "O":
                continue
            if len(atom["bonds"]) != 2:
                continue
            Ox = 0
            for bond in atom["bonds"]:
                bond = self._findbyID(bond)
                search = self._BFS(bond, border=self.Ns)
                try:
                    while True:
                        srch = next(search)[0]
                        if srch["element"] == "O":
                            Ox += 1
                except StopIteration:
                    pass
                finally:
                    del search
                    self._refresh()
            if Ox == 2:
                atom["type"] = "OT"
                self.CT = atom

    def _findNT(self):
        """Find NT atomtype. This is N-terminus of macromolecule.
        Define backbone of macromolecule.
        Found by BTS walker

        1. Check if N element
        2. Start first instances of BFS starting from bonded atoms.
           Count N elements within border defined N elements.
            The count of N elements has to be equeal to 4.
        3. Start second instances if BSF.
            Get path from NT candidate to CT atom.
        4. Sat candidate with longest path as NT.

        It presumes there is no N atoms in protecting group. First resudie
        has no additional N atoms.
        """
        candidates = []
        for atom in self.atoms:
            if atom["element"] != "N":
                continue

            N = 0
            for bond in atom["bonds"]:
                bond = self._findbyID(bond)
                search = self._BFS(bond, border=self.Ns)
                try:
                    while True:
                        srch = next(search)[0]
                        if srch["element"] == "N":
                            N += 1
                except StopIteration:
                    pass
                finally:
                    del search
                    self._refresh()

            if N == 4:
                path = self._BFS(atom)
                while True:
                    p = next(path)
                    if p[0]["type"] == "OT":
                        p[0] = atom
                        candidates.append(p)
                        del path
                        self._refresh()
                        break

        d = [x[1] for x in candidates]
        d = list(map(len, d))
        self.NT = candidates[d.index(max(d))][0]
        self.NT["type"] = "NT"
        self.backbone = candidates[d.index(max(d))][1]
        return 0

    def _countELE(self, atom):
        """Helper funtion to distiguish the elements
        It counts which and how many elements are bond to atom"""
        h = 0
        c = 0
        n = 0
        o = 0
        for bond in atom["bonds"]:
            bond = self._findbyID(bond)
            # print('\t', bond)
            if bond["element"] == "H":
                h += 1
            elif bond["element"] == "C":
                c += 1
            elif bond["element"] == "N":
                n += 1
            elif bond["element"] == "O":
                o += 1
            else:
                # print('error:', bond)
                return False
        return h, c, n, o

    def _typeC(self, atom):
        if not atom["type"] and atom["id"] in self.backbone:
            if len(atom["bonds"]) == 3:
                h, c, n, o = self._countELE(atom)
                # print(h, c, n, o)
                if o == 2 and n == 1:
                    # print('===TYPE C===')
                    atom["type"] = "C"
            if len(atom["bonds"]) == 4:
                h, c, n, o = self._countELE(atom)
                # print(h, c, n, o)
                if n == 1:
                    # print('===TYPE CG===')
                    atom["type"] = "CG"
                if o == 1:
                    # print('===TYPE CB===')
                    atom["type"] = "CB"
        elif not atom["type"]:
            atom["type"] = "CZ"

    def _typeN(self, atom):
        if not atom["type"]:
            if len(atom["bonds"]) == 3:
                h, c, n, o = self._countELE(atom)
                if c >= 2:
                    # print('===TYPE N===')
                    atom["type"] = "N"

    def _typeO(self, atom):
        if not atom["type"]:
            bonds = list(atom["bonds"])
            if atom["id"] in self.backbone:
                atom["type"] = "OA"
            elif len(bonds) == 1 and bonds[0] in self.backbone:
                atom["type"] = "O"
            else:
                atom["type"] = "OX"

    def _refresh(self):
        """Helper function. Resets the state of graph to run BSF"""
        for atom in self.atoms:
            atom["visited"] = False

    def _findbyID(self, id):
        """Helper function. Searches list of dictionaries to find dict by id value"""
        for atom in self.atoms:
            if atom["id"] == id:
                return atom

    def _BFS(self, node, border=list()):
        """Breath First Search walker. A generator.
        A border is an atom not added to queue,
        it makes walker stop searching beyond it.

        Parameters:
            node (dict): an atom
            border: an list of atoms

        Yields:
            node_path (list): [atom, [list of atoms ids]] -> [node, path]
        """
        path = [node["id"]]
        node_path = [node, path]
        queue = [node_path]
        while queue:
            node, path = queue.pop(0)
            node["visited"] = True
            for neighbor in node["bonds"]:
                neighbor = self._findbyID(neighbor)
                node_path = [neighbor, path + [neighbor["id"]]]
                if not neighbor["visited"]:
                    if neighbor["id"] in border:
                        # print('found border', neighbor)
                        neighbor["visited"] = True
                        yield node_path
                        continue
                    yield node_path

                    queue.append(node_path)


if __name__ == "__main__":
    start_pdb = sys.argv[1]
    cluster_pdb = sys.argv[2]
    with open(start_pdb, "r") as file:
        file = file.readlines()
    mol = Molecule()
    for line in file:
        line = prep_line(line)
        if "ATOM" in line or "HETATM" in line:
            mol.add_atom(line[1], line[11])
        if line[0] == "CONECT":
            mol.add_bond(line[1:])
    mol.run()
    mol.newPDB(cluster_pdb)
