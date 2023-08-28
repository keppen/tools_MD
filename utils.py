greek_letters = {
        "phi": "\u03c6",
        "psi": "\u03c8",
        "xi": "\u03be",
        "chi": "\u03c7",
        "alpha": "\u03b1",
        "theta": "\u03b8",
        }
other_characters = {
        "_1": "\u2081",
        "_2": "\u2082",
        }


def latin_2_greek(symbols):
    import re
    new_list = []
    for symbol in symbols:
        symbol = symbol.lower()
        if bool(re.search(r"\d", symbol)):
            symbol1 = symbol[:-1]
            subscript = f"_{symbol[-1]}"
            new_list.append(f"{greek_letters[symbol1]}{other_characters[subscript]}")
        else:
            new_list.append(greek_letters[symbol])
    return new_list


def split_name(file):
    """Split name of file and get the information about analyzed file.

    Parameters:
        file(str): a file to be analyzed. e.g. 10ns-a1.pdb001.pdb

    Returns:
        name(str): a na of a system
        cluster(str): a number of a cluster
    """
    import re

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
