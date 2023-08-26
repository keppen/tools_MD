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
