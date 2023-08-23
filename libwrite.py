

def keep_duplicate(file_name):
    import os
    if os.path.exists(file_name):
        import datetime
        now = datetime.datetime.now()
        time = now.strftime("%Y-%m-%d_%H:%M:%S")
        # path, file_name = os.path.basename(file_path)
        basename, extension = file_name.split('.')
        os.rename(file_name, f'{basename}-{time}.{extension}')


def write_xyz(coordinates, names, file_name, **options):
    atoms_number = len(coordinates)
    step = options['step']
    with open(file_name, 'a') as file:
        msg = f"{atoms_number:>5}\n"
        file.write(msg)
        msg = f"STEP={step:>9}\n"
        file.write(msg)
        for i in range(atoms_number):
            i = str(i + 1)
            msg = f"{names[i]:>2}{coordinates[i][0]:14}{coordinates[i][1]:14}{coordinates[i][2]:14}\n"
            file.write(msg)


def write_pdb(content, file_name, **options):
    model = options['model']
    if model:
        model_count = 0
    file = open(file_name, 'w')
    msg = []
    for record in content:
        msg.append("TITLE Generated by keppen scripts\n")
        if model:
            model_count += 1
            msg.append(f'MODEL{model_count:>8}\n')
        for line in record:
            # print(line)
            msg.append(
                f"ATOM  "  # 'ATOM  '
                f"{line[0]:>5}"  # Atom serial number
                f'{" " * 2}'  # BLANK
                f"{line[1]:<4}"  # Atom name
                f"{line[2]}"  # Alternate location indicator
                f"{line[3]:>3}"  # Residue name
                f'{" " * 1}'  # BLANK
                f'{line[4]}'  # Chain location
                f"{line[5]:>4}"  # Residue sequence number
                f'{line[6]}'  # Code for insertion of residues
                f'{" " * 4}'  # BLANK
                f"{line[7]:>8.3f}"  # x Coord
                f"{line[8]:>8.3f}"  # y Coord
                f"{line[9]:>8.3f}"  # z Coord
                f"{line[10]:>6.2f}"  # Occupancy
                f"{line[11]:>6.2f}"  # Temperature factor
                f'{" " * 6}'   # BLANK
                f"{line[12]:<4}"    # Segment infdentifier
                f"{line[13]:>2}"    # Element symbol
                f"{line[14]:2}\n"   # Charge
            )
        msg.append("TER\n")
        if model:
            msg.append("ENDMDL\n")
    # print(msg)
    file.writelines(msg)
    file.close()


def write_log(name, message):
    with open(name, 'a') as file:
        file.write(message + '\n')
