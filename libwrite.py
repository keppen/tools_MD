

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
