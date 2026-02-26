def get_version():
    import tomllib
    from pathlib import Path
    dir = Path(__file__).resolve().parent.parent
    pyproject_path = f"{dir}/pyproject.toml"
    with open(pyproject_path, "rb") as f:
        version = tomllib.load(f)["project"]["version"]
    return version


VERSION = get_version()
MOTIF_COLUMN_NAME = 'motif'
MOTIF_COLUMN_ID = 'name'
MOTIF_COLUMN_N_MODS = 'n_modules'
MOTIF_COLUMN_MOD_CLASS = 'module_classes'

DANTE_DESCRIPTION = '''
DANTE = Da Amazing NucleoTide Exposer (Remastered)
--------------------------------------------------
'''
MAX_REPETITIONS = 40
