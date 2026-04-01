def get_version():
    import tomllib
    from pathlib import Path
    dir = Path(__file__).resolve().parent.parent
    pyproject_path = f"{dir}/pyproject.toml"
    with open(pyproject_path, "rb") as f:
        version = tomllib.load(f)["project"]["version"]
    return version


VERSION = get_version()

MOTIF_COLUMN_ID = "name"
MOTIF_COLUMN_NAME = "motif"
MOTIF_COLUMN_READ_ID = "read_id"
MOTIF_COLUMN_MODULES = "modules"
MOTIF_COLUMN_N_MODS = "n_modules"
MOTIF_COLUMN_MOD_CLASS = "module_classes"
MOTIF_COLUMN_MISMATCHES_STR = "mismatches_str"
MOTIF_COLUMN_MODULE_REPETITIONS = "module_repetitions"
MOTIF_COLUMN_MODULE_NOMENCLATURES = "module_nomenclatures"
MOTIF_COLUMN_MODULE_SEQUENCES = "module_sequences"

DANTE_DESCRIPTION = '''
DANTE = Da Amazing NucleoTide Exposer (Remastered)
--------------------------------------------------
'''
MAX_REPETITIONS = 40
