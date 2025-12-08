from os import PathLike
from pathlib import Path

def file_exists(path: PathLike[str], non_empty: bool = True) -> bool:

    file_exist = Path(path).exists()
    is_file = Path(path).is_file()
    file_non_empty = Path(path).stat().st_size > 0

    if not is_file:
        return False
    if non_empty:
        return file_exist and file_non_empty
    else:
        return file_exist