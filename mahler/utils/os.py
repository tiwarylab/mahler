from os import PathLike
from pathlib import Path

def file_exists(path: PathLike[str], non_empty: bool = True) -> bool:

    file_exist = Path(path).exists()
    is_file = Path(path).is_file()

    if not is_file:
        return False
    if non_empty:
        file_non_empty = Path(path).stat().st_size > 0
        return file_exist and file_non_empty
    else:
        return file_exist
    
def list_files(directory: PathLike[str], non_empty: bool = True) -> list[Path]:
    """Return a sorted list of files in a directory matching a given pattern. 
    Ignore directories."""

    dir_path = Path(directory)
    if not dir_path.is_dir():
        raise NotADirectoryError(f"The path {directory} is not a valid directory.")
    items = list(dir_path.glob("*"))
    files = [item for item in items if file_exists(item, non_empty=non_empty)]  

    return sorted(files)