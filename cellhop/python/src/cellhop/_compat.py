import os
from typing import Optional


def get_r_home_from_env(default: Optional[str] = None) -> Optional[str]:
    return os.environ.get("R_HOME", default)
