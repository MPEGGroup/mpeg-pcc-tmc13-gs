import subprocess

from rich.console import Console

CONSOLE = Console(width=120)


def run_cmd(cmd: str, verbose: bool = True) -> subprocess.CompletedProcess:
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        CONSOLE.log(f"[red]Command failed (exit {e.returncode}):[/red] {cmd}")
        CONSOLE.log(f"[red]stdout:[/red]\n{e.stdout}")
        CONSOLE.log(f"[red]stderr:[/red]\n{e.stderr}")
        raise

    if verbose:
        CONSOLE.log(f"[red]stdout:[/red]\n{result.stdout}")
        CONSOLE.log(f"[red]stderr:[/red]\n{result.stderr}")

    return result
