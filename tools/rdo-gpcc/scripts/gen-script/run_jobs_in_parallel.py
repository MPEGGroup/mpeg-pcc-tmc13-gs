import concurrent.futures
import subprocess
import sys
from dataclasses import dataclass
from multiprocessing import get_context
from pathlib import Path

import rich
import tyro
from rich.progress import BarColumn, MofNCompleteColumn, Progress, TextColumn, TimeElapsedColumn, TimeRemainingColumn

# CONSOLE = rich.console.Console(width=120)
CONSOLE = rich.console.Console(width=120, force_terminal=True)


@dataclass
class Config:
    job_file: Path
    log_dir: Path
    n_workers: int
    exit_on_fail: bool = False
    skip_completed: bool = False


@dataclass
class Status:
    UNKNOWN = "UNKNOWN"
    RUNNING = "RUNNING"
    FAILED = "FAILED"
    SKIPPED = "SKIPPED"
    COMPLETED = "COMPLETED"


def run_job(cmd: str, log_dir: Path, skip_completed: bool = False):
    log_file = log_dir / f"{str(Path(cmd).stem)}.log"
    status_file = log_dir / f"{str(Path(cmd).stem)}.status"
    details = {"Job": cmd, "Log file": str(log_file), "Status": Status.UNKNOWN}

    if skip_completed and status_file.exists():
        status = status_file.read_text()
        if status == Status.COMPLETED:
            details["Status"] = Status.SKIPPED
            return 0, details

    status_file.write_text(Status.RUNNING)

    with open(log_file, "w") as out:
        process = subprocess.Popen(cmd, shell=True, stdout=out, stderr=out)
        process.wait()

    if process.returncode == 0:
        status_file.write_text(Status.COMPLETED)
        details["Status"] = Status.COMPLETED
    else:
        status_file.write_text(Status.FAILED)
        details["Status"] = Status.FAILED

    return process.returncode, details


def main(config: Config):
    CONSOLE.log(config)
    config.log_dir.mkdir(parents=True, exist_ok=True)

    CONSOLE.log("Reading job file.")
    with config.job_file.open("r") as f:
        jobs = [line.strip() for line in f if line.strip()]
    n_jobs = len(jobs)
    CONSOLE.log(f"# of jobs: {n_jobs}")

    CONSOLE.log("Executing jobs in parallel.")
    details_list = []
    with (
        concurrent.futures.ProcessPoolExecutor(
            max_workers=config.n_workers, mp_context=get_context("spawn")
        ) as executor,
        Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=CONSOLE,
            speed_estimate_period=3600,
        ) as progress,
    ):
        # Submit jobs.
        futures = {executor.submit(run_job, job, config.log_dir, config.skip_completed): job for job in jobs}
        success = True
        # Wait and check results.
        for future in progress.track(concurrent.futures.as_completed(futures), total=n_jobs):
            returncode, details = future.result()
            details_list.append(details)
            if returncode != 0:
                CONSOLE.log("Job failed:", details)
                if config.exit_on_fail:
                    CONSOLE.log("Cancelling jobs.")
                    success = False
                    for process in executor._processes.values():
                        process.kill()
                    break

    CONSOLE.log("Printing summary.")
    n_failed = 0
    n_completed = 0
    n_skipped = 0
    for details in details_list:
        if details["Status"] == Status.FAILED:
            n_failed += 1
        elif details["Status"] == Status.COMPLETED:
            n_completed += 1
        elif details["Status"] == Status.SKIPPED:
            n_skipped += 1
    CONSOLE.log(f"- FAILED: {n_failed}")
    CONSOLE.log(f"- COMPLETED: {n_completed}")
    CONSOLE.log(f"- SKIPPED: {n_skipped}")
    CONSOLE.log(f"- CANCELLED: {n_jobs - n_failed - n_completed - n_skipped}")

    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main(tyro.cli(Config))
