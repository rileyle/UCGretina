#!/usr/bin/env python3
"""UCGretina testing suite: functional tests and performance benchmarks."""

import argparse
import datetime
import json
import math
import os
import re
import subprocess
import sys
import shutil

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS_DIR = os.path.join(PROJECT_ROOT, "tests")
#MACROS_DIR = os.path.join(TESTS_DIR, "macros")
TMP_DIR = os.path.join(TESTS_DIR, "tmp")
BASELINES_FILE = os.path.join(TESTS_DIR, "baselines.json")
BENCHMARK_LOG = os.path.join(PROJECT_ROOT, "benchmark.log")

# # Geant4 data environment variable names and their expected subdirectory names
# _G4_DATA_VARS = {
#     "G4NEUTRONHPDATA":   "G4NDL4.6",
#     "G4LEDATA":          "G4EMLOW7.13",
#     "G4LEVELGAMMADATA":  "PhotonEvaporation5.7",
#     "G4RADIOACTIVEDATA": "RadioactiveDecay5.6",
#     "G4PARTICLEXSDATA":  "G4PARTICLEXS3.1.1",
#     "G4PIIDATA":         "G4PII1.3",
#     "G4REALSURFACEDATA": "RealSurface2.2",
#     "G4SAIDXSDATA":      "G4SAIDDATA2.0",
#     "G4ABLADATA":        "G4ABLA3.1",
#     "G4INCLDATA":        "G4INCL1.0",
#     "G4ENSDFSTATEDATA":  "G4ENSDFSTATE2.3",
# }

# _G4_DATA_SEARCH_ROOTS = [
#     "/opt/geant4-v10.7.4/share/Geant4-10.7.4/data",
#     "/usr/local/geant4/geant4-v10.7.4/share/Geant4-10.7.4/data",
#     "/usr/share/geant4/data",
# ]


# def _fix_geant4_data_paths():
#     """Override G4 data env vars if they point to non-existent directories.

#     The installed geant4.sh may have been built for a different prefix.
#     We search known root paths for the correct data directories.
#     """
#     # Find a data root that actually exists
#     data_root = None
#     for root in _G4_DATA_SEARCH_ROOTS:
#         if os.path.isdir(root):
#             data_root = root
#             break

#     if data_root is None:
#         return  # Can't fix; leave env as-is

#     for var, subdir in _G4_DATA_VARS.items():
#         current = os.environ.get(var, "")
#         if current and os.path.isdir(current):
#             continue  # Already valid
#         candidate = os.path.join(data_root, subdir)
#         if os.path.isdir(candidate):
#             os.environ[var] = candidate


# _fix_geant4_data_paths()


def find_binary(name):
    """Locate a UCGretina binary in $G4WORKDIR/bin/$G4SYSTEM/."""
    g4workdir = os.environ.get("G4WORKDIR")
    g4system = os.environ.get("G4SYSTEM")
    if not g4workdir or not g4system:
        print("ERROR: G4WORKDIR and G4SYSTEM must be set.", file=sys.stderr)
        sys.exit(1)
    path = os.path.join(g4workdir, "bin", g4system, name)
    if not os.path.isfile(path):
        print(f"ERROR: Binary not found: {path}", file=sys.stderr)
        sys.exit(1)
    return path


def find_binary_optional(name):
    """Locate a UCGretina binary; return None if not found (no abort)."""
    g4workdir = os.environ.get("G4WORKDIR")
    g4system = os.environ.get("G4SYSTEM")
    if not g4workdir or not g4system:
        return None
    path = os.path.join(g4workdir, "bin", g4system, name)
    return path if os.path.isfile(path) else None


def get_project_root():
    """Return the absolute path to the project root."""
    return PROJECT_ROOT


def setup_workdir(test_name, example_path, support_files):
    """Create a per-test working directory with support files, rereating
    it if it exists.
    Returns the working directory path.
    """
    workdir = os.path.join(TMP_DIR, test_name)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir, exist_ok=True)
    for file in support_files:
        src = os.path.join(PROJECT_ROOT, os.path.dirname(example_path))
        shutil.copy(os.path.join(src, file), workdir)
    return workdir


def write_base_macro(base_macro_path, example_path, output_command, workdir):
    """Write the body of the macro file without the /run/beamOn N"""

    lines = ""
    with open(example_path, "r") as f:
        lines = f.readlines()
    
    with open(os.path.join(workdir, base_macro_path), "w") as f:
        for line in lines:
            # Fix CAD model path
            if "/CADModelPath" in line:
                line = "/ScanningTable/CADModelPath ../../../cadModels\n"
            # Fix the cache filenames
            if "/Cache/Output" in line:
                line = "/Cache/Output cache_gen.cache\n"
            if "/Cache/Input" in line:
                line = "/Cache/Input cache_gen.cache\n"
            # Omit the output file and beamOn commands.
            if ("/Output/Filename" not in line) \
               and ("/Mode2/Filename" not in line) \
               and ("/run/beamOn" not in line):
                f.write(line)
        f.write(output_command)


def write_run_macro(base_macro_path, n_events, wrapper_path):
    """Write a wrapper macro that executes base_macro then /run/beamOn N."""
    with open(wrapper_path, "w") as f:
        f.write(f"/control/execute {base_macro_path}\n")
        f.write(f"/run/beamOn {n_events}\n")

def run_sim(binary, macro_path, workdir):
    """Run a simulation, return (stdout, stderr, returncode)."""
    result = subprocess.run(
        [binary, macro_path],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    return result.stdout, result.stderr, result.returncode


def parse_events_per_sec(stdout):
    """Extract events/sec from the RunAction end-of-run summary line.

    The final summary line looks like:
      Real time: ...   System time: ...   User time: ...   4523 events/s
    We use the last match to skip intermediate progress-update values.
    Returns float or None.
    """
    matches = re.findall(r"([\d.]+)\s+events/s", stdout)
    if matches:
        return float(matches[-1])
    return None


def check_fatal(stdout, stderr):
    """Return True if fatal error indicators are present in output."""
    fatal_patterns = [
        "Fatal Exception",
        "Segmentation fault",
        "FatalException",
        "G4Exception : Fatal",
    ]
    combined = stdout + stderr
    return any(p in combined for p in fatal_patterns)


def count_lines(filepath):
    """Return line count of a file using wc -l."""
    result = subprocess.run(["wc", "-l", filepath], capture_output=True, text=True)
    if result.returncode != 0:
        return 0
    return int(result.stdout.strip().split()[0])


def load_baselines():
    """Load baselines.json; return {} if absent.

    Strips the '_meta' key so callers only see test-name → count entries.
    """
    if not os.path.isfile(BASELINES_FILE):
        return {}
    with open(BASELINES_FILE) as f:
        data = json.load(f)
    data.pop("_meta", None)
    return data


def save_baselines(data):
    """Write baselines.json with pretty-print JSON, including provenance metadata."""
    git_hash, git_branch = get_git_info()
    cpu = get_cpu_info()
    out = {
        "_meta": {
            "git_hash": git_hash,
            "git_branch": git_branch,
            "cpu": cpu,
        }
    }
    out.update(data)
    with open(BASELINES_FILE, "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")


def check_baseline(name, observed, baselines):
    """Check observed line count against baseline with 2*sqrt(N) tolerance.

    Returns (passed: bool, message: str).
    If no baseline exists, sets it and returns (True, '[BASELINE SET] ...').
    """
    if name not in baselines:
        baselines[name] = observed
        return True, f"[BASELINE SET] {name}: {observed} lines"

    baseline = baselines[name]
    tolerance = 2 * math.sqrt(baseline)
    if abs(observed - baseline) <= tolerance:
        return True, (f"[PASS] {name:<30} output lines={observed}  "
                      f"baseline={baseline}  tolerance=±{tolerance:.0f}")
    else:
        return False, (f"[FAIL] {name:<30} output lines={observed}  "
                       f"baseline={baseline}  tolerance=±{tolerance:.0f}  "
                       f"** out of range **")


def get_git_info():
    """Return (6-char hash, branch) from git."""
    try:
        hash_ = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=PROJECT_ROOT, text=True
        ).strip()[:6]
        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            cwd=PROJECT_ROOT, text=True
        ).strip()
    except subprocess.CalledProcessError:
        hash_, branch = "unknown", "unknown"
    return hash_, branch


def get_cpu_info():
    """Return a compact CPU identifier string.

    Tries, in order:
      1. /proc/cpuinfo 'model name' (Linux)
      2. sysctl machdep.cpu.brand_string (macOS)
      3. hostname as a last resort
    """
    # Linux
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass

    # macOS
    try:
        result = subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"],
            capture_output=True, text=True
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except FileNotFoundError:
        pass

    import socket
    return socket.gethostname()


def append_benchmark_log(rows):
    """Append rows to benchmark.log (TSV). Write header if file is new."""
    header = "date\tgit_hash\tgit_branch\tcpu\tvariant\tevents\tevents_per_sec\n"
    write_header = not os.path.isfile(BENCHMARK_LOG)
    with open(BENCHMARK_LOG, "a") as f:
        if write_header:
            f.write(header)
        for row in rows:
            f.write("\t".join(str(v) for v in row) + "\n")


def main():
    parser = argparse.ArgumentParser(description="UCGretina test and benchmark driver")
    parser.add_argument("--mode", required=False, default=None,
                        choices=["smoke", "sources", "inbeam", "scanning",
                                 "background", "benchmark"],
                        help="Test mode to run")
    parser.add_argument("--events", type=int, default=10000,
                        help="Event count for benchmark mode (default: 10000)")
    parser.add_argument("--update-baselines", action="store_true",
                        help="Reset all baselines to current observed values")
    args = parser.parse_args()

    if args.update_baselines:
        update_baselines()
        return

    if args.mode is None:
        parser.error("--mode is required unless --update-baselines is specified")

    os.makedirs(TMP_DIR, exist_ok=True)

    if args.mode == "smoke":
        run_smoke()
    elif args.mode == "sources":
        run_functional("sources")
    elif args.mode == "inbeam":
        run_functional("inbeam")
    elif args.mode == "scanning":
        run_functional("scanning")
    elif args.mode == "background":
        run_functional("background")
    elif args.mode == "benchmark":
        run_benchmark(args.events)


SMOKE_EVENTS = 100

# (test_name, binary_name, macro_file, geometry_prefix)
SMOKE_CASES = [
    ("smoke_standard",  "UCGretina",      "func_inbeam_standard.mac",
     "examples/inbeam/fit/s44_1329.mac",
     ["aclust", "aeuler", "aslice", "asolid", "awalls", 
      "crmat.LINUX", "z16.a44.lvldata"]),
    ("smoke_lh",        "UCGretina_LH",   "func_inbeam_lh.mac",
     "examples/inbeam/fitLH/s44_1329.mac",
     ["aclust", "aeuler", "aslice", "asolid", "awalls", 
      "crmat.LINUX", "z16.a44.lvldata"]),
    ("smoke_pol",       "UCGretina_Pol",  "func_inbeam_pol.mac",
     "examples/inbeam/angdist/s44_1329.mac",
     ["aclust", "aeuler", "aslice", "asolid", "awalls", 
      "crmat.LINUX", "z16.a44.lvldata"]),
    ("smoke_scan",      "UCGretina_Scan", "func_scanning.mac",
     "examples/scan/scan.mac",
     ["aclust", "aeuler", "aslice", "asolid", "awalls"]),
    ("smoke_sources",   "UCGretina",      "func_sources_eu152.mac",
     "examples/sources/eu152/eu152.mac",
     ["aclust", "aeuler", "aslice", "asolid", "awalls",
      "crmat.LINUX", "z62.a152.lvldata", "z64.a152.lvldata"]),
    ("smoke_background","UCGretina",      "func_background.mac",
     "examples/background/background.mac",
     ["aclust", "aeuler", "aslice", "asolid", "awalls", 
      "crmat.LINUX"]),
]


def run_smoke():
    print(f"\n=== test-smoke ({SMOKE_EVENTS} events) ===")
    failures = 0
    for test_name, binary_name, macro_file, example_path, support_files in SMOKE_CASES:
        binary = find_binary_optional(binary_name)
        if binary is None:
            print(f"[SKIP] {test_name:<30} binary {binary_name} not found")
            continue
        workdir = setup_workdir(test_name, example_path, support_files)
        write_base_macro(macro_file, example_path,
                         "/Output/Filename output.out", workdir)
        wrapper = os.path.join(workdir, "run.mac")

        write_run_macro(macro_file, SMOKE_EVENTS, wrapper)

        stdout, stderr, returncode = run_sim(binary, wrapper, workdir)

        ok, msg = _check_run_criteria(test_name, stdout, stderr, returncode)
        print(msg)
        if not ok:
            failures += 1

    if failures:
        print(f"\n{failures} FAILED")
        sys.exit(1)
    else:
        print("\nAll smoke tests passed.\n")


def _check_run_criteria(test_name, stdout, stderr, returncode):
    """Check exit code, no fatals, end-of-run line, events/sec > 0."""
    if returncode != 0:
        return False, f"[FAIL] {test_name:<30} exit code {returncode}"
    if check_fatal(stdout, stderr):
        return False, f"[FAIL] {test_name:<30} fatal error in output"
    eps = parse_events_per_sec(stdout)
    if eps is None:
        return False, f"[FAIL] {test_name:<30} end-of-run line not found"
    if eps <= 0:
        return False, f"[FAIL] {test_name:<30} events/sec = {eps}"
    return True, f"[PASS] {test_name:<30} {eps:.0f} events/s"


FUNCTIONAL_EVENTS = 1000

# Maps mode -> list of (test_name, binary_name, macro_file, example_path, support_files, output_filename)
FUNCTIONAL_CASES = {
    "sources": [
        ("sources_eu152",  "UCGretina", "func_sources_eu152.mac",
         "examples/sources/eu152/eu152.mac",
         ["aclust", "aeuler", "aslice", "asolid", "awalls", 
          "crmat.LINUX", "z62.a152.lvldata", "z64.a152.lvldata"],
         "output.out"),
        ("sources_co60",   "UCGretina", "func_sources_co60.mac",
         "examples/sources/co60/co60.mac",
         ["aclust", "aeuler", "aslice", "asolid", "awalls", 
          "crmat.LINUX", "z28.a60.lvldata"],
         "output.out"),
        ("sources_ho166",  "UCGretina", "func_sources_ho166.mac",
         "examples/sources/ho166/ho166.mac",
         ["aclust", "aeuler", "aslice", "asolid", "awalls", 
          "crmat.LINUX", "z67.a166.decaydata", "z68.a166.lvldata"],
         "output.out"),
    ],
    "inbeam": [
        ("inbeam_standard",  "UCGretina",     "func_inbeam_standard.mac",
         "examples/inbeam/fit/s44_1329.mac",
         ["aclust", "aeuler", "aslice", "asolid", "awalls", 
          "crmat.LINUX", "z16.a44.lvldata"],
         "output.out"),
        ("inbeam_lh",        "UCGretina_LH",  "func_inbeam_lh.mac",
         "examples/inbeam/fitLH/s44_1329.mac",
         ["aclust", "aeuler", "aslice", "asolid", "awalls", 
          "crmat.LINUX", "z16.a44.lvldata"],
         "output.out"),
        ("inbeam_pol",       "UCGretina_Pol", "func_inbeam_pol.mac",
         "examples/inbeam/angdist/s44_1329.mac",
         ["aclust", "aeuler", "aslice", "asolid", "awalls", 
          "crmat.LINUX", "z16.a44.lvldata"],
         "output.out"),
        # cache pipeline handled separately — see run_cache_pipeline()
    ],
    "scanning": [
        ("scanning",  "UCGretina_Scan", "func_scanning.mac",
         "examples/scan/scan.mac",
         ["aclust", "aeuler", "aslice", "asolid", "awalls"],
         "output.out"),
    ],
    "background": [
        ("background", "UCGretina", "func_background.mac",
         "examples/background/background.mac",
         ["aclust", "aeuler", "aslice", "asolid", "awalls", 
          "crmat.LINUX"],
         "output.out"),
    ],
}


def run_functional(mode):
    print(f"\n=== test-{mode} ({FUNCTIONAL_EVENTS} events) ===")
    baselines = load_baselines()
    failures = 0

    for test_name, binary_name, macro_file, example_path, support_files, out_file in FUNCTIONAL_CASES[mode]:
        binary = find_binary_optional(binary_name)
        if binary is None:
            print(f"[SKIP] {test_name:<30} binary {binary_name} not found")
            continue
        workdir = setup_workdir(test_name, example_path, support_files)
        output_command = f"/Output/Filename {out_file}"
        if ".dat" in out_file:
            output_command = f"/Mode2/Filename {out_file}"
        write_base_macro(macro_file, example_path, output_command, workdir)
        wrapper = os.path.join(workdir, "run.mac")
        write_run_macro(macro_file, FUNCTIONAL_EVENTS, wrapper)
        stdout, stderr, returncode = run_sim(binary, wrapper, workdir)
        ok, msg = _check_run_criteria(test_name, stdout, stderr, returncode)
        if not ok:
            print(msg)
            failures += 1
            continue

        output_path = os.path.join(workdir, out_file)
        if not os.path.isfile(output_path):
            msg = f"[FAIL] {test_name:<30} output file not created: {output_path}"
            print(msg)
            failures += 1
            continue
        observed = count_lines(output_path)
        passed, msg = check_baseline(test_name, observed, baselines)
        print(msg)
        if not passed:
            failures += 1

    # Run cache pipeline as part of inbeam mode
    if mode == "inbeam":
        cache_failures = run_cache_pipeline(baselines)
        failures += cache_failures

    save_baselines(baselines)

    if failures:
        print(f"\n{failures} FAILED")
        sys.exit(1)
    else:
        print(f"\nAll {mode} tests passed.\n")


def run_cache_pipeline(baselines):
    """Run two-step cache pipeline: generation then playback. Returns failure count."""
    failures = 0

    # Step 1: Cache generation — requires UCGretina_LH; skip if unavailable
    gen_name = "inbeam_cache_gen"
    gen_binary = find_binary_optional("UCGretina_LH")
    if gen_binary is None:
        print(f"[SKIP] {gen_name:<30} binary UCGretina_LH not found")
        print(f"[SKIP] {'inbeam_cache_run':<30} skipped (no cache file: step 1 skipped)")
        return failures

    example_path = "examples/inbeam/cache/s44_1329_cache.mac"
    support_files = ["aclust", "aeuler", "aslice", "asolid", "awalls", 
                     "crmat.LINUX", "z16.a44.lvldata"]
    gen_workdir = setup_workdir(gen_name, example_path, support_files)
    macro_file = "func_inbeam_cache_gen.mac"
    write_base_macro(macro_file, example_path, "", gen_workdir)
    gen_wrapper = os.path.join(gen_workdir, "run.mac")
    write_run_macro(macro_file, FUNCTIONAL_EVENTS, gen_wrapper)
    stdout, stderr, returncode = run_sim(gen_binary, gen_wrapper, gen_workdir)
    ok, msg = _check_run_criteria(gen_name, stdout, stderr, returncode)
    print(msg)
    if not ok:
        return failures + 1

    cache_file = os.path.join(gen_workdir, "cache_gen.cache")
    if not os.path.isfile(cache_file) or os.path.getsize(cache_file) == 0:
        print(f"[FAIL] {gen_name:<30} cache file missing or empty: {cache_file}")
        return failures + 1

    # Step 2: Cache playback — uses UCGretina_LH (same binary as gen;
    # LH target geometry commands in the macro require the LH binary)
    run_name = "inbeam_cache_run"
    run_binary = find_binary_optional("UCGretina_LH")
    if run_binary is None:
        print(f"[SKIP] {run_name:<30} binary UCGretina_LH not found")
        return failures

    #    run_geom = os.path.join(PROJECT_ROOT, "Geometry/GretinaLH/G120C4")
    example_path = "examples/inbeam/cache/s44_1329.mac"
    support_files = ["aclust", "aeuler", "aslice", "asolid", "awalls", 
                     "crmat.LINUX", "z16.a44.lvldata"]
    run_workdir = setup_workdir(run_name, example_path, support_files)

    # Playback needs the generated cache file
    shutil.move(os.path.join(gen_workdir, "cache_gen.cache"), run_workdir)
    
    write_base_macro("func_inbeam_cache_run.mac", example_path,
                     "/Output/Filename output.out", run_workdir)
    run_wrapper = os.path.join(run_workdir, "run.mac")

    # Write wrapper: execute base macro, inject Cache/Input, then beamOn
    with open(run_wrapper, "w") as f:
        f.write("/control/execute func_inbeam_cache_run.mac\n")
        #        f.write(f"/Cache/Input {cache_file}\n")
        f.write(f"/run/beamOn {FUNCTIONAL_EVENTS}\n")

    stdout, stderr, returncode = run_sim(run_binary, run_wrapper, run_workdir)
    ok, msg = _check_run_criteria(run_name, stdout, stderr, returncode)
    if not ok:
        print(msg)
        return failures + 1

    output_path = os.path.join(run_workdir, "output.out")
    if not os.path.isfile(output_path):
        msg = f"[FAIL] {run_name:<30} output file not created: {output_path}"
        print(msg)
        return failures + 1
    observed = count_lines(output_path)
    passed, msg = check_baseline(run_name, observed, baselines)
    print(msg)
    if not passed:
        failures += 1

    return failures


# (variant_label, binary_name, macro_file, geometry_prefix)
BENCHMARK_CASES = [
    ("UCGretina",      "UCGretina",      "bench_standard.mac",
     "examples/inbeam/fit/s44_1329.mac",
     ["aclust", "aeuler", "aslice", "asolid", "awalls", 
      "crmat.LINUX", "z16.a44.lvldata"]),
    ("UCGretina_LH",   "UCGretina_LH",   "bench_lh.mac",
     "examples/inbeam/fitLH/s44_1329.mac",
     ["aclust", "aeuler", "aslice", "asolid", "awalls", 
      "crmat.LINUX", "z16.a44.lvldata"]),
    ("UCGretina_Pol",  "UCGretina_Pol",  "bench_pol.mac",
     "examples/inbeam/angdist/s44_1329.mac",
     ["aclust", "aeuler", "aslice", "asolid", "awalls", 
      "crmat.LINUX", "z16.a44.lvldata"]),
    ("UCGretina_Scan", "UCGretina_Scan", "bench_scan.mac",
     "examples/scan/scan.mac",
     ["aclust", "aeuler", "aslice", "asolid", "awalls"]),
]


def run_benchmark(n_events):
    git_hash, git_branch = get_git_info()
    cpu = get_cpu_info()
    today = datetime.date.today().isoformat()
    print(f"\n=== Benchmark ({n_events} events, commit {git_hash}, branch {git_branch}) ===")
    print(f"CPU: {cpu}")
    print(f"{'Variant':<20} {'Events/sec':>12}")
    print("-" * 34)

    rows = []

    # Standard variants
    for variant, binary_name, macro_file, example_path, support_files in BENCHMARK_CASES:
        binary = find_binary_optional(binary_name)
        if binary is None:
            print(f"  {variant:<20} {'[SKIP]':>12}")
            rows.append((today, git_hash, git_branch, cpu, variant, n_events, "SKIP"))
            continue
        workdir = setup_workdir(f"bench_{variant}", example_path, support_files)
        write_base_macro(macro_file, example_path,
                         "/Mode2/Filename output.dat", workdir)
        wrapper = os.path.join(workdir, "run.mac")
        write_run_macro(macro_file, n_events, wrapper)
        stdout, stderr, returncode = run_sim(binary, wrapper, workdir)

        eps = parse_events_per_sec(stdout)
        if eps is None or returncode != 0:
            print(f"  {variant:<20} {'ERROR':>12}")
            eps = 0
        else:
            print(f"  {variant:<20} {eps:>12.0f}")

        rows.append((today, git_hash, git_branch, cpu, variant, n_events, f"{eps:.0f}"))

    # Cache pipeline benchmark
    cache_gen_eps, cache_run_eps = _run_cache_benchmark(n_events, today,
                                                        git_hash, git_branch,
                                                        cpu, rows)

    # Cache speedup factor
    if cache_gen_eps and cache_gen_eps > 0:
        speedup = cache_run_eps / cache_gen_eps
        print(f"  {'cache_speedup':<20} {speedup:>11.1f}x")
        rows.append((today, git_hash, git_branch, cpu, "cache_speedup", "\u2014", f"{speedup:.1f}x"))
    else:
        print(f"  {'cache_speedup':<20} {'N/A':>12}")
        # cache rows already written as SKIP/N/A in _run_cache_benchmark

    append_benchmark_log(rows)
    print(f"\nResults appended to benchmark.log")


def _run_cache_benchmark(n_events, today, git_hash, git_branch, cpu, rows):
    """Run cache gen + cache run benchmarks. Appends rows. Returns (gen_eps, run_eps)."""
    # Cache generation — requires UCGretina_LH; skip gracefully if unavailable
    gen_binary = find_binary_optional("UCGretina_LH")
    if gen_binary is None:
        print(f"  {'cache_gen':<20} {'[SKIP]':>12}")
        print(f"  {'cache_run':<20} {'[SKIP]':>12}")
        rows.append((today, git_hash, git_branch, cpu, "cache_gen", n_events,
                     "SKIP"))
        rows.append((today, git_hash, git_branch, cpu, "cache_run", n_events,
                     "SKIP"))
        rows.append((today, git_hash, git_branch, cpu, "cache_speedup",
                     "\u2014", "N/A"))
        return 0, 0

    # gen_geom = os.path.join(PROJECT_ROOT, "Geometry/GretinaLH/G120C4")
    example_path = "examples/inbeam/cache/s44_1329_cache.mac"
    support_files = ["aclust", "aeuler", "aslice", "asolid", "awalls", 
                     "crmat.LINUX", "z16.a44.lvldata"]
    gen_workdir = setup_workdir("bench_cache_gen", example_path, support_files)
    write_base_macro("bench_cache_gen.mac", example_path,
                     "", gen_workdir)
    gen_wrapper = os.path.join(gen_workdir, "run.mac")
    write_run_macro("bench_cache_gen.mac", n_events, gen_wrapper)
    stdout, stderr, returncode = run_sim(gen_binary, gen_wrapper, gen_workdir)
    gen_eps = parse_events_per_sec(stdout)
    if gen_eps is None or returncode != 0:
        print(f"  {'cache_gen':<20} {'ERROR':>12}")
        return 0, 0
    print(f"  {'cache_gen':<20} {gen_eps:>12.0f}")
    rows.append((today, git_hash, git_branch, cpu, "cache_gen", n_events, f"{gen_eps:.0f}"))

    cache_file = os.path.join(gen_workdir, "cache_gen.cache")

    # Cache playback — uses UCGretina_LH (LH target geometry commands in macro)
    run_binary = find_binary_optional("UCGretina_LH")
    if run_binary is None:
        print(f"  {'cache_run':<20} {'[SKIP]':>12}")
        return gen_eps, 0

    example_path = "examples/inbeam/cache/s44_1329.mac"
    support_files = ["aclust", "aeuler", "aslice", "asolid", "awalls", 
                     "crmat.LINUX", "z16.a44.lvldata"]
    run_workdir = setup_workdir("bench_cache_run", example_path, support_files)

    # Playback needs the generated cache file
    shutil.move(os.path.join(gen_workdir, "cache_gen.cache"), run_workdir)
    
    write_base_macro("bench_cache_run.mac", example_path,
                     "/Mode2/Filename output.dat", run_workdir)
    run_wrapper = os.path.join(run_workdir, "run.mac")
    with open(run_wrapper, "w") as f:
        f.write(f"/control/execute bench_cache_run.mac\n")
        #f.write(f"/Cache/Input {cache_file}\n")
        f.write(f"/run/beamOn {n_events}\n")

    stdout, stderr, returncode = run_sim(run_binary, run_wrapper, run_workdir)
    run_eps = parse_events_per_sec(stdout)
    if run_eps is None or returncode != 0:
        print(f"  {'cache_run':<20} {'ERROR':>12}")
        return gen_eps, 0
    print(f"  {'cache_run':<20} {run_eps:>12.0f}")
    rows.append((today, git_hash, git_branch, cpu, "cache_run", n_events, f"{run_eps:.0f}"))

    return gen_eps, run_eps


def update_baselines():
    """Re-run all functional tests and reset baselines to observed values."""
    print("Resetting all baselines...")
    # Save old baselines in case we need to restore on failure
    old_content = None
    if os.path.isfile(BASELINES_FILE):
        with open(BASELINES_FILE) as f:
            old_content = f.read()
    # Wipe existing baselines so every test triggers [BASELINE SET]
    if os.path.isfile(BASELINES_FILE):
        os.remove(BASELINES_FILE)
    try:
        for mode in ["sources", "inbeam", "scanning", "background"]:
            run_functional(mode)
    except SystemExit:
        # Restore old baselines if update failed partway through
        if old_content is not None:
            with open(BASELINES_FILE, "w") as f:
                f.write(old_content)
            print("Baseline update failed — original baselines restored.")
        raise
    print("Baselines updated.")


if __name__ == "__main__":
    main()
