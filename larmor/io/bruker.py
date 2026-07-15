"""Read-only Bruker TopSpin EXPNO reader built on nmrglue.

Design guarantees:
  - never writes into the EXPNO folder; `snapshot`/`verify_untouched` make the
    guarantee checkable, and `read_expno(..., verify=True)` enforces it.
  - surfaces metadata conflicts instead of silently trusting one source
    (e.g. acqus MASR vs the operator-typed title -- they disagree in real
    NMRFAM data).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import nmrglue as ng


def snapshot(root: Path) -> dict[str, tuple[int, int]]:
    return {
        str(p): (p.stat().st_mtime_ns, p.stat().st_size)
        for p in sorted(Path(root).rglob("*")) if p.is_file()
    }


def verify_untouched(root: Path, before: dict) -> None:
    after = snapshot(root)
    if after != before:
        changed = {k for k in set(before) | set(after) if before.get(k) != after.get(k)}
        raise RuntimeError(f"read-only violation in {root}: {changed}")


@dataclass
class BrukerExperiment:
    path: str
    nucleus: str
    sfo1_MHz: float
    pulse_program: str
    td: int
    sw_Hz: float
    masr_Hz: float | None
    title: str
    fid: np.ndarray | None
    processed: np.ndarray | None          # pdata/1/1r if present
    processed_ppm: np.ndarray | None      # ppm axis for `processed`
    conflicts: list[str] = field(default_factory=list)

    @property
    def summary(self) -> str:
        lines = [
            f"EXPNO: {self.path}",
            f"nucleus: {self.nucleus}   SFO1: {self.sfo1_MHz:.4f} MHz",
            f"pulse program: {self.pulse_program}   TD: {self.td}   SW: {self.sw_Hz:.0f} Hz",
            f"MASR (acqus): {self.masr_Hz} Hz",
            f"title: {self.title.splitlines()[0] if self.title else '(none)'}",
        ]
        for c in self.conflicts:
            lines.append(f"CONFLICT: {c}")
        return "\n".join(lines)


def is_expno(path: str | Path) -> bool:
    p = Path(path)
    return p.is_dir() and (p / "acqus").exists()


def read_expno(path: str | Path, procno: int = 1, verify: bool = True) -> BrukerExperiment:
    """Read one EXPNO folder (acqus + fid/ser + pdata/<procno> when present)."""
    path = Path(path)
    before = snapshot(path) if verify else None

    dic, fid = ng.bruker.read(str(path))
    acqus = dic["acqus"]

    processed = processed_ppm = None
    pdata_dir = path / "pdata" / str(procno)
    if (pdata_dir / "procs").exists() and (pdata_dir / "1r").exists():
        dic_p, processed = ng.bruker.read_pdata(str(pdata_dir))
        procs = dic_p["procs"]
        si, sf = int(procs["SI"]), float(procs["SF"])
        offset_ppm, sw_p = float(procs["OFFSET"]), float(procs["SW_p"])
        processed_ppm = offset_ppm - np.arange(si) * (sw_p / sf / si)

    title = ""
    title_path = pdata_dir / "title"
    if title_path.exists():
        title = title_path.read_text(encoding="utf-8", errors="replace")

    exp = BrukerExperiment(
        path=str(path),
        nucleus=str(acqus.get("NUC1", "")),
        sfo1_MHz=float(acqus.get("SFO1", 0.0)),
        pulse_program=str(acqus.get("PULPROG", "")),
        td=int(acqus.get("TD", 0)),
        sw_Hz=float(acqus.get("SW_h", 0.0)),
        masr_Hz=float(acqus["MASR"]) if acqus.get("MASR") is not None else None,
        title=title,
        fid=fid,
        processed=processed,
        processed_ppm=processed_ppm,
    )
    exp.conflicts = _find_conflicts(exp)

    if verify:
        verify_untouched(path, before)
    return exp


def _find_conflicts(exp: BrukerExperiment) -> list[str]:
    conflicts = []
    # MAS rate: instrument parameter vs operator-typed title
    m = re.search(r"MASR?\s*[=:]?\s*([\d.]+)\s*kHz", exp.title, re.IGNORECASE)
    if m and exp.masr_Hz is not None:
        title_hz = float(m.group(1)) * 1000.0
        if abs(title_hz - exp.masr_Hz) > 0.02 * max(title_hz, exp.masr_Hz):
            conflicts.append(
                f"MAS rate: acqus says {exp.masr_Hz:.0f} Hz but the title says "
                f"{title_hz:.0f} Hz -- confirm which is right before fitting"
            )
    return conflicts
